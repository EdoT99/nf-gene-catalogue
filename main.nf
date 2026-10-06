#!/usr/bin/env nextflow

include { PREPROCESSING  } from './subworkflows/preprocessing.nf'
include { ASSEMBLY       } from './subworkflows/assembly.nf'
include { ORF_PREDICTION } from './modules/orf_prediction.nf'

include { HMMSEARCH                } from './modules/hmmsearch.nf'        // per sample, all profiles combined
include { FILTER_HITS              } from './modules/filter_hits.nf'      // E-value + bit score + best hit per ORF
include { EXTRACT_GENES            } from './modules/extract_genes.nf'    // per sample, with its own contigs
include { POOL_ORFS as POOL_HITS   } from './modules/pool_orfs.nf'        // same prefixing logic, reused
include { POOL_ORFS as POOL_GENES  } from './modules/pool_orfs.nf'
include { DEREPLICATION            } from './modules/mmseqs.nf'
include { EXTRACT_REP_GENES        } from './modules/extract_rep_genes.nf'  // catalogue genes (nucleotide)
include { ABUNDANCE                } from './subworkflows/abundance.nf'     // Bowtie2 index + mapping + CoverM

params {
    // ---- inputs ----
    input_sample_table: Path?   = null      // CSV: sample,r1,r2         reads (assembly and/or abundance)
    input_orf_table:    Path?   = null      // CSV: sample,orfs,contigs  existing Prodigal output (skips assembly)
    reference:          Path?   = null      // existing nucleotide catalogue (rep_genes.fna):
                                            //   ONLY the abundance step runs (needs --input_sample_table)

    preprocessing:      String  = 'fastp'   // 'fastp' or 'none' (reads already trimmed)
    assembler:          String  = 'megahit' // 'megahit' or 'metaspades'

    // ---- HMM annotation ----
    hmm_db:             Path    = "${projectDir}/assets/hmm_db"
    hmm_evalue:         Float   = 1e-5      // loose reporting threshold in hmmsearch
    hmm_z:              Integer? = null     // database size for E-values; null = total ORFs of all samples
    filter_evalue:      Float   = 1e-10     // E-value cutoff in FILTER_HITS
    min_bitscore:       Float   = 30        // bit-score cutoff in FILTER_HITS (comparable across samples)

    batch:              String  = 'batch01'

    // ---- abundance ----
    skip_abundance:     Boolean = false     // stop after the dereplicated catalogue (no mapping, no CoverM)
}

workflow {
    main:

    // =========================================================================
    // 0. Parameter checks: collect every problem, stop once
    // =========================================================================
    problems = []

    if( params.reference ) {
        // abundance-only mode
        if( !params.input_sample_table ) {
            problems << "--reference needs --input_sample_table: the reads to map against the reference."
        }
        if( params.skip_abundance ) {
            problems << "--reference and --skip_abundance together would run nothing. Remove one of them."
        }
    }
    else if( !params.input_sample_table && !params.input_orf_table ) {
        problems << "Provide --input_sample_table (reads) and/or --input_orf_table (Prodigal output), or --reference."
    }

    // [ value, parameter name, required columns ]
    sheets = [
        [ params.input_sample_table, '--input_sample_table', [ 'sample', 'r1', 'r2' ] ],
        [ params.input_orf_table,    '--input_orf_table',    [ 'sample', 'orfs', 'contigs' ] ]
    ]
    sheets.each { sheet, name, required_cols ->
        if( !sheet ) {
            return
        }
        def sheet_file = file(sheet)
        if( !sheet_file.exists() ) {
            problems << "${name}: file not found: ${sheet}"
            return
        }
        def header  = sheet_file.readLines()[0].split(',').collect { it.trim() }
        def missing = required_cols - header
        if( missing ) {
            problems << "${name}: ${sheet} is missing column(s): ${missing.join(', ')}. " +
                        "Found: ${header.join(', ')}. Expected: ${required_cols.join(', ')}."
        }
    }

    if( !(params.preprocessing in ['fastp', 'none']) ) {
        problems << "--preprocessing: '${params.preprocessing}' is not valid. Must be 'fastp' or 'none'."
    }
    if( !params.reference ) {
        // only needed when the catalogue is built
        if( !(params.assembler in ['megahit', 'metaspades']) ) {
            problems << "--assembler: '${params.assembler}' is not valid. Must be 'megahit' or 'metaspades'."
        }
        if( !files("${params.hmm_db}/*.hmm") ) {
            problems << "--hmm_db: no .hmm profiles found in ${params.hmm_db}"
        }
    }

    if( problems ) {
        error "Found ${problems.size()} problem(s) with the input:\n  - " + problems.join('\n  - ')
    }

    if( params.reference && params.input_orf_table ) {
        log.warn "--reference is set: --input_orf_table is ignored (annotation and dereplication are skipped)."
    }
    if( !params.reference && params.skip_abundance && params.input_orf_table && params.input_sample_table ) {
        log.warn "--skip_abundance is set and --input_orf_table is given: the reads in " +
                 "--input_sample_table will not be used."
    }

    // channels filled only on some routes start empty
    ch_trimmed_reads  = channel.empty()     // used for assembly and mapping
    ch_trimmed_output = channel.empty()     // published only if fastp ran
    ch_fastp_reports  = channel.empty()
    ch_contigs_out    = channel.empty()
    ch_assembly_logs  = channel.empty()
    ch_contig_map     = channel.empty()

    ch_hmm_tblout     = channel.empty()
    ch_hmm_hits       = channel.empty()
    ch_hit_proteins   = channel.empty()
    ch_hit_genes      = channel.empty()
    ch_hit_coords     = channel.empty()
    ch_all_hits       = channel.empty()
    ch_pooled_hits    = channel.empty()
    ch_pooled_genes   = channel.empty()
    ch_clusters       = channel.empty()
    ch_rep_seqs       = channel.empty()
    ch_all_seqs       = channel.empty()
    ch_kept_ids       = channel.empty()
    ch_rep_genes      = channel.empty()   // the catalogue genes: built (steps 1b-3) or given (--reference)

    ch_map_logs       = channel.empty()
    ch_abundance      = channel.empty()

    // =========================================================================
    // 1a. Reads: prepared whenever given (needed for assembly AND abundance)
    // =========================================================================
    if( params.input_sample_table ) {
        ch_fastq_pairs = channel
            .fromPath(params.input_sample_table)
            .splitCsv(header: true)
            .map { row -> [ row.sample, file(row.r1, checkIfExists: true), file(row.r2, checkIfExists: true) ] }

        if( params.preprocessing == 'none' ) {
            ch_trimmed_reads = ch_fastq_pairs                  // already trimmed: use as they are
        }
        else {
            PREPROCESSING(ch_fastq_pairs, params.preprocessing)
            ch_trimmed_reads  = PREPROCESSING.out.trimmed_reads
            ch_trimmed_output = PREPROCESSING.out.trimmed_reads
            ch_fastp_reports  = PREPROCESSING.out.reports
        }
    }

    // =========================================================================
    // Reference: either given (--reference) or built (steps 1b-3)
    // =========================================================================
    if( params.reference ) {
        // ---- abundance-only mode: everything up to the catalogue is skipped ----
        ch_rep_genes = channel.value( [ params.batch, file(params.reference, checkIfExists: true) ] )
    }
    else {
        // =====================================================================
        // 1b. ONE per-sample channel: [ sample, orfs.faa, contigs.fa ]
        // =====================================================================
        if( params.input_orf_table ) {
            ch_samples = channel
                .fromPath(params.input_orf_table, checkIfExists: true)
                .splitCsv(header: true)
                .map { row -> [ row.sample, file(row.orfs, checkIfExists: true), file(row.contigs, checkIfExists: true) ] }
        }
        else {
            ASSEMBLY(ch_trimmed_reads, params.assembler)
            ch_contigs_out   = ASSEMBLY.out.contigs                 // [ sample, tool, contigs ]  (for publishing)
            ch_assembly_logs = ASSEMBLY.out.logs
            ch_contig_map    = ASSEMBLY.out.mapping                 // [ sample, tool, mapping.tsv ]

            ch_contigs = ASSEMBLY.out.contigs.map { sample, tool, f -> [ sample, f ] }   // [ sample, contigs ]

            ORF_PREDICTION(ch_contigs)                              // emits faa: [ sample, orfs.faa ]
            ch_samples = ORF_PREDICTION.out.faa.join(ch_contigs)    // [ sample, orfs.faa, contigs.fa ]
        }

        ch_orfs        = ch_samples.map { sample, faa, contigs -> [ sample, faa ] }       // [ sample, orfs.faa ]
        ch_sample_ctgs = ch_samples.map { sample, faa, contigs -> [ sample, contigs ] }   // [ sample, contigs.fa ]

        // =====================================================================
        // 2. Per-sample HMM annotation (runs in parallel, one task per sample)
        // =====================================================================
        ch_hmm_db = channel.value( file(params.hmm_db, checkIfExists: true) )

        // same database size for every sample -> comparable E-values
        ch_hmm_z = params.hmm_z
            ? channel.value(params.hmm_z)
            : ch_orfs.map { sample, faa -> faa.countFasta() }.sum()  // total ORFs over all samples

        HMMSEARCH(ch_orfs, ch_hmm_db, params.hmm_evalue, ch_hmm_z)   // tblout: [ sample, tblout ]

        FILTER_HITS(ch_orfs.join(HMMSEARCH.out.tblout), params.filter_evalue, params.min_bitscore)

        EXTRACT_GENES(FILTER_HITS.out.faa.join(ch_sample_ctgs))     // [ sample, hits.faa, contigs.fa ]

        // =====================================================================
        // 3. Gather all samples -> pool hits -> dereplicate
        // =====================================================================
        ch_gathered_hit_faa = FILTER_HITS.out.faa
            .toSortedList { a, b -> a[0] <=> b[0] }
            .map { pairs -> [ pairs.collect { it[0] }, pairs.collect { it[1] } ] }   // [ [names], [hits.faa] ]

        ch_gathered_hit_fna = EXTRACT_GENES.out.fna
            .toSortedList { a, b -> a[0] <=> b[0] }
            .map { pairs -> [ pairs.collect { it[0] }, pairs.collect { it[1] } ] }   // [ [names], [hits.fna] ]

        POOL_HITS(ch_gathered_hit_faa)       // headers become >SAMPLE_contig_12_3
        POOL_GENES(ch_gathered_hit_fna)

        DEREPLICATION(POOL_HITS.out.faa, params.batch)

        // nucleotide sequences of the representatives = nucleotide catalogue (= mapping reference)
        EXTRACT_REP_GENES(DEREPLICATION.out.kept_ids, POOL_GENES.out.faa)

        // outputs for publishing
        ch_hmm_tblout   = HMMSEARCH.out.tblout
        ch_hmm_hits     = FILTER_HITS.out.hits
        ch_hit_proteins = FILTER_HITS.out.faa
        ch_hit_genes    = EXTRACT_GENES.out.fna
        ch_hit_coords   = EXTRACT_GENES.out.coords
        ch_all_hits     = FILTER_HITS.out.hits
            .map { sample, tsv -> tsv }
            .collectFile(name: 'all_samples_hmm_hits.tsv', keepHeader: true, sort: true)
        ch_pooled_hits  = POOL_HITS.out.faa
        ch_pooled_genes = POOL_GENES.out.faa
        ch_clusters     = DEREPLICATION.out.clusters
        ch_rep_seqs     = DEREPLICATION.out.rep_seqs
        ch_all_seqs     = DEREPLICATION.out.all_seqs
        ch_kept_ids     = DEREPLICATION.out.kept_ids
        ch_rep_genes    = EXTRACT_REP_GENES.out.fna                  // [ batch, rep_genes.fna ]
    }

    // =========================================================================
    // 4. Abundance: map every sample's reads to the reference (built or given)
    // =========================================================================
    if( params.input_sample_table && !params.skip_abundance ) {
        ABUNDANCE(ch_trimmed_reads, ch_rep_genes, params.batch)
        ch_abundance = ABUNDANCE.out.table                             // [ batch, abundance.tsv ]
        ch_map_logs  = ABUNDANCE.out.logs                              // [ sample, bowtie2.log ]
    }

    publish:
    trimmed        = ch_trimmed_output
    fastp_reports  = ch_fastp_reports
    contigs        = ch_contigs_out
    assembly_logs  = ch_assembly_logs
    contig_mapping = ch_contig_map

    hmm_tblout     = ch_hmm_tblout
    hmm_hits       = ch_hmm_hits
    hit_proteins   = ch_hit_proteins
    hit_genes      = ch_hit_genes
    hit_coords     = ch_hit_coords
    all_hits       = ch_all_hits

    pooled_hits    = ch_pooled_hits
    pooled_genes   = ch_pooled_genes

    clusters       = ch_clusters
    rep_seqs       = ch_rep_seqs
    all_seqs       = ch_all_seqs
    kept_ids       = ch_kept_ids
    rep_genes      = params.reference ? channel.empty() : ch_rep_genes   // don't re-publish a given reference

    mapping_logs   = ch_map_logs
    abundance      = ch_abundance
}

output {
    // per sample
    trimmed        { path { sample, r1, r2 -> "preprocessing/${params.preprocessing}/${sample}" } ; mode 'copy' }
    fastp_reports  { path { sample, json, html -> "preprocessing/${params.preprocessing}/${sample}" } ; mode 'copy' }
    contigs        { path { sample, tool, f -> "assembly/${tool}/${sample}" } ; mode 'copy' }
    assembly_logs  { path { sample, tool, f -> "assembly/${tool}/${sample}" } ; mode 'copy' }
    contig_mapping { path { sample, tool, f -> "assembly/${tool}/${sample}" } ; mode 'copy' }

    hmm_tblout     { path { sample, f -> "hmm_annotation/${sample}" } ; mode 'copy' }
    hmm_hits       { path { sample, f -> "hmm_annotation/${sample}" } ; mode 'copy' }
    hit_proteins   { path { sample, f -> "hmm_annotation/${sample}" } ; mode 'copy' }
    hit_genes      { path { sample, f -> "hmm_annotation/${sample}" } ; mode 'copy' }
    hit_coords     { path { sample, f -> "hmm_annotation/${sample}" } ; mode 'copy' }

    // all samples
    all_hits       { path "hmm_annotation" ; mode 'copy' }
    pooled_hits    { path "pooled_hits/proteins" ; mode 'copy' }
    pooled_genes   { path "pooled_hits/genes" ; mode 'copy' }

    clusters       { path "gene_catalog/${params.batch}" ; mode 'copy' }
    rep_seqs       { path "gene_catalog/${params.batch}" ; mode 'copy' }
    all_seqs       { path "gene_catalog/${params.batch}" ; mode 'copy' }
    kept_ids       { path "gene_catalog/${params.batch}" ; mode 'copy' }
    rep_genes      { path "gene_catalog/${params.batch}" ; mode 'copy' }

    // abundance
    mapping_logs   { path "abundance/${params.batch}/mapping_logs" ; mode 'copy' }
    abundance      { path "abundance/${params.batch}" ; mode 'copy' }
}
