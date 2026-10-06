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
include { BOWTIE2_BUILD            } from './modules/bowtie2_build.nf'      // one index, built once
include { BOWTIE2_ALIGN            } from './modules/bowtie2_align.nf'      // per sample -> sorted BAM
include { COVERM_CONTIG            } from './modules/coverm_contig.nf'      // all BAMs -> abundance table

params {
    // ---- entry points: give ONE of the two ----
    input_sample_table: Path?   = null      // CSV: sample,r1,r2               (start from reads)
    input_orf_table:    Path?   = null      // CSV: sample,orfs,contigs        (start from Prodigal output)

    preprocessing:      String  = 'fastp'
    assembler:          String  = 'megahit' // 'megahit' or 'metaspades'

    // ---- HMM annotation ----
    hmm_db:             Path    = "${projectDir}/assets/hmm_db"
    hmm_evalue:         Float   = 1e-5      // loose reporting threshold in hmmsearch
    hmm_z:              Integer? = null     // database size for E-values; null = total ORFs of all samples
    filter_evalue:      Float   = 1e-10     // E-value cutoff in FILTER_HITS
    min_bitscore:       Float   = 30        // bit-score cutoff in FILTER_HITS (comparable across samples)

    batch:              String  = 'batch01'

}

workflow {
    main:

    // =========================================================================
    // 0. Parameter checks
    // =========================================================================
    // 
    problems = []

    if( !params.input_sample_table && !params.input_orf_table ) {
        error "Provide --input_sample_table (reads) AND/or --input_orf_table (Prodigal output)."
    }
    if( params.input_sample_table ) {
        // ---- route A: from reads ----
        required_cols = [ 'sample', 'r1', 'r2' ]
        header  = file(params.input_sample_table, checkIfExists: true).readLines()[0].split(',').collect { it.trim() }
        missing = required_cols - header
        if( missing ) {
            problems << "Samplesheet ${params.input_sample_table} is missing column(s): ${missing.join(', ')}. " +
                  "Found: ${header.join(', ')}. Expected: ${required_cols.join(', ')}."
        }
    }
    if( params.input_orf_table ) {
        // ---- route B: from existing Prodigal output ----
        required_cols = [ 'sample', 'orfs', 'contigs' ]
        header  = file(params.input_orf_table, checkIfExists: true).readLines()[0].split(',').collect { it.trim() }
        missing = required_cols - header
        if( missing ) {
            problems << "Samplesheet ${params.input_orf_table} is missing column(s): ${missing.join(', ')}. " +
                  "Found: ${header.join(', ')}. Expected: ${required_cols.join(', ')}."
        }
    }
    if( !(params.assembler in ['megahit', 'metaspades']) ) {
        error "Wrong assembler specified: ${params.assembler}. Must be either 'megahit' or 'metaspades'."
    }
    profiles = files("${params.hmm_db}/*.hmm")
    if( !profiles ) {
        error "No .hmm profiles found in: ${params.hmm_db}"
    }
    if( problems ) {
        error "Found ${problems.size()} problem(s) with the input:\n  - " + problems.join('\n  - ')
    }

    ch_trimmed_reads = channel.empty()
    ch_fastp_reports = channel.empty()
    ch_contigs_out   = channel.empty()
    ch_assembly_logs = channel.empty()
    ch_contig_map    = channel.empty()

    // =========================================================================
    // 1. ONE per-sample channel: [ sample, orfs.faa, contigs.fa ]
    // =========================================================================
    if( params.input_sample_table ) {
        // ---- route A: from reads ----
        ch_fastq_pairs = channel
            .fromPath(params.input_sample_table)
            .splitCsv(header: true)
            .map { row -> [ row.sample, file(row.r1, checkIfExists: true), file(row.r2, checkIfExists: true) ] }

        PREPROCESSING(ch_fastq_pairs, params.preprocessing)
        ch_trimmed_reads = PREPROCESSING.out.trimmed_reads          // [ sample, R1, R2 ]
        ch_fastp_reports = PREPROCESSING.out.reports

        ASSEMBLY(ch_trimmed_reads, params.assembler)
        ch_contigs_out   = ASSEMBLY.out.contigs                     // [ sample, tool, contigs ]  (for publishing)
        ch_assembly_logs = ASSEMBLY.out.logs
        ch_contig_map    = ASSEMBLY.out.mapping                     // [ sample, tool, mapping.tsv ]

        ch_contigs = ASSEMBLY.out.contigs.map { sample, tool, f -> [ sample, f ] }   // [ sample, contigs ]

        ORF_PREDICTION(ch_contigs)                                  // must emit faa: [ sample, orfs.faa ]

        // join by sample name → [ sample, orfs.faa, contigs.fa ]
        ch_samples = ORF_PREDICTION.out.faa.join(ch_contigs)
    }
    else if( params.input_sample_table && params.input_orf_table) {
        // ---- route B: from existing Prodigal output ----
        ch_samples = channel
            .fromPath(params.input_orf_table, checkIfExists: true)
            .splitCsv(header: true)
            .map { row -> [ row.sample, file(row.orfs, checkIfExists: true), file(row.contigs, checkIfExists: true) ] }
    }

    // split the per-sample channel into the shapes the next steps need
    ch_orfs        = ch_samples.map { sample, faa, contigs -> [ sample, faa ] }       // [ sample, orfs.faa ]
    ch_sample_ctgs = ch_samples.map { sample, faa, contigs -> [ sample, contigs ] }   // [ sample, contigs.fa ]

    // =========================================================================
    // 2. Per-sample HMM annotation (runs in parallel, one task per sample)
    // =========================================================================
    ch_hmm_db = channel.value( file(params.hmm_db, checkIfExists: true) )

    // same database size for every sample → comparable E-values
    ch_hmm_z = params.hmm_z
        ? channel.value(params.hmm_z)
        : ch_orfs.map { sample, faa -> faa.countFasta() }.sum()      // total ORFs over all samples

    HMMSEARCH(ch_orfs, ch_hmm_db, params.hmm_evalue, ch_hmm_z)       // tblout: [ sample, tblout ]

    // pair each sample's proteins with its own hits table -> [ sample, orfs.faa, tblout ]
    ch_filter_in = ch_orfs.join(HMMSEARCH.out.tblout)

    FILTER_HITS(ch_filter_in, params.filter_evalue, params.min_bitscore)
    // FILTER_HITS.out.faa  : [ sample, hits.faa ]
    // FILTER_HITS.out.hits : [ sample, hits.tsv ]

    // pair each sample's hits with its own contigs -> [ sample, hits.faa, contigs.fa ]
    ch_extract_in = FILTER_HITS.out.faa.join(ch_sample_ctgs)

    EXTRACT_GENES(ch_extract_in)
    // EXTRACT_GENES.out.fna    : [ sample, hits.fna ]
    // EXTRACT_GENES.out.coords : [ sample, coords.tsv ]

    // =========================================================================
    // 3. Gather all samples -> pool hits -> dereplicate
    // =========================================================================
    ch_gathered_hit_faa = FILTER_HITS.out.faa
        .toSortedList { a, b -> a[0] <=> b[0] }
        .map { pairs -> [ pairs.collect { it[0] }, pairs.collect { it[1] } ] }   // [ [names], [hits.faa] ]

    ch_gathered_hit_fna = EXTRACT_GENES.out.fna
        .toSortedList { a, b -> a[0] <=> b[0] }
        .map { pairs -> [ pairs.collect { it[0] }, pairs.collect { it[1] } ] }   // [ [names], [hits.fna] ]

    POOL_HITS(ch_gathered_hit_faa)       // headers become >SAMPLE_contig_12_3
    POOL_GENES(ch_gathered_hit_fna)

    DEREPLICATION(POOL_HITS.out.faa, params.batch)

    // one combined hits table for all samples (keeps the header of the first file only)
    ch_all_hits = FILTER_HITS.out.hits
        .map { sample, tsv -> tsv }
        .collectFile(name: 'all_samples_hmm_hits.tsv', keepHeader: true, sort: true)

    // =========================================================================
    // 4. Abundance of the catalogue genes in every sample
    //    (needs reads -> only when starting from --input_sample_table)
    // =========================================================================
    ch_rep_genes   = channel.empty()
    ch_map_logs    = channel.empty()
    ch_abundance   = channel.empty()

    if( params.input_sample_table ) {

        // nucleotide sequences of the MMseqs2 representatives = the shared reference
        EXTRACT_REP_GENES(DEREPLICATION.out.kept_ids, POOL_GENES.out.faa)
        ch_rep_genes = EXTRACT_REP_GENES.out.fna                       // [ batch, rep_genes.fna ]

        // ONE index for all samples. Every input upstream is a single-value channel,
        // so the index is a value channel too and is reused by every mapping task.
        BOWTIE2_BUILD(ch_rep_genes)                                    // [ batch, bowtie2_index/ ]

        // map every sample's trimmed reads against the same index, in parallel
        BOWTIE2_ALIGN(ch_trimmed_reads, BOWTIE2_BUILD.out.index)       // bam: [ sample, sample.bam ]
        ch_map_logs = BOWTIE2_ALIGN.out.log

        // gather all BAMs -> one CoverM run -> one table (genes x samples)
        ch_all_bams = BOWTIE2_ALIGN.out.bam
            .map { sample, bam -> bam }
            .collect()

        COVERM_CONTIG(ch_all_bams, params.batch)
        ch_abundance = COVERM_CONTIG.out.table                         // [ batch, abundance.tsv ]
    }

    publish:
    trimmed        = ch_trimmed_reads
    fastp_reports  = ch_fastp_reports
    contigs        = ch_contigs_out
    assembly_logs  = ch_assembly_logs
    contig_mapping = ch_contig_map

    hmm_tblout     = HMMSEARCH.out.tblout
    hmm_hits       = FILTER_HITS.out.hits
    hit_proteins   = FILTER_HITS.out.faa
    hit_genes      = EXTRACT_GENES.out.fna
    hit_coords     = EXTRACT_GENES.out.coords
    all_hits       = ch_all_hits

    pooled_hits    = POOL_HITS.out.faa
    pooled_genes   = POOL_GENES.out.faa

    clusters       = DEREPLICATION.out.clusters
    rep_seqs       = DEREPLICATION.out.rep_seqs
    all_seqs       = DEREPLICATION.out.all_seqs
    kept_ids       = DEREPLICATION.out.kept_ids

    rep_genes      = ch_rep_genes
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
    mapping_logs   { path "abundance/mapping_logs" ; mode 'copy' }
    abundance      { path "abundance" ; mode 'copy' }
}
