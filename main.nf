#!/usr/bin/env nextflow

include { PREPROCESSING } from './subworkflows/preprocessing.nf'
include { ASSEMBLY } from './subworkflows/assembly.nf'

include { ORF_PREDICTION } from './subworkflows/orf_prediction.nf'

include { POOL_ORFS } from './modules/pool_orfs.nf'
include { POOL_ASSEMBLIES } from './modules/pool_assemblies.nf'
include { DEREPLICATION } from './modules/mmseqs.nf'
include { HMMSEARCH } from './modules/hmmsearch_per_profile.nf'
//include { FILTER_PROTEINS } from './modules/filter_proteins.nf'

//include { BOWTIE2 } from './modules/bowtie2.nf'
//include { SAMTOOLS } from './modules/samtools.nf'
//include { COVERM } from './modules/coverm.nf'

params {
    input_sample_table: Path                  // CSV with columns: sample,fastq_1,fastq_2
    batch: String = 'batch01'    // name of the batch for MMSEQ dereplication

    preprocessing = 'fastp'
    assembler = 'megahit'  // default assembler, can be overwritten via command line or yaml config
    
    input_orf_table: Path                  // CSV with columns: sample,contigs,orfs
    pooled_faa: Path? = null      // path to pooled ORFs (optional, if not provided, will pool from input CSV)
    pooled_contigs: Path? = null  // path to pooled assemblies (optional, if not provided, will pool from input CSV)
    
    hmm_db: Path = 'assets/hmm_db'                // path to HMM database for HMMSEARCH
    hmm_evalue: Float = 1e-5      // e-value threshold for HMMSEARCH, can be overwirtten via command line or yaml config
}

workflow {
    main:

    // check the samplesheet header
    required_cols = [ 'sample', 'r1', 'r2' ]

    samplesheet = file(params.input_sample_table, checkIfExists: true)
    header      = samplesheet.readLines()[0].split(',').collect { it.trim() }
    missing     = required_cols - header

    if( missing ) {
        error "Samplesheet ${params.input_sample_table} is missing column(s): ${missing.join(', ')}. " +
              "Found: ${header.join(', ')}. Expected: ${required_cols.join(', ')}."
    }
    // Get the forward and reverse reads from the samplesheet
    ch_fastq_pairs = channel
        .fromPath(params.input_sample_table)
        .splitCsv(header: true)
        .map {row ->  [ row.sample,
          file(row.r1, checkIfExists: true),
          file(row.r2, checkIfExists: true) ]
    }}

    PREPROCESSING(ch_fastq_pairs)
    ch_trimmed_reads = PREPROCESSING.out.trimmed_reads

    // Perfomr assembly with the specified assembler (megahit or metaspades)
    if( !(params.assembler in ['megahit', 'metaspades']) ) {
        error "Wrong assembler specified: ${params.assembler}. Must be either 'megahit' or 'metaspades'."
    }
    ASSEMBLY(ch_trimmed_reads, params.assembler)

    // pooling ORFs
    ch_rows = channel
        .fromPath(params.input_orf_table)
        .splitCsv(header: true)

    if( params.pooled_faa ) {
        ch_pooled_orfs = channel.value( file(params.pooled_faa, checkIfExists: true) )
    }
    else {
        ch_gathered_orfs = ch_rows
            .map { row -> [ row.sample, file(row.orfs, checkIfExists: true) ] }
            .toSortedList { a, b -> a[0] <=> b[0] }
            .map { pairs -> [ pairs.collect { it[0] }, pairs.collect { it[1] } ] }
        
        POOL_ORFS(ch_gathered_orfs)
        ch_pooled_orfs = POOL_ORFS.out.faa
    }
    // pooling assemblies
    if(params.pooled_contigs) {
        ch_pooled_assemblies = channel.value( file(params.pooled_contigs, checkIfExists: true) )
    }
    else {
        ch_gathered_assemblies = ch_rows
            .map { row -> [ row.sample, file(row.contigs, checkIfExists: true) ] }
            .toSortedList { a, b -> a[0] <=> b[0] }
            .map { pairs -> [ pairs.collect { it[0] }, pairs.collect { it[1] } ] }
        
        POOL_ASSEMBLIES(ch_gathered_assemblies)
        ch_pooled_assemblies = POOL_ASSEMBLIES.out.fna
    }
    
    // ## execute HMM-search process on the pooled ORFs
    if( params.hmm_db ) {

        hmm_db = file(params.hmm_db, checkIfExists: true)
        profiles = files("${params.hmm_db}/*.hmm")
        
        if( !profiles ) {
            error "No .hmm profiles found in: ${params.hmm_db}"
        }
        ch_proteins = ch_pooled_orfs
        ch_profiles = channel.fromList(profiles)

        HMMSEARCH(ch_profiles, ch_proteins, params.hmm_evalue)
    }

    ch_all_hits = HMMSEARCH.out.tblout
        .map { name, tbl -> tbl }
        .collectFile(name: 'all_profiles.tblout', sort: true)
    
    // ## execute FILTER_PROTEINS process to filter the pooled ORFs based on HMMSEARCH results
    //FILTER_PROTEINS(ch_pooled, ch_all_hits)
    // ## execute MMSEQ process to dereplicate sequences
    DEREPLICATION(ch_proteins, params.batch)
    

    publish:
    trimmed = PREPROCESSING.out.reads
    reports = PREPROCESSING.out.reports

    contigs = ASSEMBLY.out.contigs
    logs    = ASSEMBLY.out.logs

    collected_orfs = ch_pooled_orfs
    collected_contigs = ch_pooled_assemblies
    // HMM selection
    hmm_table = HMMSEARCH.out.tblout
    domtblout = HMMSEARCH.out.domtblout
    out = HMMSEARCH.out.out
    // filter contigs & proteins
    // dereplication_output
    clusters = DEREPLICATION.out.clusters
    rep_seqs = DEREPLICATION.out.rep_seqs
    all_seqs = DEREPLICATION.out.all_seqs
    kept_ids = DEREPLICATION.out.kept_ids
}

output {
    // preprocessing outputs
    trimmed      { path { sample, tool, f -> "preprocessing/${tool}/${sample}" } ; mode 'copy' }
    reports      { path { sample, tool, f -> "preprocessing/${tool}/${sample}" } ; mode 'copy' }
    // assebly outputs
    contigs      { path { sample, tool, f -> "assembly/${tool}/${sample}" } ; mode 'copy' }
    logs         { path { sample, tool, f -> "assembly/${tool}/${sample}" } ; mode 'copy' }

    collected_orfs { path "pooled_orfs/" ; mode 'copy' }
    collected_contigs { path "pooled_contigs/" ; mode 'copy' }

    // HMM selection
    hmm_table     { path "hmmsearch/" ; mode 'copy' }
    domtblout     { path "hmmsearch/" ; mode 'copy' }
    out           { path "hmmsearch/" ; mode 'copy' }
    // dereplication_output
    clusters       { path "gene_catalog/${params.batch}" ; mode 'copy' } 
    rep_seqs       { path "gene_catalog/${params.batch}" ; mode 'copy' } 
    all_seqs       { path "gene_catalog/${params.batch}" ; mode 'copy' } 
    kept_ids       { path "gene_catalog/${params.batch}" ; mode 'copy' }


    
}

