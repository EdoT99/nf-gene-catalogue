#!/usr/bin/env nextflow

//include { FASTP } from './modules/fastp.nf'
//include { MEGAHIT } from './modules/megahit.nf'

include { POOL_ORFS } from './modules/pool_orfs.nf'
include { DEREPLICATION } from './modules/mmseqs.nf'
include { HMMSEARCH } from './modules/hmmsearch_per_profile.nf'
include { FILTER_PROTEINS } from './modules/filter_proteins.nf'

//include { BOWTIE2 } from './modules/bowtie2.nf'
//include { SAMTOOLS } from './modules/samtools.nf'
//include { COVERM } from './modules/coverm.nf'

params {
    input: Path                  // CSV with columns: sample,fastq_1,fastq_2
    batch: String = 'batch01'    // name of the batch for MMSEQ dereplication
    pooled_faa: Path = null      // path to pooled ORFs (optional, if not provided, will pool from input CSV)
    hmm_db: Path = 'assets/hmm_db'                // path to HMM database for HMMSEARCH
    hmm_evalue: Float = 1e-5      // e-value threshold for HMMSEARCH, can be overwirtten via command line or yaml config
}

workflow {
    main:
   
    // ## execute POOL_ORFS process to pool ORFs from multiple samples into a single FASTA file
    if( params.pooled_faa ) {
        // Already pooled orfs, so start from MMSEQ with the file you provide
        ch_pooled = channel.value( file(params.pooled_faa, checkIfExists: true) )
    }
    else {
        // gather prodigal results from the CSV file
        ch_gathered = channel
            .fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> [ row.sample, file(row.orfs, checkIfExists: true) ] }
            .toSortedList { a, b -> a[0] <=> b[0] }
            .map { pairs -> [ pairs.collect { it[0] }, pairs.collect { it[1] } ] }

        POOL_ORFS(ch_gathered)
        ch_pooled = POOL_ORFS.out.faa
    }

    // ## execute HMM-search process on the pooled ORFs
    if( params.hmm_db ) {

        hmm_db = file(params.hmm_db, checkIfExists: true)
        profiles = files("${params.hmm_db}/*.hmm")
        
        if( !profiles ) {
            error "No .hmm profiles found in: ${params.hmm_db}"
        }
        ch_proteins = ch_pooled
        ch_profiles = channel.fromList(profiles)

        HMMSEARCH(ch_profiles, ch_proteins, params.hmm_evalue)
    }

    ch_all_hits = HMMSEARCH.out.tblout
        .map { name, tbl -> tbl }
        .collectFile(name: 'all_profiles.tblout', sort: true)
    
    // ## execute FILTER_PROTEINS process to filter the pooled ORFs based on HMMSEARCH results
    //FILTER_PROTEINS(ch_pooled, ch_all_hits)
    // ## execute MMSEQ process to dereplicate sequences
    DEREPLICATION(ch_pooled, params.batch)
    

    publish:
    //trimmed = FASTP.out.reads
    //reports = FASTP.out.reports
    collected_orfs = ch_pooled
    // HMM selection
    hmm_table = HMMSEARCH.out.tblout
    domtblout = HMMSEARCH.out.domtblout
    out = HMMSEARCH.out.out
    // dereplication_output
    clusters = DEREPLICATION.out.clusters
    rep_seqs = DEREPLICATION.out.rep_seqs
    all_seqs = DEREPLICATION.out.all_seqs
    kept_ids = DEREPLICATION.out.kept_ids
}

output {
    //trimmed      { path { sample, r1, r2 -> "fastp/${sample}" } }
    //reports      { path { sample, json, html -> "fastp/${sample}" } }
    collected_orfs { path "pooled_orfs/" ; mode 'copy' }
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

