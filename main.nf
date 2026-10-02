#!/usr/bin/env nextflow

//include { FASTP } from './modules/fastp.nf'
//include { MEGAHIT } from './modules/megahit.nf'

include { POOL_ORFS } from './modules/pool_orfs.nf'
include { MMSEQ } from './modules/mmseqs.nf'



//include { BOWTIE2 } from './modules/bowtie2.nf'
//include { SAMTOOLS } from './modules/samtools.nf'
//include { COVERM } from './modules/coverm.nf'

params {
    input: Path                  // CSV with columns: sample,fastq_1,fastq_2
    batch: String = 'batch01'    // name of the batch for MMSEQ dereplication
    hmm_db: Path                 // path to HMM database for HMMSEARCH
    hmm_folder: Path             // path to folder for HMMSEARCH output
    hmm_evlue: Float = 1e-5      // e-value threshold for HMMSEARCH, can be overwirtten via command line or yaml config
}

workflow {
    main:
    //execute FASTP process to trim reads
    //FASTP(ch_reads)
    // ch_reads = channel
    //     .fromPath(params.input)
    //     .splitCsv(header: true)
    //     .map { row ->
    //         [ row.sample, file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true) ]
    //     }.view()
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
    // execute HMM-search process on the pooled ORFs
    HMMSEARCH(ch_pooled, params.hmm_db, params.hmm_folder, params.hmm_evlue)
    // execute MMSEQ process to dereplicate sequences
    MMSEQ(ch_pooled, params.batch)

    //execute  GATHER_PRODIGAL process to gather prodigal results
    POOL_ORFS(ch_gathered)
    // execute HMM-search process on the pooled ORFs
    HMMSEARCH(POOL_ORFS.out.faa, params.hmm_db, params.hmm_folder, params.hmm_evlue)
    // execute MMSEQ process to dereplicate sequences
    MMSEQ(POOL_ORFS.out.faa, params.batch)

    publish:
    //trimmed = FASTP.out.reads
    //reports = FASTP.out.reports
    collected_orfs = ch_pooled
    hmm_table = HMMSEARCH.out.hmm_table
    clusters = MMSEQ.out.clusters
    rep_seqs = MMSEQ.out.rep_seqs
    all_seqs = MMSEQ.out.all_seqs
    kept_ids = MMSEQ.out.kept_ids
}

output {
    //trimmed      { path { sample, r1, r2 -> "fastp/${sample}" } }
    //reports      { path { sample, json, html -> "fastp/${sample}" } }
    collected_orfs { path "pooled_orfs/" ; mode 'copy' }
    // HMM selection
    hmm_table     { path "hmmsearch/" ; mode 'copy' }
    // dereplication_output
    clusters       { path "gene_catalog/${params.batch}" ; mode 'copy' } 
    rep_seqs       { path "gene_catalog/${params.batch}" ; mode 'copy' } 
    all_seqs       { path "gene_catalog/${params.batch}" ; mode 'copy' } 
    kept_ids       { path "gene_catalog/${params.batch}" ; mode 'copy' }


    
}

