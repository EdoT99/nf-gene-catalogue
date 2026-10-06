// subworkflows/pre_processing.nf
include { FASTP } from '../modules/fastp.nf'

workflow PREPROCESSING {
    take:
    ch_fastq_pairs        // [ sample, R1, R2 ]
    preprocessing

    main:
    if( preprocessing == 'fastp' ) {
        FASTP(ch_fastq_pairs)
        trimmed_reads = FASTP.out.reads      // [ sample, R1_trimmed, R2_trimmed ]
        reports       = FASTP.out.reports
    }

    emit:
    ch_trimmed_reads = trimmed_reads.map  { sample, f -> [ sample, preprocessing, f ] }
    reports          = reports.map        { sample, f -> [ sample, preprocessing, f ] } 
}