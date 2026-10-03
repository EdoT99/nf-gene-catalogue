// subworkflows/pre_processing.nf
include { FASTP } from '../modules/fastp.nf'

workflow PREPROCESSING {
    take:
    ch_fastq_pairs        // [ sample, R1, R2 ]

    main:
    FASTP(ch_fastq_pairs)

    emit:
    trimmed_reads = FASTP.out.reads      // [ sample, R1_trimmed, R2_trimmed ]
    reports       = FASTP.out.reports
}