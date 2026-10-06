include { BOWTIE2_BUILD } from '../modules/bowtie2_build.nf'   // one index, built once
include { BOWTIE2_ALIGN } from '../modules/bowtie2_align.nf'   // per sample -> sorted BAM
include { COVERM_CONTIG } from '../modules/coverm_contig.nf'   // all BAMs -> abundance table

workflow ABUNDANCE {
    take:
    ch_reads        // [ sample, R1, R2 ]           trimmed reads, one item per sample
    ch_reference    // [ batch, reference.fna ]     catalogue genes (built or given), single value
    batch           // name for the output table

    main:
    // ONE index for all samples (single-value channel, reused by every mapping task)
    BOWTIE2_BUILD(ch_reference)                                    // [ batch, bowtie2_index/ ]

    // map every sample's reads against the same index, in parallel
    BOWTIE2_ALIGN(ch_reads, BOWTIE2_BUILD.out.index)               // bam: [ sample, sample.bam ]

    // gather all BAMs -> one CoverM run -> one table (genes x samples)
    ch_all_bams = BOWTIE2_ALIGN.out.bam
        .map { sample, bam -> bam }
        .collect()

    COVERM_CONTIG(ch_all_bams, batch)

    emit:
    table = COVERM_CONTIG.out.table     // [ batch, abundance.tsv ]
    logs  = BOWTIE2_ALIGN.out.log       // [ sample, bowtie2.log ]
    bams  = BOWTIE2_ALIGN.out.bam       // [ sample, sample.bam ]   (not published by default)
    index = BOWTIE2_BUILD.out.index     // [ batch, bowtie2_index/ ] (not published by default)
}
