process BOWTIE2_ALIGN {
    tag "${sample}"
    cpus 8
    memory '16 GB'
    time '12h'

    // bowtie2 + samtools in one image (same mulled container used by nf-core bowtie2/align)
    conda 'bioconda::bowtie2=2.4.4 bioconda::samtools=1.15.1 conda-forge::pigz=2.6'
    container 'quay.io/biocontainers/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:1744f68fe955578c63054b55309e05b41c37a80d-0'

    input:
    tuple val(sample), path(read1), path(read2)     // trimmed reads of ONE sample
    tuple val(ref_id), path(index)                  // shared index: [ batch, bowtie2_index/ ]

    output:
    tuple val(sample), path("${sample}.bam"),         emit: bam   // sorted by coordinate (required by CoverM)
    tuple val(sample), path("${sample}.bowtie2.log"), emit: log   // alignment rate summary

    script:
    def args = task.ext.args ?: ''                  // e.g. '--very-sensitive'
    """
    bowtie2 \\
        -x ${index}/${ref_id} \\
        -1 ${read1} \\
        -2 ${read2} \\
        --threads ${task.cpus} \\
        ${args} \\
        2> ${sample}.bowtie2.log \\
        | samtools sort -@ ${task.cpus} -o ${sample}.bam -
    """

    stub:
    """
    touch ${sample}.bam ${sample}.bowtie2.log
    """
}
