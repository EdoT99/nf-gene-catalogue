process MEGAHIT {
    tag "${sample}"
    cpus 16
    memory '64 GB'
    time '24h'

    conda 'bioconda::megahit=1.2.9'
    container 'ghcr.io/voutcn/megahit/megahit:1.2.9'

    input:
    tuple val(sample), path(read1), path(read2)       // trimmed paired-end reads

    output:
    tuple val(sample), path("${sample}.contigs.faa"),  emit: contigs
    tuple val(sample), path("${sample}.megahit.log"), emit: log

    script:
    def args = task.ext.args ?: ''                     // e.g. '--min-contig-len 1000' or '--presets meta-sensitive'
    """
    megahit \\
        -1 ${read1} \\
        -2 ${read2} \\
        -o megahit_out \\
        --out-prefix ${sample} \\
        -t ${task.cpus} \\
        --memory ${task.memory.toBytes()} \\
        ${args}

    # common output name, whichever assembler ran
    mv megahit_out/${sample}.contigs.faa ${sample}.contigs.faa
    cp megahit_out/${sample}.log        ${sample}.megahit.log
    """

    stub:
    """
    touch ${sample}.contigs.fa ${sample}.megahit.log
    """
}
