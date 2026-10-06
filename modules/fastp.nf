process FASTP {
    tag "${sample}"
    cpus 4
    memory '8 GB'

    conda 'bioconda::fastp=0.23.4'
    container 'quay.io/biocontainers/fastp:0.23.4--hadf994f_2'

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    tuple val(sample), path("${sample}_R1.trimmed.fastq.gz"), path("${sample}_R2.trimmed.fastq.gz"), emit: reads
    tuple val(sample), path("${sample}.fastp.json"), path("${sample}.fastp.html"),                   emit: reports

    script:
    def args = task.ext.args ?: ''
    """
    fastp \\
        --in1 ${read1} \\
        --in2 ${read2} \\
        --out1 ${sample}_R1.trimmed.fastq.gz \\
        --out2 ${sample}_R2.trimmed.fastq.gz \\
        --detect_adapter_for_pe \\
        --thread ${task.cpus} \\
        --json ${sample}.fastp.json \\
        --html ${sample}.fastp.html \\
        ${args}
    """

    stub:
    """
    echo "" | gzip > ${sample}_R1.trimmed.fastq.gz
    echo "" | gzip > ${sample}_R2.trimmed.fastq.gz
    touch ${sample}.fastp.json ${sample}.fastp.html
    """
}
