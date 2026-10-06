process PRODIGAL {
    tag "${sample}"
    cpus 1                      // Prodigal is single-threaded
    memory '4 GB'

    conda 'bioconda::prodigal=2.6.3'
    container 'quay.io/biocontainers/prodigal:2.6.3--h7b50bb2_10'

    input:
    tuple val(sample), path(contigs)          // FASTA, plain or .gz

    output:
    tuple val(sample), path("${sample}.faa"), emit: faa     // proteins: >contig_12_3 # start # end # strand # ...
    tuple val(sample), path("${sample}.fna"), emit: fna     // genes
    tuple val(sample), path("${sample}.gff"), emit: gff     // coordinates

    script:
    def args = task.ext.args ?: '-p meta'
    """
    gzip -cdf ${contigs} | prodigal \\
        -a ${sample}.faa \\
        -d ${sample}.fna \\
        -o ${sample}.gff \\
        -f gff \\
        ${args}
    """

    stub:
    """
    touch ${sample}.faa ${sample}.fna ${sample}.gff
    """
}
