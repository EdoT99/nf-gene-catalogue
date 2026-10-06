process COVERM_CONTIG {
    tag "${batch}"
    cpus 8
    memory '16 GB'

    conda 'bioconda::coverm=0.7.0'
    container 'quay.io/biocontainers/coverm:0.7.0--TAG'   // TODO: replace TAG, see README / quay.io/biocontainers/coverm tags

    input:
    path(bams, stageAs: 'bams/*')    // sorted BAMs of ALL samples, named <sample>.bam
    val(batch)

    output:
    tuple val(batch), path("${batch}_abundance.tsv"), emit: table   // one row per catalogue gene, columns per sample x method

    script:
    def args = task.ext.args ?: '-m count trimmed_mean tpm --min-read-percent-identity 95 --min-read-aligned-percent 50'
    """
    coverm contig \\
        --bam-files bams/*.bam \\
        --threads ${task.cpus} \\
        ${args} \\
        --output-file ${batch}_abundance.tsv
    """

    stub:
    """
    touch ${batch}_abundance.tsv
    """
}
