process METASPADES {
    tag "${sample}"
    cpus 16
    memory '128 GB'                                    // metaSPAdes usually needs much more memory than MEGAHIT
    time '48h'

    conda 'bioconda::spades=4.0.0'
    container 'staphb/spades:4.0.0'

    input:
    tuple val(sample), path(read1), path(read2)       // trimmed paired-end reads (metaSPAdes requires paired-end)

    output:
    tuple val(sample), path("${sample}.contigs.fa"),     emit: contigs
    tuple val(sample), path("${sample}.metaspades.log"), emit: log

    script:
    def args = task.ext.args ?: ''                     // e.g. '-k 21,33,55,77'
    """
    # same as Geomosaic: metagenomic mode, assembly only (reads are already trimmed)
    spades.py \\
        --meta \\
        --only-assembler \\
        -1 ${read1} \\
        -2 ${read2} \\
        -o spades_out \\
        -t ${task.cpus} \\
        -m ${task.memory.toGiga()} \\
        ${args}

    # common output name, whichever assembler ran
    mv spades_out/contigs.fasta ${sample}.contigs.fa
    cp spades_out/spades.log    ${sample}.metaspades.log
    """

    stub:
    """
    touch ${sample}.contigs.fa ${sample}.metaspades.log
    """
}
