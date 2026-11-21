process BOWTIE2_BUILD {
    tag "${id}"
    cpus 8
    memory '16 GB'
    time '12h'

    conda 'bioconda::bowtie2=2.5.4'
    container 'quay.io/biocontainers/bowtie2:2.5.4--he96a11b_6'

    input:
    tuple val(id), path(fasta)              // [ batch, rep_genes.fna ]: plain or .gz

    output:
    tuple val(id), path("bowtie2_index"), emit: index    // folder with ${id}.*.bt2 (or .bt2l for large references)

    script:
    def args   = task.ext.args ?: ''        // e.g. '--seed 42'
    // decide in Nextflow, not with 'gzip -f' (minimal gzip in containers does not pass plain files through)
    def unpack = fasta.name.endsWith('.gz')
        ? "gzip -cd ${fasta} > reference.fa"     // compressed: decompress
        : "ln -s ${fasta} reference.fa"          // plain: link, no copy
    """
    mkdir bowtie2_index

    ${unpack}

    bowtie2-build \\
        --threads ${task.cpus} \\
        ${args} \\
        reference.fa \\
        bowtie2_index/${id}

    rm reference.fa
    """

    stub:
    """
    mkdir bowtie2_index
    touch bowtie2_index/${id}.1.bt2 bowtie2_index/${id}.2.bt2 bowtie2_index/${id}.3.bt2 \\
          bowtie2_index/${id}.4.bt2 bowtie2_index/${id}.rev.1.bt2 bowtie2_index/${id}.rev.2.bt2
    """
}
