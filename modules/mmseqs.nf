process DEREPLICATION {
    tag "${batch}"
    cpus 8
    memory '16 GB'

    conda 'bioconda::mmseqs2=18.8cc5c'
    container 'quay.io/biocontainers/mmseqs2:18.8cc5c--hd6d6fdc_0'

    input:
    path(pooled_faa)      // pooled hit proteins
    val(batch)

    output:
    tuple val(batch), path("${batch}_rep_seq.fasta"),  emit: rep_seqs
    tuple val(batch), path("${batch}_cluster.tsv"),    emit: clusters
    tuple val(batch), path("${batch}_all_seqs.fasta"), emit: all_seqs
    tuple val(batch), path("${batch}_kept_ids.txt"),   emit: kept_ids

    script:
    def args = task.ext.args ?: '--min-seq-id 0.95 -c 0.9 --cov-mode 1'
    """
    mmseqs easy-linclust \\
        ${pooled_faa} \\
        ${batch} \\
        tmp \\
        --threads ${task.cpus} \\
        ${args}

    grep '^>' ${batch}_rep_seq.fasta | sed 's/^>//; s/ .*//' > ${batch}_kept_ids.txt
    rm -rf tmp
    """

    stub:
    """
    touch ${batch}_rep_seq.fasta ${batch}_cluster.tsv ${batch}_all_seqs.fasta ${batch}_kept_ids.txt
    """
}
