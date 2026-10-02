process DEREPLICATION {

    tag "${batch}"

    conda 'bioconda::mmseqs2=18.8cc5c'
    container 'quay.io/biocontainers/mmseqs2:18.8cc5c--hd6d6fdc_0'

    input:
    path(pooled_faa)
    val(batch)
    
    output:
    tuple val(batch), path("${batch}_rep_seq.fasta"), emit: rep_seqs   // dereplicated sequences
    tuple val(batch), path("${batch}_cluster.tsv"),   emit: clusters   // representative <tab> member
    tuple val(batch), path("${batch}_all_seqs.fasta"), emit: all_seqs
    tuple val(batch), path("${batch}_kept_ids.txt"), emit: kept_ids

    script:
    // default = 100% identity over 100% length (exact dereplication); override via ext.args
    def args = task.ext.args ?: '--min-seq-id 1.0 -c 1.0 --cov-mode 1'
    """
    mmseqs easy-linclust \\
         ${pooled_faa} \\
         ${batch} \\
         tmp \\
         --threads ${task.cpus} \\
         ${args}

     grep '^>' ${batch}_clu_rep_seq.fasta | sed 's/^>//; s/ .*//' > ${batch}_kept_ids.txt
     rm -rf tmp
    """
}
