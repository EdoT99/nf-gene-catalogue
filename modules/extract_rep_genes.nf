process EXTRACT_REP_GENES {
    tag "${batch}"
    cpus 1
    memory '4 GB'

    conda 'bioconda::seqkit=2.9.0'
    container 'quay.io/biocontainers/seqkit:2.9.0--h9ee0642_0'

    input:
    tuple val(batch), path(kept_ids)       // IDs of the MMseqs2 representatives (DEREPLICATION.out.kept_ids)
    path(pooled_genes)                     // nucleotide sequences of all hits, same IDs (POOL_GENES output)

    output:
    tuple val(batch), path("${batch}_rep_genes.fna"), emit: fna   // nucleotide catalogue = mapping reference

    script:
    """
    seqkit grep -f ${kept_ids} ${pooled_genes} > ${batch}_rep_genes.fna

    # every representative protein must have its gene sequence
    expected=\$(grep -c . ${kept_ids} || true)
    found=\$(grep -c '^>' ${batch}_rep_genes.fna || true)
    if [ "\$expected" -ne "\$found" ]; then
        echo "ERROR: \$expected representatives but \$found gene sequences found - check that proteins and genes share IDs" >&2
        exit 1
    fi
    """

    stub:
    """
    touch ${batch}_rep_genes.fna
    """
}
