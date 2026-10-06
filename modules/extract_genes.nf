process EXTRACT_GENES {
    tag "${sample}"
    cpus 1
    memory '8 GB'

    conda 'bioconda::seqkit=2.9.0'
    container 'quay.io/biocontainers/seqkit:2.9.0--h9ee0642_0'

    input:
    tuple val(sample), path(hits_faa), path(contigs)   // ONE sample's hit proteins + its own contigs

    output:
    tuple val(sample), path("${sample}.hits.fna"),        emit: fna      // gene sequence of each hit, same IDs
    tuple val(sample), path("${sample}.hits_coords.tsv"), emit: coords   // contig, start, end, strand, orf_id

    script:
    """
    # 1. Coordinates from Prodigal headers:  >ORF_ID # start # end # strand(1/-1) # ...
    grep '^>' ${hits_faa} | sed 's/^>//' \\
        | awk -F ' # ' 'BEGIN{OFS="\\t"} {
              contig = \$1; sub(/_[0-9]+\$/, "", contig)
              print contig, \$2, \$3, (\$4 == 1 ? "+" : "-"), \$1
          }' > ${sample}.hits_coords.tsv || true

    # 2. Stream the contigs and cut out each gene (reverse-complement on - strand)
    seqkit fx2tab -i ${contigs} \\
        | awk -F '\\t' '
            function revcomp(s,   i, c, r) {
                r = ""
                for (i = length(s); i > 0; i--) { c = substr(s, i, 1); r = r ((c in comp) ? comp[c] : "N") }
                return r
            }
            BEGIN {
                split("A C G T N a c g t n", f, " "); split("T G C A N t g c a n", t, " ")
                for (i in f) comp[f[i]] = t[i]
            }
            NR == FNR { n[\$1]++; region[\$1, n[\$1]] = \$0; next }
            (\$1 in n) {
                for (k = 1; k <= n[\$1]; k++) {
                    split(region[\$1, k], r, "\\t")
                    seq = substr(\$2, r[2], r[3] - r[2] + 1)
                    if (r[4] == "-") seq = revcomp(seq)
                    print ">" r[5] "\\n" seq
                }
            }' ${sample}.hits_coords.tsv - > ${sample}.hits.fna

    # 3. Every hit must have been found in this sample's contigs
    expected=\$(wc -l < ${sample}.hits_coords.tsv)
    found=\$(grep -c '^>' ${sample}.hits.fna || true)
    if [ "\$expected" -ne "\$found" ]; then
        echo "ERROR: ${sample}: \$expected hits but \$found genes extracted - check contig names" >&2
        exit 1
    fi
    """

    stub:
    """
    touch ${sample}.hits.fna ${sample}.hits_coords.tsv
    """
}
