process EXTRACT_GENES {
    tag "${sample}"
    cpus 2
    memory '8 GB'

    conda 'bioconda::seqkit=2.9.0'
    container 'quay.io/biocontainers/seqkit:2.9.0--h9ee0642_0'

    input:
    // ONE sample: its hit proteins, its own contigs, and (optional) an ORF->coordinates table.
    // orf_map = [] when not given: coordinates are then read from Prodigal's protein headers.
    tuple val(sample), path(hits_faa), path(contigs), path(orf_map)

    output:
    tuple val(sample), path("${sample}.hits.fna"),        emit: fna      // gene sequence of each hit, same IDs
    tuple val(sample), path("${sample}.hits_coords.tsv"), emit: coords   // contig, start, end, strand, orf_id

    script:
    def use_map = orf_map ? true : false
    """
    export GOMAXPROCS=${task.cpus}

    # IDs of the hits (first word of each header)
    grep '^>' ${hits_faa} | sed 's/^>//; s/[[:space:]].*//' > hit_ids.txt || true

    # 1. Coordinates of every hit -> contig, start, end, strand(+/-), orf_id   (1-based, inclusive)
    if ${use_map}; then
        # from an ORF table (e.g. Geomosaic orf_contig_mapping.tsv):
        #   orf_id  contig  leftmost_coord  rightmost_coord  strand(1/-1)  ...
        awk -F '\\t' 'BEGIN{OFS="\\t"}
            NR == FNR { if (FNR > 1) m[\$1] = \$2 OFS \$3 OFS \$4 OFS (\$5 == 1 ? "+" : "-"); next }
            (\$1 in m) { print m[\$1], \$1 }
        ' ${orf_map} hit_ids.txt > ${sample}.hits_coords.tsv
    else
        # from Prodigal headers:  >ORF_ID # start # end # strand(1/-1) # ...
        grep '^>' ${hits_faa} | sed 's/^>//' \\
            | awk -F ' # ' 'BEGIN{OFS="\\t"} NF >= 4 {
                  contig = \$1; sub(/_[0-9]+\$/, "", contig)
                  print contig, \$2, \$3, (\$4 == 1 ? "+" : "-"), \$1
              }' > ${sample}.hits_coords.tsv || true
    fi

    # 2. Stream the contigs and cut out each gene (reverse-complement on - strand)
    seqkit fx2tab -j ${task.cpus} -i ${contigs} \\
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

    # 3. Every hit must have coordinates AND a gene sequence
    hits=\$(wc -l < hit_ids.txt)
    coords=\$(wc -l < ${sample}.hits_coords.tsv)
    found=\$(grep -c '^>' ${sample}.hits.fna || true)
    if [ "\$hits" -ne "\$coords" ]; then
        echo "ERROR: ${sample}: \$hits hits but coordinates for only \$coords (check the ORF table / Prodigal headers)" >&2
        exit 1
    fi
    if [ "\$hits" -ne "\$found" ]; then
        echo "ERROR: ${sample}: \$hits hits but \$found genes extracted - check contig names" >&2
        exit 1
    fi
    """

    stub:
    """
    touch ${sample}.hits.fna ${sample}.hits_coords.tsv
    """
}
