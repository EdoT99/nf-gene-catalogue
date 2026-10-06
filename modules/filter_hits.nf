process FILTER_HITS {
    tag "${sample}"
    cpus 1
    memory '4 GB'

    conda 'bioconda::seqkit=2.9.0'
    container 'quay.io/biocontainers/seqkit:2.9.0--h9ee0642_0'

    input:
    tuple val(sample), path(proteins), path(tblout)   // ONE sample's proteins + its hmmsearch table
    val(evalue)                                       // keep full-sequence E-value <= this
    val(min_bitscore)                                 // keep full-sequence bit score >= this

    output:
    tuple val(sample), path("${sample}.hmm_hits.tsv"), emit: hits   // best profile per ORF, passing both cutoffs
    tuple val(sample), path("${sample}.hits.faa"),     emit: faa    // protein sequences of those ORFs

    script:
    """
    # tblout columns: 1 = ORF id, 3 = profile, 5 = full-seq E-value, 6 = full-seq bit score
    # 1. both cutoffs  2. sort by ORF, then bit score (highest first)  3. keep the best profile per ORF
    printf "sample\\torf_id\\tcontig_id\\tprofile\\tevalue\\tbitscore\\n" > ${sample}.hmm_hits.tsv
    grep -v '^#' ${tblout} \\
        | awk -v e=${evalue} -v b=${min_bitscore} -v s=${sample} 'BEGIN{OFS="\\t"}
              NF >= 6 && (\$5+0) <= (e+0) && (\$6+0) >= (b+0) {
                  contig = \$1; sub(/_[0-9]+\$/, "", contig)
                  print s, \$1, contig, \$3, \$5, \$6
              }' \\
        | sort -t \$'\\t' -k2,2 -k6,6gr \\
        | awk -F '\\t' '!seen[\$2]++' \\
        >> ${sample}.hmm_hits.tsv

    # 4. protein sequences of the kept ORFs (empty file if a sample has no hits)
    tail -n +2 ${sample}.hmm_hits.tsv | cut -f2 > hit_ids.txt
    seqkit grep -f hit_ids.txt ${proteins} > ${sample}.hits.faa
    """

    stub:
    """
    printf "sample\\torf_id\\tcontig_id\\tprofile\\tevalue\\tbitscore\\n" > ${sample}.hmm_hits.tsv
    touch ${sample}.hits.faa
    """
}
