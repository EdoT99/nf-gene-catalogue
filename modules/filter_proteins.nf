process FILTER_PROTEINS {
    tag "evalue <= ${evalue}"
    cpus 1
    memory '4 GB'

    conda 'bioconda::seqkit=2.9.0'
    container 'quay.io/biocontainers/seqkit:2.9.0--h9ee0642_0'

    input:
    path(proteins)                              // pooled ORFs (.faa)
    path(tblouts, stageAs: 'tblout/*')          // all HMMSEARCH .tblout files
    val(evalue)                                 // keep hits with full-sequence E-value <= this

    output:
    path "hmm_hits.tsv",      emit: hits        // every ORF-profile hit passing the threshold
    path "hmm_best_hits.tsv", emit: best_hits   // best (lowest E-value) profile per ORF
    path "hmm_hits.faa",      emit: faa         // protein sequences of the annotated ORFs

    script:
    """
    # 1. All hits passing the E-value threshold
    #    tblout columns: 1 = ORF id, 3 = profile name, 5 = full-sequence E-value, 6 = score
    printf "orf_id\\tcontig_id\\tprofile\\tevalue\\tscore\\n" > hmm_hits.tsv
    cat tblout/* \\
        | grep -v '^#' \\
        | awk -v e=${evalue} 'BEGIN{OFS="\\t"} (\$5+0) <= (e+0) {
              contig = \$1; sub(/_[0-9]+\$/, "", contig)    # SAMPLE_k141_1_3 -> SAMPLE_k141_1
              print \$1, contig, \$3, \$5, \$6
          }' \\
        | sort -k1,1 -k4,4g \\
        >> hmm_hits.tsv

    # 2. Best hit per ORF (rows are sorted by ORF, then E-value)
    head -n 1 hmm_hits.tsv > hmm_best_hits.tsv
    tail -n +2 hmm_hits.tsv | awk '!seen[\$1]++' >> hmm_best_hits.tsv

    # 3. Extract the annotated ORF sequences
    tail -n +2 hmm_best_hits.tsv | cut -f1 > hit_ids.txt
    seqkit grep -f hit_ids.txt ${proteins} > hmm_hits.faa
    """
}
