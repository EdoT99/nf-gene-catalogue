process HMMSEARCH {
    tag "${proteins.baseName}"
    cpus 8
    memory '16 GB'

    conda 'bioconda::hmmer=3.4'
    container 'quay.io/biocontainers/hmmer:3.4--hdbdd923_1'

    input:
    path(proteins)          // protein FASTA to search (e.g. pooled ORFs or gene catalogue)
    path(hmm_db)            // folder containing the .hmm profiles
    val(evalue)             // e-value threshold, e.g. 1e-5

    output:
    path "hmmsearch.tblout",    emit: tblout      // one line per protein-profile hit
    path "hmmsearch.domtblout", emit: domtblout   // one line per domain hit
    path "hmmsearch.out",       emit: out         // full human-readable report

    script:
    def args = task.ext.args ?: ''
    """
    # combine all profiles into one database for a single search
    cat ${hmm_db}/*.hmm > all_profiles.hmm

    hmmsearch \\
              --cpu ${task.cpus} \\
        -E ${evalue} \\
        --tblout hmmsearch.tblout \\
        --domtblout hmmsearch.domtblout \\
        --noali \\
        -o hmmsearch.out \\
        ${args} \\
        all_profiles.hmm \\
        ${proteins}

    rm all_profiles.hmm
    """
}
