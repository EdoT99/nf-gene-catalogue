process HMMSEARCH {
    tag "${sample}"
    cpus 4
    memory '8 GB'

    conda 'bioconda::hmmer=3.4'
    container 'quay.io/biocontainers/hmmer:3.4--hdbdd923_1'

    input:
    tuple val(sample), path(proteins)    // ONE sample's predicted proteins
    path(hmm_db)                         // folder with .hmm profiles (shared)
    val(evalue)                          // loose reporting threshold
    val(z)                               // database size for E-values, same for every sample

    output:
    tuple val(sample), path("${sample}.tblout"),    emit: tblout      // one line per protein-profile hit
    tuple val(sample), path("${sample}.domtblout"), emit: domtblout   // one line per domain hit
    tuple val(sample), path("${sample}.hmmsearch.out"), emit: out     // human-readable report

    script:
    def args = task.ext.args ?: ''
    """
    # all profiles in one database → one search per sample
    cat ${hmm_db}/*.hmm > all_profiles.hmm

    hmmsearch \\
        --cpu ${task.cpus} \\
        -E ${evalue} \\
        -Z ${z} \\
        --tblout ${sample}.tblout \\
        --domtblout ${sample}.domtblout \\
        --noali \\
        -o ${sample}.hmmsearch.out \\
        ${args} \\
        all_profiles.hmm \\
        ${proteins}

    rm all_profiles.hmm
    """

    stub:
    """
    touch ${sample}.tblout ${sample}.domtblout ${sample}.hmmsearch.out
    """
}
