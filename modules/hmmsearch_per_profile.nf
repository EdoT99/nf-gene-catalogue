process HMMSEARCH {
    tag "${profile.baseName}"
    cpus 2
    memory '4 GB'

    conda 'bioconda::hmmer=3.4'
    container 'quay.io/biocontainers/hmmer:3.4--hdbdd923_1'

    input:
    //path(profile)           // ONE .hmm profile per task
    //path(proteins)          // protein FASTA (same file for every task)
    //val(evalue)             // e-value threshold, e.g. 1e-5
    tuple val(sample), path(sample_proteins), path(profile) val(evalue)       // trimmed paired-end reads

    output:
    tuple val(profile.baseName), path("${profile.baseName}.tblout"),    emit: tblout
    tuple val(profile.baseName), path("${profile.baseName}.domtblout"), emit: domtblout
    tuple val(profile.baseName), path("${profile.baseName}.out"),       emit: out

    script:
    def args   = task.ext.args ?: ''
    def prefix = profile.baseName
    """
    hmmsearch \\
        --cpu ${task.cpus} \\
        -E ${evalue} \\
        --tblout ${prefix}.tblout \\
        --domtblout ${prefix}.domtblout \\
        --noali \\
        -o ${prefix}.out \\
        ${args} \\
        ${profile} \\
        ${sample_proteins}
    """
}
