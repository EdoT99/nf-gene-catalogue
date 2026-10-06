process POOL_ORFS {
    tag "${samples.size()} samples"
    cpus 1
    memory '2 GB'

    input:
    // samples and files are two lists in the SAME order; each file is staged in its own folder
    tuple val(samples), path(files, stageAs: 'input?/*')

    output:
    path "pooled.fasta", emit: faa

    script:
    """
    samples=(${samples.join(' ')})
    files=(${files.join(' ')})

    : > pooled.fasta
    for i in "\${!files[@]}"; do
        # >contig_12_3 -> >SAMPLE_contig_12_3
        awk -v s="\${samples[\$i]}" '/^>/{sub(/^>/, ">" s "_"); print; next} {print}' "\${files[\$i]}" \\
            >> pooled.fasta
    done
    """

    stub:
    """
    touch pooled.fasta
    """
}
