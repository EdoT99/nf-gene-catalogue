process POOL_ORFS {
    
    tag "${samples.size()} samples"

    input:
    // samples and faa_files are two lists in the SAME order.
    // All inputs are named orf_predicted.faa, so each is staged in its own folder (input1/, input2/, ...).
    // script prefix every header with the sample name: >k141_1_1 -> >SAMPLE_k141_1_1
    tuple val(samples), path(faa_files, stageAs: 'input?/*')

    output:
    path "pooled_orf_predicted.faa", emit: faa

    script:
    """
    samples=(${samples.join(' ')})
    files=(${faa_files.join(' ')})

    : > pooled_orf_predicted.faa
    for i in "\${!files[@]}"; do
        awk -v s="\${samples[\$i]}" '/^>/{sub(/^>/, ">" s "_"); print; next} {print}' "\${files[\$i]}" \\
            >> pooled_orf_predicted.faa
    done
    """
}
