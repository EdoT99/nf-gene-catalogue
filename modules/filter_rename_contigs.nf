process FILTER_RENAME_CONTIGS {
    tag "${sample}"
    cpus 2
    memory '4 GB'

    conda 'bioconda::seqkit=2.9.0'
    container 'quay.io/biocontainers/seqkit:2.9.0--h9ee0642_0'

    input:
    tuple val(sample), path(contigs)                   // raw contigs from MEGAHIT or metaSPAdes

    output:
    tuple val(sample), path("${sample}.final_contigs.fa"),   emit: contigs   // filtered, renamed contig_1, contig_2, ...
    tuple val(sample), path("${sample}.contig_mapping.tsv"), emit: mapping   // old_header <tab> new_header

    script:
    def args = task.ext.args ?: '-m 1000'              // seqkit seq options; default: keep contigs >= 1000 bp
    """
    # 1. Length filter (as Geomosaic: seqkit seq with user options)
    seqkit seq ${args} -j ${task.cpus} ${contigs} > filtered.fa

    # stop with a clear message if nothing passed the filter
    if [ ! -s filtered.fa ]; then
        echo "ERROR: ${sample}: the assembler produced no contigs passing the filter (${args}). Try lowering the minimum length." >&2
        exit 1
    fi

    # 2. Rename to contig_1, contig_2, ... and keep a mapping to the original headers
    printf "old_header\\tnew_header\\n" > ${sample}.contig_mapping.tsv
    awk -v map=${sample}.contig_mapping.tsv '
        /^>/ { n++; old = substr(\$0, 2); print old "\\t" "contig_" n >> map; print ">contig_" n; next }
        { print }
    ' filtered.fa > ${sample}.final_contigs.fa

    rm filtered.fa
    """

    stub:
    """
    touch ${sample}.final_contigs.fa ${sample}.contig_mapping.tsv
    """
}
