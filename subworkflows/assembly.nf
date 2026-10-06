// subworkflows/assembly.nf
include { MEGAHIT    } from '../modules/megahit.nf'
include { METASPADES } from '../modules/metaspades.nf'
include { FILTER_RENAME_CONTIGS } from '../modules/filter_rename_contigs.nf'

workflow ASSEMBLY {
    take:
    ch_reads      // [ sample, [r1, r2] ]
    assembler     // 'megahit' or 'metaspades'

    main:
    ch_raw_contigs = channel.empty()
    ch_logs    = channel.empty()

    if( assembler == 'megahit' ) {
        MEGAHIT(ch_reads)
        ch_raw_contigs = MEGAHIT.out.contigs
        ch_logs    = MEGAHIT.out.log
    }
    if( assembler == 'metaspades' ) {
        METASPADES(ch_reads)
        ch_raw_contigs = METASPADES.out.contigs
        ch_logs    = METASPADES.out.log
    }

    FILTER_RENAME_CONTIGS(ch_raw_contigs)

    emit:
    contigs = FILTER_RENAME_CONTIGS.out.contigs.map { sample, f -> [ sample, assembler, f ] }   // [ sample, tool, final_contigs.fa ]
    mapping = FILTER_RENAME_CONTIGS.out.mapping.map { sample, f -> [ sample, assembler, f ] }   // [ sample, tool, mapping.tsv ]
    logs    = ch_logs.map                        { sample, f -> [ sample, assembler, f ] }   // [ sample, tool, log ]
}
