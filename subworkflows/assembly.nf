// subworkflows/assembly.nf
include { MEGAHIT    } from '../modules/megahit.nf'
include { METASPADES } from '../modules/metaspades.nf'

workflow ASSEMBLY {
    take:
    ch_reads      // [ sample, [r1, r2] ]
    assembler     // 'megahit' or 'metaspades'

    main:
    ch_contigs = channel.empty()
    ch_logs    = channel.empty()

    if( assembler == 'megahit' ) {
        MEGAHIT(ch_reads)
        ch_contigs = MEGAHIT.out.contigs
        ch_logs    = MEGAHIT.out.log
    }
    if( assembler == 'metaspades' ) {
        METASPADES(ch_reads)
        ch_contigs = METASPADES.out.contigs
        ch_logs    = METASPADES.out.log
    }

    emit:
    ccontigs = ch_contigs.map { sample, f -> [ sample, assembler, f ] }   // [ sample, tool, contigs.fa ]
    logs    = ch_logs.map    { sample, f -> [ sample, assembler, f ] }   // [ sample, tool, log ]
}
