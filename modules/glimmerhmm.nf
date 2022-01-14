process glimmerhmm {
    label 'glimmerhmm'

    input:
    path(fasta_hardmasked)

    output:
    path "*.out.gff"

    script:
    
    """
    glimmerhmm ${fasta_hardmasked} > ${fasta_hardmasked.simpleName}.fa.out.gff
    """
}

process train_glimmerhmm {
    label 'glimmerhmm'

    input:
    path(fasta)
    path(exon)

    output:
    
    script:

    """
    trainGlimmerHMM ${fasta} ${exon} ${params.train_glimmerhmm_additional}
    """
}