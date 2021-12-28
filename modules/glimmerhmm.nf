process glimmerhmm {
    label 'glimmer'

    input:
    path(fasta_hardmasked)

    output:
    path "*.out.gff"

    script:
    
    """
    glimmerhmm ${fasta_hardmasked} > ${fasta_hardmasked.simpleName}.fa.out.gff
    """
}

process split_for_glimmerhmm {

}

process merge_result_glimmerhmm {

}

process convert_format_glimmerhmm {

}