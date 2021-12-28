process seqkit_to_uppercase {
    label 'seqkit'
    publishDir "${params.output}/${params.seqkit_dir}", mode: 'copy', pattern: "*"

    input:
    path(fasta_lowercase)

    output:
    path "*_uppercase.fa"
    
    """
    seqkit seq -u ${fasta_lowercase} ${params.seqkit_additional_params} > ${fasta_lowercase.simpleName}_uppercase.fa
    """
}

process seqkit_to_singleline {
    label 'seqkit'
    publishDir "${params.output}/${params.seqkit_dir}", mode: 'copy', pattern: "*"

    input:
    path(fasta_multiline)
    
    output:
    path "*_singleline.fa"
    
    """
    seqkit seq ${fasta_multiline} -w 0 ${params.seqkit_additional_params} > ${fasta_multiline.simpleName}_singleline.fa
    """

}

process seqkit_length_filter {
    label 'seqkit'
    publishDir "${params.output}/${params.seqkit_dir}", mode: 'copy', pattern: "*"

    input:
    path(fasta)
    
    output:
    path "*_filted.fa"

    """
    seqkit seq ${fasta} -m 500 -g ${params.seqkit_additional_params} > ${fasta.simpleName}_filted.fa
    """

}

process seqkit_sliding_trimming {
    label 'seqkit'
    publishDir "${params.output}/${params.seqkit_dir}", mode: 'copy', pattern: "*"

    input:
    path(fasta)

    output:
    path '*_trimmed.fa'

    """
    seqkit sliding -g -s 900000 -W 1000000 ${params.seqkit_additional_params} ${fasta} | seqkit seq -w 0 ${params.seqkit_additional_params} > ${fasta.simpleName}_trimmed.fa
    """
}
