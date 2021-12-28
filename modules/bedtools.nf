process bedtools_mask {
    label 'bedtools'
    publishDir "${params.output}/${params.repeatannotation_dir}", mode: 'copy', pattern: "*"

    input:
    path(fasta)
    path(gffs)

    output:
    path '*.softmasked', emit: softmasked_file
    path '*.hardmasked', emit: hardmasked_file
    path 'mergedmask.bed', emit: repeat_annotation

    script:
    """
    cat ${gffs} | bedtools sort | bedtools merge > mergedmask.bed
    bedtools maskfasta -soft -fi ${fasta} -bed mergedmask.bed -fo ${fasta.baseName}.softmasked
    sed 'y/atcg/NNNN/' ${fasta.baseName}.softmasked > ${fasta.baseName}.hardmasked
    """

}