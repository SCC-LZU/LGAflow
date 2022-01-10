process  cdhit{
    label cdhit
    
    publishDir "${params.output}/${params.cdhit_dir}", mode: 'copy', pattern: "merged_prot.fa"

    input:
    path(proteins)

    output:
    path 'merged_prot.fa'

    script:
    MEM=task.memory * 1000
    """
    cat ${proteins} > all.fa
    cd-hit -i all.fa -o merged_prot.fa -T ${task.cpus} -M ${MEM} ${params.cdhit_additional_params} 
    """
}