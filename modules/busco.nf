process busco {
    label 'busco'
   
    publishDir "${params.output}/${params.busco_dir}", mode: 'copy', pattern: "*"

    input:
    path(fasta)
    val(model_name)
    path(busco_db)
    val(species)

    output:
    path "run_${model_name}/augustus_output", emit:busco_model

    script:

    """
    busco -i ${fasta} -o ${model_name} -l ${busco_db} --long -m geno --species ${species} -c ${task.cpus} ${params.busco_additional_params}
    """
}
