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
    if ( workflow.containerEngine == 'docker' ) {
        """
        bash -c 'source /root/.bashrc && source activate busco && run_busco -i ${fasta} -o ${model_name} -l ${busco_db} --long -m geno --species ${species} -c ${task.cpus} ${params.busco_additional_params}'
        """
    } else {
        """
        run_busco -i ${fasta} -o ${model_name} -l ${busco_db} --long -m geno --species ${species} -c ${task.cpus} ${params.busco_additional_params}
        """
    }
}
