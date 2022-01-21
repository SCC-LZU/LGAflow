process busco {
    label 'busco'
   
    publishDir "${params.output}/${params.busco_dir}", mode: 'copy', pattern: "*"

    input:
    path(fasta)
    path(busco_db)
    path(ausgutus_config)
    val(sepecies_name)
    val(augustus_species)
    
    output:
    path "run_${model_name}/augustus_output", emit:busco_model

    script:
    if (params.useConda) {
        run_busco = 'run_busco'
    } else {
        run_busco = 'run_BUSCO.py'
    }
    """
    export AUGUSTUS_CONIF_PATH=${ausgutus_config}
    ${run_busco} -i ${fasta} -o ${sepecies_name} -l ${busco_db} --long -m geno --species ${augustus_species} -c ${task.cpus} ${params.busco_additional_params}
    unset AUGUSTUS_CONIF_PATH
    """
}
