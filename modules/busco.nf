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
    path "BUSCO_${sepecies_name}", emit:busco_model

    script:
    if (params.useConda) {
        run_busco = 'run_busco'
    } else {
        run_busco = 'run_BUSCO.py'
    }
    """
    export AUGUSTUS_CONFIG_PATH=${ausgutus_config} && \
    ${run_busco} -i ${fasta} -o ${sepecies_name} -l ${busco_db} --long -m geno --species ${augustus_species} -c ${task.cpus} ${params.busco_additional_params} && \
    unset AUGUSTUS_CONIF_PATH && \
    mkdir BUSCO_${sepecies_name} && \
    cp run_${sepecies_name}/augustus_output/retraining_parameters/* ./BUSCO_${sepecies_name}/
    """
}
