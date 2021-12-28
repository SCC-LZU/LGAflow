process busco {
    label 'busco'

    input:
    path(fasta)
    val(model_name)
    path(busco_db)
    val(species)

    output:
    path 'run_'+model_name+'/augustus', emit:busco_model

    script:
    """
    tar -zxf ${busco_db}
    run_busco -i ${fasta} -o ${model_name} -l ${busco_db} --long -m geno --species ${species} -c ${task.cpus} ${params.busco_additional_params}
    """
}

process download_busco_db {
    label 'basic_tools'
    label 'smalljob'
    storeDir "${params.permanentCacheDir}/databases/busco/${params.busco_db}"

    output:
    file("${params.busco_db}.tar.gz")

    script:
    """
    wget https://busco-data.ezlab.org/v4/data/lineages/${params.busco_db}.tar.gz 
    """
}