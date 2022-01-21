
process augustus { 
    label 'augustus'

    input:
    path(fasta)
    path(config)
    val(species)

    output:
    path "*.gff"

    script:
    """
    export AUGUSTUS_CONIF_PATH=${config}
    augustus --softmasking=1 --species=${species} --UTR=off ${params.augustus_additional_params} ${fasta} > ${fasta.name}.augustus.out.gff
    unset AUGUSTUS_CONIF_PATH
    """

}

process augustus_partition {
    label 'augustus'

    input:
    path(fasta)
    path(split_genome_script)
    //val(slice_length)

    output:
    path "Split-*.fa"

    script:
    """
    perl ${split_genome_script} 3e6 ${fasta}
    """
}


process copy_augustus_model {
    label 'augustus'
    input:
    path(augustus_model)
    path(augustus_config_path)

    script:
    """
    cp -r -f -u ${augustus_model} ${augustus_config_path}/species
    """

}

process create_augustus_config {
    label 'augustus'

    input:
    path(config_tarball)

    output:
    path './config'

    script:
    """
    tar -xzvf ${config_tarball}
    """
}

process merge_result_augustus {
    label 'augustus'
    publishDir "${params.output}/${params.augustus_dir}", pattern: "*", mode: "copy"

    input:
    path(gffs)
    val(genome_name)

    output:
    path "${genome_name}.fa.augustus_out.gff"

    """
    cat ${gffs} > ${genome_name}.fa.augustus_out.gff
    """
}

process convert_format_augustus {
    label 'augustus'
    publishDir "${params.output}/${params.augustus_dir}", pattern: "*", mode: "copy"

    input:
    path(convert_format_script)
    each path(fasta)

    output:
    path "*.fa.augustus_out.convertd.gff"

    script:
    """
    perl ${convert_format_script} ${fasta} ${fasta.simpleName}.fa.augustus_out.convertd.gff
    """

}