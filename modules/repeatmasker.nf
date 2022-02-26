process build_database {
    label 'repeatmodeler'
    publishDir "${params.output}/${params.repeatmodeler_dir}", mode: 'copy', pattern: "*"

    input:
    path(genome_file)
    val(db_name)
    val(engine)

    output:
    path "RM_*/consensi.fa.classified", emit: model

    script:
    parallel = params.cores/4 
    """
    BuildDatabase -name ${db_name} ${genome_file}
    RepeatModeler -pa ${parallel} -database ${db_name}
    """
}

process repeatmasker {
    label 'repeatmasker'

    publishDir "${params.output}/${params.repeatmasker_dir}", mode: 'copy', pattern: "*"

    input:
    path(genome_file)
    val(species)
    val(engine)//ncbi(RMBLASTN) by default

    output:
    path "*.fa.out.gff", emit: mask_gff
    path "*.out", emit: results
    path "*.tbl", emit: tbl_file

    script: 
    parallel = params.cores/4 
    """
    RepeatMasker -pa ${parallel} -e ${engine} -species '${species}' -gff -dir . ${params.repeatmasker_additional_params} ${genome_file}
    mv ${genome_file.name}.out.gff ${genome_file.simpleName}_masker.fa.out.gff
    """
}

process repeatmodeler {
    label 'repeatmasker'

    publishDir "${params.output}/${params.repeatmodeler_dir}", mode: 'copy', pattern: "*"

    input:
    path(genome_file)
    path(model)
    val(engine)//ncbi(RMBLASTN) by default

    output:
    path "*.fa.out.gff", emit: mask_gff
    path "*.out", emit: results
    path "*.tbl", emit: tbl_file

    script:
    parallel = params.cores/4
    """
    RepeatMasker -pa ${parallel} -e ${engine} -lib ${model} -gff -dir . ${params.repeatmasker_additional_params} ${genome_file}
    mv ${genome_file.name}.out.gff ${genome_file.simpleName}_modeler.fa.out.gff
    """
}

process repeatproteinmasker {
    label 'repeat'

    publishDir "${params.output}/${params.repeatproteinmasker_dir}", mode: 'copy', pattern: "*"

    input:
    path(genome_file)
    val(engine)

    output:
    path "*.fa.out.gff", emit: mask_gff
    path "*.out", emit: results
    path "*.tbl", emit: tbl_file

    script:
    parallel = params.cores/4
    """
    RepeatProteinMask -pa ${parallel} -engine ${engine} -dir . ${repeatproteinmasker_additional_params} ${genome_file}
    """

}

