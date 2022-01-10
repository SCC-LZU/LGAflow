process pasa_seq_clean {
    label 'pasa'

    input:
    path(transcripts)
    path(univec)

    output:
    path "*.fasta.cln", emit: cln
    path "*.fasta.clean", emit: clean
    script:
    """
    \$PASAHOME/bin/seqclean ${transcripts} -v ${univec}
    """

}

process pasa_concat {
    label 'smalljobs'
    input:
    path(assembly_gg)
    path(assembly_denovo)
    
    output:
    path 'assembly_conacted.fasta'
    
    script:
    """
    cat ${assembly_gg} ${assembly_denovo} > assembly_conacted.fasta
    """
}

process pasa_mysql_config {
    label 'smalljobs'
    
    input:
    path(conf_template)
    path(align_conf_template)
    val(host)
    val(db_name)
    val(username)
    val(password)
    
    output:
    path 'conf.txt',emit:pasa_conf
    path 'alignAssembly.txt', emit:pasa_alignassembly_conf

    script:
    """
    sed -e 's/MYSQLSERVER=localhost/MYSQLSERVER=${host}/1; s/MYSQL_RW_USER=xxxxxx/MYSQL_RW_USER=${username}/1; s/MYSQL_RW_PASSWORD=xxxxx/MYSQL_RW_PASSWORD=${password}/1' ${conf_template} > conf.txt
    sed 's!DATABASE=<__DATABASE__>!DATABASE=${db_name}!1' ${align_conf_template} > alignAssembly.txt
    """
    
}

process pasa_sqlite_config {
    label 'smalljobs'

    input:
    path(conf_template)
    path(align_conf_template)
    val(db_path)
    
    output:
    path 'conf.txt',emit:pasa_conf
    path 'alignAssembly.txt', emit:pasa_alignassembly_conf
    
    script:
    """
    cp ${conf_template} conf.txt
    sed 's!DATABASE=<__DATABASE__>!DATABASE=${db_path}!1' ${align_conf_template} > alignAssembly.txt
    
    """
}

process pasa_create_tdn {
    label 'pasa'
    publishDir "${params.output}/${params.pasa_dir}", mode: 'copy', pattern: "*.tdn.accs"

    input:
    path(assembly_denovo)
 
    output:
    path 'tdn.accs'

    script:
    """
    \$PASAHOME/misc_utilities/accession_extractor.pl < ${assembly_denovo} > tdn.accs
    """
}

process pasa {
    label 'pasa'
    publishDir "${params.output}/${params.pasa_dir}", mode: 'copy', pattern: "*.gff3"

    input:
    path(align_config)
    path(pasa_config)
    path(ref_genome_masked)
    path(transcript)
    path(transcript_clean)
    path(transcript_cln)
    path(tdn)

    output:
    //path '*.transdecoder.genome.gff3',emit: trans_gff
    path '*.pasa_assemblies.gff3',emit: assemblies_gff

    script:
    if ( workflow.containerEngine != null ) {
        params.pasa_additional_params = '--ALIGNERS blat,minimap2 -I 600000'
    }

    """
    rm -f ${params.pasa_sqlite_path}
    export PASACONF=${pasa_config}
    \$PASAHOME/Launch_PASA_pipeline.pl -c ${align_config} -C -R -g ${ref_genome_masked} -T -u ${transcript} -t ${transcript_clean} --CPU ${task.cpus} --TDN ${tdn} ${params.pasa_additional_params}
    unset PASACONF
    """
    
}