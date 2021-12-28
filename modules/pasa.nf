process pasa_seq_clean {
    label 'pasa'

    input:
    path(transcripts)
    //path(univec)

    output:
    path "*.fasta.cln", emit: cln
    path "*.fasta.clean", emit: clean
    script:
    /*
    if (univec) {
        vector_sequences_option = "-v ${univec}"
    } else {
        vector_sequences_option = ''
    }
    */
    """
    \$PASAHOME/bin/seqclean ${transcripts}
    """

}

process pasa_concat {
    label 'pasa'
    
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
    path(config)
    path(ref_genome_masked)
    path(transcript)
    path(transcript_clean)
    path(transcript_cln)
    path(tdn)

    output:
    path '*.transdecoder.genome.gff3',emit: trans_gff
    path '*.pasa_assemblies.gff3',emit: assemble_gff

    script:
    
    """
    rm -f /tmp/mydb.sqlite
    \$PASAHOME/Launch_PASA_pipeline.pl -c ${config} -C -R -g ${ref_genome_masked} -T -u ${transcript} -t ${transcript_clean} --CPU ${task.cpus} --TRANSDECODER ${params.pasa_additional_options}
    """
    
}