process run_evm{
    label 'evm'
    cpus 1

    input:
    path(command)

    output:
    path 'evm_*.log'

    shell:
    '''
    RUN_ID=$(pwd | awk -F/ '{print $NF}')
    $EVM_HOME/EvmUtils/execute_EVM_commands.pl !{command} > evm_${RUN_ID}.log
    '''
}

process evm_partition {
    label 'evm'
    publishDir "${params.output}/${params.evidencemodeler_dir}", mode: 'copy', pattern: "*"
    
    input:
    path(genome_file)
    path(de_novo_gff)
    path(protein_gff)
    path(transcript_gff)
    path(weights)
    val(segment_size)
    val(overlap_size)

    output:
    path 'partitions_list.out', emit:partition_list
    path 'commands_list.out', emit:commands_list

    script:

    gene_predictions = "--gene_predictions ${de_novo_gff}"
    transcript_alignments = "--transcript_alignments ${transcript_gff}"
    protein_alignments = "--protein_alignments ${protein_gff}"

    """
    \$EVM_HOME/EvmUtils/partition_EVM_inputs.pl --genome ${genome_file} ${gene_predictions} ${transcript_alignments} ${protein_alignments} --segmentSize ${segment_size} --overlapSize ${overlap_size} --partition_listing partitions_list.out
    \$EVM_HOME/EvmUtils/write_EVM_commands.pl --genome ${genome_file} ${gene_predictions} ${transcript_alignments} ${protein_alignments}  --weights `pwd`/${weights} --partitions partitions_list.out --output_file_name evm.out > commands_list.out
    """

}

process evm_convert_and_merge_result {
    label 'evm'
    publishDir "${params.output}/${params.evidencemodeler_dir}", mode: 'copy'

    input:
    path(genome_file)
    path(partition)
    path(results)
    path(merge_evm_gff_script)

    output:
    path 'evm.out.gff'

    script:
    """
    \$EVM_HOME/EvmUtils/recombine_EVM_partial_outputs.pl --partitions ${partition} --output_file_name evm.out
    \$EVM_HOME/EvmUtils/convert_EVM_outputs_to_GFF3.pl  --partitions ${partition} --output evm.out  --genome ${genome_file}
    perl ${merge_evm_gff_script} --partitions ${partition} --gff evm.out.gff
    """

}

process augustus_to_evm {
    label 'evm'
    publishDir "${params.output}/${params.augustus_dir}", pattern: "*", mode: "copy"

    input:
    path(augustus_out_gff)
    path(relocate_script)

    output:
    path "*.fa.augustus_out.convertd.gff"

    script:
    
    """
    \$EVM_HOME/EvmUtils/misc/augustus_GTF_to_EVM_GFF3.pl ${augustus_out_gff} | awk -f ${relocate_script} > ${augustus_out_gff.simpleName}.fa.augustus_out.convertd.gff
    """
}

process exonerate_to_evm {
    label 'evm'
    publishDir "${params.output}/${params.augustus_dir}", pattern: "*", mode: "copy"

    input:
    path(exonerate_out_gff)
    path(relocate_script)

    output:
    path "*.fa.exonerate_out.convertd.gff"

    script:
    
    """
    \$EVM_HOME/EvmUtils/misc/exonerate_gff_to_alignment_gff3.pl ${exonerate_out_gff} | awk -f ${relocate_script} > ${exonerate_out_gff.simpleName}.fa.exonerate_out.convertd.gff
    """
}

process evm_merge_input {
    label 'smalljob'

    input:
    path(de_novo_gff)
    path(transdecoder_gff)

    output:
    path 'gene_pred.gff'

    script:
    """
    cat ${de_novo_gff} ${transdecoder_gff} > gene_pred.gff
    """
}