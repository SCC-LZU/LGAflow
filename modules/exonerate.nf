process exonerate {
    label 'exonerate'

    input:
    each path(fasta)
    path(uniprot_fasta)

    output:
    path '*.fa.exonerate_out.gff'

    script:
    """
    exonerate -q ${uniprot_fasta} -Q protein --softmaskquery no -t ${fasta} -T dna --softmasktarget yes --showquerygff no --showtargetgff yes --ryo 'AveragePercentIdentity: %pi\n' --showvulgar no --showalignment no ${params.exonerate_additional_params} | tee ${fasta.simpleName}.fa.exonerate_out.gff
    """
}

process split_for_exonerate{
    label 'exonerate'

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

process merge_result_exonerate {
    label 'exonerate'
    publishDir "${params.output}/${params.exonerate_dir}", pattern: "*", mode: "copy"

    input:
    path(gffs)
    val(genome_name)

    output:
    path "${genome_name}.fa.exonerate_out.gff"

    """
    cat ${gffs} > ${genome_name}.fa.exonerate_out.gff
    """
}

process convert_format_exonerate {
    label 'exonerate'
    publishDir "${params.output}/${params.exonerate_dir}", pattern: "*",  mode: "copy"

    input:
    path(convert_format_script)
    path(relocate_script)
    each path(fasta)

    output:
    path "*.fa.exonerate_out.convertd.gff"

    shell:
    '''
    perl !{convert_format_script} !{fasta} | awk \'BEGIN{i=0};{if($0~/alignment_id 1\\s/){i++};print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\t"$6"\\t"$7"\\t"$8"\\t" "ID=match_NO_" i ";AveragePercentIdentity="$19 }\' | grep -v "line" | awk -f !{relocate_script} > !{fasta.simpleName}.fa.exonerate_out.convertd.gff
    '''
}