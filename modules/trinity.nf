process trinity_de_novo_assembly {
    label 'trinity'
    publishDir "${params.output}/${params.trinity_dir}/denovo", mode: 'copy', pattern: "trinity_out_dir"

    input:
    path(reads)

    output:
    path "trinity.fasta"

    script:
    """
    MEM=\$(echo ${task.memory} | awk '{print \$1}')
    Trinity --seqType fq --left ${reads[0]} --right ${reads[1]} --CPU ${task.cpus} --max_memory \${MEM}G
    cp trinity_out_dir.Trinity.fasta trinity.fasta
    """

}

process tririty_genome_guided_assembly {
    label 'trinity'
    publishDir "${params.output}/${params.trinity_dir}/genome_guided", mode: 'copy', pattern: "trinity_out_dir"

    input:
    path(genome_bam)
    val(max_intron)

    output:
    path "trinity-gg.fasta"

    script:
    """
    MEM=\$(echo ${task.memory} | awk '{print \$1}')
    Trinity --genome_guided_bam ${genome_bam} --genome_guided_max_intron ${max_intron} --CPU ${task.cpus} --max_memory \${MEM}G
    cp trinity_out_dir.Trinity-GG.fasta trinity-gg.fasta
    """

}