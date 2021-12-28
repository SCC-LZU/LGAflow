process hisat2index {
    label 'hisat2'
    input:
    path(ref_genome)

    output:
    path "${ref_genome.baseName}*.ht2"

    script:
    """
    hisat2-build -p ${task.cpus} ${ref_genome} ${ref_genome.baseName}
    """
}

process hisat2 {
    label 'hisat2'
    publishDir "${params.output}/${params.hisat2_dir}", mode: 'copy', pattern: "*.sorted.bam"

    input:
    val(sample_name)
    path(reads)
    path(reference)
    path(index)

    output:
    path "${sample_name}.sorted.bam", emit: sample_bam 
    path "${sample_name}_summary.log", emit: log

    script:
    """
    hisat2 -x ${reference.baseName} -1 ${reads[0]} -2 ${reads[1]} -p ${task.cpus} --new-summary --summary-file ${sample_name}_summary.log | samtools view -bS | samtools sort -o ${sample_name}.sorted.bam -T tmp --threads ${task.cpus}
    """

}

process index_bam {
    label 'samtools'
    label 'smallTask'
    
    publishDir "${params.output}/${params.hisat2_dir}", mode: 'copy', pattern: "*.bai"

    input:
    tuple val(sample_name), path(bam_file)

    output:
    path("${bam_file}.bai")

    script:
    """
    samtools index ${bam_file}
    """
}