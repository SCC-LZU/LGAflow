process fastp {
    label 'fastp'
    publishDir "${params.output}/${params.fastp_dir}", mode: 'copy', pattern: "*.trimmed.fastq.gz" 

    input:
    val(name)
    path(reads)

    output:
    path("${name}*.trimmed.fastq.gz"), emit: sample_trimmed
    path "${name}_fastp.json", emit: json_report

    script:
    """
    fastp -i ${reads[0]} -I ${reads[1]} -o ${name}.R1.trimmed.fastq.gz -O ${name}.R2.trimmed.fastq.gz --thread ${task.cpus} --json ${name}_fastp.json ${params.fastp_additional_params}
    """
}
