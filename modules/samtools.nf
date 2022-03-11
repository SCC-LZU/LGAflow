process samtools_sam_to_bam {
    label 'samtools'

    input:
    path(sam)

    output:
    path "${sam.baseName}.bam"

    script:
    """
    samtools view -@ ${task.cpus} -b -S ${sam} -o ${sam.baseName}.bam
    """
}

