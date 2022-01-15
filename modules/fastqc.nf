process fastqc {
    label 'fastqc'

    input:
    path(reads) 

    output:
    path "*_fastqc.zip", emit: zip

    script:
    """
    fastqc --noextract -t ${task.cpus} ${reads}
    """
}