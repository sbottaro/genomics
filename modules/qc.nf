
process multiqc{

    input:
    path(fastqc_zip)
    val(stringa)
    
    output:
    publishDir("${params.project_name}/multiqc_${stringa}/", mode:'copy')
    path("multiqc*")
    
    script:
    """
    echo ${fastqc_zip}
    multiqc .  -f -n multiqc_report_${stringa}.html
    """
    
}

process fastqc{

    cpus { 2 * task.attempt }
    memory { 2.GB * task.attempt }
    maxRetries 3

    input:
    tuple val(sample_ID),path(reads)
    val(stringa)
    
    output:
    publishDir("${params.project_name}/fastqc_${stringa}/", mode:'copy')
    path("${sample_ID}*fastqc*"),emit:fastqc
    
    script:
    """
    echo ${reads}
    fastqc -f fastq -q ${reads} -t ${task.cpus} 
    """

}
    
