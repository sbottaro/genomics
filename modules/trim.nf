
process trim_illumina{

    cpus 8
    memory { 10.GB * task.attempt }
    maxRetries 3

    input:
    tuple val(sample_ID),path(reads)
    path(trimmomatic_file)

    output:
    publishDir("${params.project_name}/fastq_trim/", mode:'copy')
    tuple val(sample_ID), path("${sample_ID}_{1,2}P.fastq.gz"), emit:trimmed_reads
    path("${sample_ID}.summary"), emit: stats
    
    script:
    """
    # check if trimmomatic file exists. 
    # Trimmomatic does not complain if absent 
    if [ -f "${trimmomatic_file}" ]; then
    echo "${trimmomatic_file} exists."
    else 
	echo "${trimmomatic_file} does not exist."
    exit 1
    fi
    trimmomatic PE -threads ${task.cpus} ${reads} \
	${sample_ID}_1P.fastq.gz ${sample_ID}_1U.fastq.gz \
	${sample_ID}_2P.fastq.gz ${sample_ID}_2U.fastq.gz \
	ILLUMINACLIP:${trimmomatic_file}${params.trimmomatic_options}\
	-summary ${sample_ID}.summary
    """

    
}
