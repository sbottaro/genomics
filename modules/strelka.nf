process strelka{

    cpus 2
    memory '20 GB'


    input:
    tuple val(sample_name), path(tumor_bam), path(normal_bam)
    path(interval_list)
    path(resource_dir)

    publishDir "${params.project_name}/strelka/", mode: "copy"
    output:
    //path("StrelkaSomaticWorkflow/results/*")
    path("${sample_name}.strelka.{snvs,indels}.vcf.gz*"), emit: vcfs
    //path("${sample_name}.somatic.indels.snvs.gz*"), emit: indels
    
    script:
    """
    flags=""
    if [[ "${params.is_exome_strelka}" == "yes" ]]
    then
    flags=\${flags}" --exome"
    fi 

    picard IntervalListToBed -I ${params.interval_list} -O interval_list.bed
    sort -V -k1,1 -k2,2 interval_list.bed > sorted.bed
    bgzip sorted.bed 
    tabix sorted.bed.gz  

    tumor_bamfile=`ls ${tumor_bam} | grep -v .bai`
    normal_bamfile=`ls ${normal_bam} | grep -v .bai`

    /opt/strelka/bin/configureStrelkaSomaticWorkflow.py \
	--normalBam \${normal_bamfile} --tumorBam \${tumor_bamfile} \
	--referenceFasta ${params.reference_genome} \
	\${flags} --callRegions sorted.bed.gz 


    python2.7 StrelkaSomaticWorkflow/runWorkflow.py -m local -j ${task.cpus} -g 20


    mv StrelkaSomaticWorkflow/results/variants/somatic.indels.vcf.gz ${sample_name}.strelka.indels.vcf.gz 
    mv StrelkaSomaticWorkflow/results/variants/somatic.snvs.vcf.gz ${sample_name}.strelka.snvs.vcf.gz 
    mv StrelkaSomaticWorkflow/results/variants/somatic.indels.vcf.gz.tbi ${sample_name}.strelka.indels.vcf.gz.tbi 
    mv StrelkaSomaticWorkflow/results/variants/somatic.snvs.vcf.gz.tbi ${sample_name}.strelka.snvs.vcf.gz.tbi 

    """
    //	#--indelCandidates ${indel_candidates_vcf} \
}
