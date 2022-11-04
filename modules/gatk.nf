process split_align_mark{

    cpus 4
    memory '40 GB'

    input:
    tuple val(sample_ID),path(reads)
    path(resource_dir)
    
    output:
    publishDir "${params.project_name}/dedup_bam/",mode: 'copy'
    tuple val(sample_ID),path("${sample_ID}_merged_dedup_sorted_tagged.ba*"), emit: bam_bai

    
    script:
    """
    # This script is UGLY. But avoids lots of R/W between scratch on the local node and 
    # /home. Therefore the splitting, BWAMEM and mark duplicates is done in one go. 

    ############################################3
    # Split reads according to read groups

    zcat  < ${reads[0]} | awk 'BEGIN {FS = ":| "} {
    fc=\$3; lane=\$4; 
    name="${sample_ID}."fc".L"lane".R1.fastq"; 
    print > name;  
    for (i = 1; i <= 3; i++) {getline; print > name }}' & 


    zcat  < ${reads[1]} | awk 'BEGIN {FS = ":| "} {
    fc=\$3; lane=\$4; 
    name="${sample_ID}."fc".L"lane".R2.fastq"; 
    print > name;  
    for (i = 1; i <= 3; i++) {getline; print > name }}' &

    wait 

    ###################################
    # BWA MEM
    for r1 in *.R1.fastq

    do

    r2=\${r1//".R1."/".R2."}

    basename=`echo \$r1 | awk -F"/" '{print \$NF}' | awk -F".R1." '{print \$1}'`
    fc_name=`echo \$basename | awk -F"." '{print \$(NF-1)}'`
    lane=`echo \$basename | awk -F"." '{print \$(NF)}'`
    id="\$fc_name.\${lane}"
    echo \${basename} \${fc_name} \$lane \${id}

    bwa mem -M -t ${task.cpus} -R "@RG\\tID:\${id}\\tPL:Illumina\\tSM:${sample_ID}\\tLB:${sample_ID}" \
	${params.reference_genome} \${r1} \${r2} | \
	picard SortSam -I /dev/stdin -SORT_ORDER queryname \
	-O ${sample_ID}.\${id}_map.bam
    
    done
    
m    ###################################
    # MARK DUPLICATES

    ls -latr ${sample_ID}.*_map.bam | awk '{print "-I "\$9}' > ${sample_ID}_mapped_bams.list.txt

    # sorted automatically by markduplicates spark
    # SB: optical-duplicate-pixel-distance 2500  not sure this number is correct, but google is also uncertain
    # SB: in my tests there is no gain in using more than 8 cores --spark-master local[8] 
    # using 16 as they are allocated anyway 

    gatk MarkDuplicatesSpark --arguments_file ${sample_ID}_mapped_bams.list.txt  ${params.markDuplicatesOptions} \
        -M ${sample_ID}_duplicate_metrics \
	-O ${sample_ID}_dedup_sorted.bam  --spark-master local[${task.cpus}]

    picard SetNmMdAndUqTags  R=${params.reference_genome} I=${sample_ID}_dedup_sorted.bam \
        O=${sample_ID}_merged_dedup_sorted_tagged.bam  
    samtools index ${sample_ID}_merged_dedup_sorted_tagged.bam

    """

}


process alignment_metrics{


    memory { 10.GB * task.attempt }
    cpus 2
    maxRetries 3
    
    input:
    tuple val(sample_ID),path(tagged_bam)
    path(resource_dir)
    
    output:
    publishDir "${params.project_name}/metrics/", mode: 'copy'
    path("${sample_ID}_merged_AlignmentMetrics.txt")
    path("${sample_ID}_validation.txt")

    script:
    """
    picard CollectAlignmentSummaryMetrics R=${params.reference_genome} \
        I=${tagged_bam[0]} \
        O=${sample_ID}_merged_AlignmentMetrics.txt
    picard ValidateSamFile -I ${tagged_bam[0]} -O ${sample_ID}_validation.txt -M SUMMARY
    """
    

}

process base_recalibrator{

    memory { 10.GB * task.attempt }
    cpus 2
    maxRetries 3

    input:
    tuple val(sample_ID), path(tagged_bam_bai)
    path(interval_list)
    path(resource_dir)

    output:
    publishDir "${params.project_name}/bqsr_bam/" , mode: 'copy'
    tuple val(sample_ID),path("${sample_ID}_bqsr.ba*"), emit: bqsr_bam_bai

    script:
    """
    
    bam_file=`ls ${tagged_bam_bai} | grep -v bai`
    gatk --java-options "-Xmx16g -Xms16g" \
        BaseRecalibrator \
        -I \${bam_file} \
        -R ${params.reference_genome} \
	-L ${interval_list} \
	--known-sites ${params.BQSR_known_sites.join(" --known-sites ")} \
        --ip ${params.interval_padding} \
        -O ${sample_ID}_recal_data.table
    
    gatk --java-options "-Xmx16g -Xms16g" \
        ApplyBQSR \
        -R ${params.reference_genome} \
	-I \${bam_file} \
        --bqsr-recal-file ${sample_ID}_recal_data.table  \
	-L ${interval_list} \
	--ip  ${params.interval_padding} -O ${sample_ID}_bqsr.bam
 
     
    """

}


process haplotype_caller_gvcf {

    cpus 6
    memory '20 GB'

    input:
    tuple val(sample_ID),path(bqsr_bam_bai)
    path(interval_list)
    path(resource_dir)
    
    output:
    publishDir "${params.project_name}/vcf/", mode: 'copy'
    path("${sample_ID}.g.vcf.gz*"), emit: gvcf
    //path("${sample_ID}.g.vcf.gz.tbi"), emit: vcf_index
    
    script:
    """
    bam_file=`ls ${bqsr_bam_bai} | grep bam`
    gatk --java-options "-Xmx16g -Xms16g" \
	HaplotypeCaller \
        -R ${params.reference_genome} \
        -I \${bam_file} --native-pair-hmm-threads ${task.cpus} \
	-L ${interval_list}\
	 --ip ${params.interval_padding} \
	-O ${sample_ID}.g.vcf.gz \
        -ERC GVCF
    """

}



process genomicsDBImport {

    cpus 1
    memory '75 GB'
    
    input:
    path(interval_list)
    path(gvcfs)
    path(resource_dir)
    
    publishDir "${params.project_name}/joint_call/", mode: 'copy'

    output:
    path("*scatter.vcf.gz*"), emit: chr_vcf
    
    script:
    """
    # list gvcf files
    ls *.gz | awk -F  ".g.vcf.gz" 'BEGIN{OFS="\\t"} {print \$1,\$0}'  >map_file.txt

    ii=`basename -s ".interval_list" ${interval_list}`
    # loop over chromosomes 
    mkdir tmp_\${ii}

    gatk --java-options "-Xmx50g -Xms30g" GenomicsDBImport \
	--sample-name-map map_file.txt \
	--genomicsdb-workspace-path genomicsDB_\${ii} \
	--batch-size 50 \
	--consolidate true \
	--tmp-dir tmp_\${ii} \
	-L ${interval_list} \
	--reader-threads ${task.cpus}


    gatk --java-options "-Xmx50g -Xms30g" GenotypeGVCFs \
	-R ${params.reference_genome} \
	-V gendb://\${PWD}/genomicsDB_\${ii} \
        -O \${ii}.scatter.vcf.gz \
	--tmp-dir tmp_\${ii}  

    """
}


process gather {

    cpus 2
    memory '200 GB'
    clusterOptions = "--partition=boost --qos=boost"
    
    input:
    path(vcfs)
    //path(vcfs_index)

    publishDir "${params.project_name}/joint_call/", mode: 'copy'
    
    output:
    path("allChr.vcf.gz*"), emit: all_vcf
    path("allChr_sitesonly.vcf.gz*"), emit: all_vcf_sites
    
    script:
    """
    for el in `ls ${vcfs} | grep -v .tbi | sort -V`; do echo -I \${el}; done > lista.txt

    mkdir tmp_gather
    picard GatherVcfs \
	--arguments_file lista.txt  \
	-O allChr.vcf.gz \
	--TMP_DIR tmp_gather

    tabix -f allChr.vcf.gz 

    gatk MakeSitesOnlyVcf \
        -I allChr.vcf.gz \
        -O allChr_sitesonly.vcf.gz \
        --TMP_DIR tmp_gather           

    tabix -f allChr_sitesonly.vcf.gz

    """
}


process variantRecalibratorSNP {

    cpus 2
    memory '80 GB'
    
    input:
    path(vcf)
    path(resource_dir)
    

    publishDir "${params.project_name}/joint_call/", mode:'copy'
    
    output:
    path("allChr_sitesonly_SNP.recal*"), emit: allChr_snp_recal
    path("allChr_sitesonly_SNP.tranches*"), emit: allChr_snp_tranches
    path("allChr_sitesonly_SNP.plots.R*"), emit: allChr_snp_plots
    
    script:
    """

    mkdir snp_tmp
    PATH=/usr/bin/:\${PATH} # otherwise R complains
    vcf_file=`ls ${vcf} | grep -v .tbi`
    gatk --java-options "-Xmx50g -Xms30g" VariantRecalibrator \
	-tranche ${params.recalibrator_tranches.join(" -tranche ")} \
	-R ${params.reference_genome} \
	-V \${vcf_file} -mode SNP \
	-resource:${params.recalibrator_snp_resources.join(" -resource:")} \
	${params.variant_recalibrator_snp_flags} \
	-O allChr_sitesonly_SNP.recal \
	--tmp-dir snp_tmp \
	--tranches-file allChr_sitesonly_SNP.tranches\
        --rscript-file allChr_sitesonly_SNP.plots.R 
    """
}


process variantRecalibratorINDEL {

    cpus 2
    memory '80 GB'
    
    input:
    path(vcf)
    path(resource_dir)
    
    publishDir "${params.project_name}/joint_call/", mode:'copy'
    output:
    path("allChr_sitesonly_INDEL.recal*"), emit: allChr_indel_recal
    path("allChr_sitesonly_INDEL.tranches*"), emit: allChr_indel_tranches
    path("allChr_sitesonly_INDEL.plots.R*"), emit: allChr_indel_plots
    
    script:
    """
    mkdir indel_tmp
    vcf_file=`ls ${vcf} | grep -v .tbi`
    PATH=/usr/bin/:\${PATH} # otherwise R complains
    gatk --java-options "-Xmx50g -Xms30g" VariantRecalibrator \
	-tranche ${params.recalibrator_tranches.join(" -tranche ")} \
	-R ${params.reference_genome} \
	-V \${vcf_file} -mode INDEL\
	-resource:${params.recalibrator_indel_resources.join(" -resource:")} \
	${params.variant_recalibrator_indel_flags} \
	-O allChr_sitesonly_INDEL.recal \
	--tmp-dir indel_tmp \
	--tranches-file allChr_sitesonly_INDEL.tranches\
        --rscript-file allChr_sitesonly_INDEL.plots.R 

    """
}

process applyVqsrINDEL {

    cpus 2
    memory '80 GB'
    
    input:
    path(vcf)
    path(recal)
    path(tranches)
    
    publishDir "${params.project_name}/joint_call/",mode:'copy'
    output:
    path("allChr_INDEL_recalibrated.vcf.gz*"), emit: allChr_indel_vcf

    
    script:
    """
    vcf_file=`ls ${vcf} | grep -v .tbi`
    tranches_file=`ls ${tranches} | grep -v .pdf`
    recal_file=`ls ${recal} | grep -v .idx`

    gatk --java-options "-Xmx80g -Xms50g" \
	ApplyVQSR \
	-V \${vcf_file} \
	--recal-file \${recal_file} \
	--tranches-file \${tranches_file} \
	--truth-sensitivity-filter-level ${params.indel_sensitivity_level} \
	--create-output-variant-index true \
	-mode INDEL \
	-O allChr_INDEL_recalibrated.vcf.gz
    """
}


process applyVqsrSNP {

    cpus 2
    memory '80 GB'
    
    input:
    path(allChr_indel)
    path(recal)
    path(tranches)
    
    publishDir "${params.project_name}/joint_call/", mode:'copy'
    
    output:
    path("allChr_SNP_recalibrated.vcf.gz*"), emit: allChr_snp_vcf
    
    script:
    """
    allChr_indel_file=`ls ${allChr_indel} | grep -v .tbi`
    tranches_file=`ls ${tranches} | grep -v .pdf`
    recal_file=`ls ${recal} | grep -v .idx`

    gatk --java-options "-Xmx80g -Xms50g" \
	ApplyVQSR \
	-V \${allChr_indel_file} \
	--recal-file \${recal_file} \
	--tranches-file \${tranches_file}\
	--truth-sensitivity-filter-level ${params.snp_sensitivity_level} \
	--create-output-variant-index true \
	-mode SNP \
	-O allChr_SNP_recalibrated.vcf.gz
    """
}


process mutect2{

    cpus 4
    memory '20 GB'

    
    input:
    tuple val(sample_name), path(tumor_bam), path(normal_bam)
    path(interval_list)
    path(resource_dir)
    
    // The germline resource is used to get the frequency of a variant allele in the population, thereby providing the prior probability that the sample carries the allele in the germline. This prior is one ingredient in a statistical model for germline variation. When an allele is missing from the germline resource, Mutect2 uses the same model with a very small imputed allele frequency.
    // Mutect2 marks variants that are found in the PON with the “PON” info field, which FilterMutectCalls then uses for filtering. Additionally, Mutect2 considers variants in the PON as inactive by default (this can be changed with the -genotype-pon-sites argument), so most will be silently pre-filtered without ending up in the output. Some PON sites are output because they may appear in an active region containing a non-PON site.
    // f1r2-tar-gz exhibit orientation bias artifacts

    publishDir "${params.project_name}/mutect2/", mode: "copy"

    output:
    tuple val(sample_name),\
    path("${sample_name}.unfiltered.vcf"),\
    path("${sample_name}.unfiltered.vcf.stats"),
    path("${sample_name}.orientation-model.tar.gz")
    
    script:
    """

    tabix -f -p vcf ${params.panel_of_normal_file}
    tumor_bamfile=`ls ${tumor_bam} | grep -v .bai`
    tumor_name=`basename -s ${params.bam_suffix} \${tumor_bamfile}`

    normal_bamfile=`ls ${normal_bam} | grep -v .bai`
    normal_name=`basename -s ${params.bam_suffix} \${normal_bamfile}`
    
    gatk Mutect2 -R ${params.reference_genome} \
        -L ${params.interval_list} --ip ${params.interval_padding}\
        -I \${tumor_bamfile} \
        -tumor \${tumor_name} \
	-I \${normal_bamfile} \
        -normal \${normal_name} \
        -germline-resource ${params.allele_frequency_file} \
        -pon ${params.panel_of_normal_file}  \
        --f1r2-tar-gz ${sample_name}.f1r2.tar.gz \
        -O ${sample_name}.unfiltered.vcf --native-pair-hmm-threads ${task.cpus}


    gatk LearnReadOrientationModel -I ${sample_name}.f1r2.tar.gz -O ${sample_name}.orientation-model.tar.gz 

    """
}


process mutect2_single{

    cpus 4
    memory '20 GB'

    
    input:
    tuple val(sample_ID),path(normal_bam)
    path(interval_list)
    path(resource_dir)
    
    // The germline resource is used to get the frequency of a variant allele in the population, thereby providing the prior probability that the sample carries the allele in the germline. This prior is one ingredient in a statistical model for germline variation. When an allele is missing from the germline resource, Mutect2 uses the same model with a very small imputed allele frequency.
    // Mutect2 marks variants that are found in the PON with the “PON” info field, which FilterMutectCalls then uses for filtering. Additionally, Mutect2 considers variants in the PON as inactive by default (this can be changed with the -genotype-pon-sites argument), so most will be silently pre-filtered without ending up in the output. Some PON sites are output because they may appear in an active region containing a non-PON site.
    // f1r2-tar-gz exhibit orientation bias artifacts

    publishDir "${params.project_name}/mutect2_single/", mode: "copy"

    output:
    path("*.single.vcf.gz"), emit: single_vcf
    path("*.single.vcf.gz.tbi"), emit: single_vcf_tbi
    
    script:
    """
    
    normal_bamfile=`ls ${normal_bam} | grep -v .bai`
    normal_name=`basename -s ${params.bam_suffix} \${normal_bamfile}`
    echo \${normal_name} \${normal_bamfile}

    gatk Mutect2 -R ${params.reference_genome} \
        -L ${params.interval_list} --ip ${params.interval_padding}\
	-I \${normal_bamfile} \
	-tumor \${normal_name} \
        -germline-resource ${params.allele_frequency_file} \
        -max-mnp-distance 0 \
        -O \${normal_name}.single.vcf.gz --native-pair-hmm-threads ${task.cpus}




    """
}


process create_pon{

    cpus 2
    memory '30 GB'

    input:
    path(normal_vcfs)
    path(normal_vcfs_tbi)
    path(interval_list)
    path(resource_dir)
    
    publishDir "${params.project_name}/pon/", mode: "copy"

    output:
    path("*.pon.vcf.gz")
    
    script:
    """

    ii=`basename -s ".interval_list" ${interval_list}`
    # loop over chromosomes 
    mkdir tmp_\${ii}
    
    gatk --java-options "-Xmx20g -Xms10g" GenomicsDBImport -R ${params.reference_genome} \
	-L ${interval_list} --ip ${params.interval_padding} \
	--genomicsdb-workspace-path pon_db_\${ii} --tmp-dir tmp_\${ii} \
	-V ${normal_vcfs.join(" -V ")}

    
    gatk CreateSomaticPanelOfNormals -R ${params.reference_genome} \
	-V gendb://\${PWD}/pon_db_\${ii} -O \${ii}.pon.vcf.gz

    """

}



process learn_orientation{

    input:
    path(f1r2)

    output:
    path("read-orientation-model.tar.gz")
    
    script:
    """
    gatk LearnReadOrientationModel -I ${f1r2} -O read-orientation-model.tar.gz
    """

}

//Run GetPileupSummaries to summarize read support for a set number of known variant sites.

process pileup{

    cpus 2
    memory '20 GB'

    input:
    tuple val(sample_name), path(tumor_bam), path(normal_bam)
    path(resource_dir)

    publishDir "${params.project_name}/mutect2/", mode: "copy"

    output:
    path("${sample_name}.tumor_pileupsummaries.table"), emit: tumor_pileup
    path("${sample_name}.normal_pileupsummaries.table"), emit: normal_pileup

    tuple val(sample_name),\
    path("${sample_name}.segments.table"),\
    path("${sample_name}.calculatecontamination.table"), emit: contamination_table
    
    script:
    """

    tumor_bamfile=`ls ${tumor_bam} | grep -v .bai`
    normal_bamfile=`ls ${normal_bam} | grep -v .bai`

    gatk GetPileupSummaries \
	-I \${tumor_bamfile} \
	-V ${params.pileup_interval} \
	-L ${params.pileup_interval} --ip ${params.interval_padding} \
	-O ${sample_name}.tumor_pileupsummaries.table  

    gatk GetPileupSummaries \
	-I \${normal_bamfile} \
	-V ${params.pileup_interval} \
	-L ${params.pileup_interval} --ip ${params.interval_padding}\
	-O ${sample_name}.normal_pileupsummaries.table  

    gatk CalculateContamination \
        -I  ${sample_name}.tumor_pileupsummaries.table \
        -tumor-segmentation ${sample_name}.segments.table \
	-matched ${sample_name}.normal_pileupsummaries.table \
        -O ${sample_name}.calculatecontamination.table

    """

}



//Finally, pass the learned read orientation model to FilterMutectCallswith the -ob-priors argument:
process filter_mutect{

    publishDir "${params.project_name}/mutect2/", mode: "copy"
    
    input:
    tuple val(sample_name),\
    path(unfiltered),\
    path(stats),\
    path(f1r2),\
    path(segmentation_table),\
    path(contamination_table)

    path(resource_dir)
    
    output:
    path("*.mutect2.filtered.vcf"), emit: filtered_vcf
    //path("${sample_name}.mutect2.unfiltered.vcf"), emit: unfiltered_vcf
    
    script:
    """
    name=`basename -s .unfiltered.vcf ${unfiltered}`
    segmentation_name=`basename -s .segments.table ${segmentation_table}`
    contamination_name=`basename -s .calculatecontamination.table ${contamination_table}`
    echo ${sample_name}
    echo \${name} \${segmentation_name} \${contamination_name}
    gatk FilterMutectCalls -V ${unfiltered} \
	-R ${params.reference_genome} \
        --tumor-segmentation ${segmentation_table} \
        --contamination-table ${contamination_table} \
        --ob-priors ${f1r2} \
	--stats ${stats}  -O ${sample_name}.mutect2.filtered.vcf
    #ouch 
    """
    
}
