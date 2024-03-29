#!/usr/bin/env nextflow
nextflow.enable.dsl=2

VERSION="0.0"

log.info "===================================================================="
log.info "FASTQC fasta files in folder                                        "
log.info "===================================================================="

params.help = ""
if (params.help) {
  log.info " "
  log.info "USAGE: "
  log.info " "
  log.info "nextflow -c somatic.config run somatic.nf --fasta list.txt"
  log.info " "
    log.info "   list.txt MUST be in the following format "
    log.info "   sampleId    tumor    normal "
    log.info "   sample_01   FULLPATH/tumor_01.ba{m,i}    FULLPATH/normal_01.ba{m,i}"
    log.info "   sample_02   FULLPATH/tumor_02.ba{m,i}    FULLPATH/normal_03.ba{m,i}"
    log.info "   sample_03   FULLPATH/tumor_03.ba{m,i}    FULLPATH/normal_03.ba{m,i}"
    log.info "   ...."
    log.info " "
    log.info " 1. the header MUST be present as described above"
    log.info " 2. all fields MUST be tab separated with no spaces"
    log.info " "
  log.info "===================================================================="
  exit 1
}


//includeConfig 'step_A.config'
include {fastqc; multiqc } from './modules/qc.nf'
include {trim_illumina } from './modules/trim.nf'
include {split_align_mark; alignment_metrics; base_recalibrator } from './modules/gatk.nf'
include {vep_annotation as vep_annotation_mutect } from './modules/annotation.nf'
include {vep_annotation as vep_annotation_strelka } from './modules/annotation.nf'
include {mutect2; learn_orientation; pileup;  filter_mutect} from './modules/gatk.nf'
include {strelka } from './modules/strelka.nf'



workflow qc_pretrimming {
    take:
    fasta_channel
    
    main:
    // fastqc_pre trimming
    fastqc(fasta_channel,"pre")
    //multiqc
    multiqc(fastqc.out.fastqc.collect(),"pre")
    
}

workflow trimming {
    
    take:
    fasta_channel
    
    main:
    // trim reads
    trim_illumina(fasta_channel,params.trimmomatic_file)
    
    // fastqc_post trimming
    fastqc(trim_illumina.out.trimmed_reads,"post")

    //multiqc post trimming 
    multiqc(fastqc.out.fastqc.collect(),"post")


    emit:
    trim_illumina.out.trimmed_reads
}



workflow align_and_calibrate{

    take:
    fasta_channel
    
    main:
    // split into readgroups, align using BWA MEM, mark duplicates, and set tags
    split_align_mark(fasta_channel,\
		     params.resource_dir)

    //split_align_mark.out.bam_bai.collect().view()
    // calculate alignment metrics
    alignment_metrics(split_align_mark.out.bam_bai,params.resource_dir)

    // base recalibrator
    base_recalibrator(split_align_mark.out.bam_bai,\
		      params.interval_list,\
		      params.resource_dir)

    emit:
    base_recalibrator.out.bqsr_bam_bai
}

workflow mutect_pipeline {

    take:
    bam_pair

    main:    
    mutect2(bam_pair,params.interval_list,params.resource_dir)

    // mutect2.stats is needed only because I do not know how to
    // handle channels properly.
    pileup(bam_pair,params.resource_dir)
    //calculate_contamination(pileup.out.tumor_pileup, pileup.out.normal_pileup)

    ch2 = mutect2.out.join(pileup.out.contamination_table,by:0)
    ch2.view()
    filter_mutect(ch2,params.resource_dir)

    emit:
    filter_mutect.out.filtered_vcf
    
}


workflow strelka_pipeline {

    take:
    bam_pair

    main:    
    strelka(bam_pair,params.interval_list,params.resource_dir)

    emit:
    strelka.out.vcfs
    
}


workflow{

    
    bam_files = Channel.fromPath( params.bam )\
	    .splitCsv( header: true, sep: '\t' )\
	    .map { row -> tuple( row.sampleId, file(row.tumor), file(row.normal) ) }
    bam_files.view()
    
    // read file with fasta names
    //bam_files =  Channel.fromList(file(params.bam).readLines())
    //bam_files.view()
    //strelka_pipeline(bam_files)
    mutect_pipeline(bam_files)
    
    //vep_annotation_strelka(strelka_pipeline.out,params.vep_resource_dir)
    //vep_annotation_mutect(mutect_pipeline.out,params.vep_resource_dir)
    
}


