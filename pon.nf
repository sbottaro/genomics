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
include {mutect2_single; create_pon} from './modules/gatk.nf'
include {strelka } from './modules/strelka.nf'



workflow mutect_pon {

    take:
    bam

    main:    
    mutect2_single(bam,params.interval_list,params.resource_dir)

    chr_channel = Channel.fromPath(params.pon_lists )
    create_pon(mutect2_single.out.single_vcf.collect(),\
	       mutect2_single.out.single_vcf_tbi.collect(),\
	       chr_channel,params.resource_dir)
    

    emit:
    create_pon.out
    
}



workflow{

    
    //bam_files = Channel.fromPath( params.bam )\
//	    .splitCsv( header: true, sep: '\t' )\
//	    .map { row -> tuple( row.sampleId, file(row.normal) ) }
    //bam_files.view()
    
    // read file with fasta names
    bam_files =   Channel.fromFilePairs(file(params.bam).readLines())
    bam_files.view()
    //strelka_pipeline(bam_files)
    mutect_pon(bam_files)
    
    //vep_annotation_strelka(strelka_pipeline.out,params.vep_resource_dir)
    //vep_annotation_mutect(mutect_pipeline.out,params.vep_resource_dir)
    
}


