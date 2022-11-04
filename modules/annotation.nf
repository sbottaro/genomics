
process vep_annotation {

    cpus 4
    memory '20 GB'
    //clusterOptions = "--nodelist=node7"

    input:
    path(vcf)
    path(vep_resource_dir)
    
    output:
    publishDir "${params.project_name}/annotated_vcf/", mode: 'copy'
    path("*_annotated.vcf.gz*"), emit: vcf
    //path("${sample_ID}.g.vcf.gz.tbi"), emit: vcf_index
    
    script:

    """

    for vcf_file in `ls ${vcf} | grep -v ".tbi"`
    do


    if [[ "\${vcf_file}" = *.gz ]]; then
    bb=`basename -s .vcf.gz \${vcf_file}`
    tabix -f \${vcf_file} 
    elif [[ "\${vcf_file}" = *.vcf ]]; then
    bb=`basename -s .vcf.gz \${vcf_file}`                       
    else 
	printf "Not a valid vcf file, exiting." 
    false 
    fi


    /opt/vep/src/ensembl-vep/vep --cache \
	--refseq \
	-i \${vcf_file}\
	-o \${bb}_annotated.vcf.gz \
	--use_transcript_ref \
	--dir_cache /opt/vep/.vep \
	--dir_plugins /opt/vep/.vep/Plugins \
	--species homo_sapiens \
	--use_transcript_ref \
	--assembly GRCh38 \
	--compress_output bgzip \
	--fork ${task.cpus} \
	--offline \
	--fasta /opt/vep/.vep/homo_sapiens_refseq/104_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \
	--buffer_size 5000 --vcf \
	--biotype --hgvs --symbol --canonical --regulatory --gene_phenotype --variant_class \
	--plugin ${params.vep_plugins.join(" --plugin ")}	

    done
    """
}

//	--vcf 
//
//
