#!/usr/bin/env nextflow

/*
* AUTHOR: CERC Genomic Medicine,  Daniel Taliun, PhD <daniel.taliun@mcgill.ca>
* VERSION: 2.0
* YEAR: 2023
*/


process chunk {
	cache "lenient"
	executor "local"
	cpus 1

	input:
	tuple path(sites_vcf), path(sites_vcf_index)

	output:
	path("*.chunk.*.txt")
	path("*.chunks.log")

	publishDir "results/logs/chunk/", pattern: "*.chunks.log", mode: "copy"

	"""
	n_chrom=`bcftools index -s ${sites_vcf} | wc -l`
	if [[ \${n_chrom} -gt 1 ]]; then
		echo "Multiple chromosomes within one reference panel VCF are not allowed." 1>&2
		exit 1
	fi
	chrom=`bcftools index -s ${sites_vcf} | cut -f1`

 	${params.chunk_exec} --input ${sites_vcf} --region \${chrom} --window-size ${params.window_size} --buffer-size ${params.buffer_size} --output \${chrom}.chunks.txt > \${chrom}.chunks.log
	split -l 1 -d --additional-suffix=.txt \${chrom}.chunks.txt \${chrom}.chunk.
	"""	
}


process reference_by_chrom {
	executor "local"
	cpus 1

	input:
	tuple path(vcf), path(vcf_index)

	output:
	tuple stdout, path(vcf), path(vcf_index)

	"""
	n_chrom=`bcftools index -s ${vcf} | wc -l`
	if [[ \${n_chrom} -gt 1 ]]; then
		echo "Multiple chromosomes within one reference panel VCF are not allowed." 1>&2
		exit 1
	fi
	chrom=`bcftools index -s ${vcf} | cut -f1`
	printf "\${chrom}"
	"""
}


process reference_sites_by_chrom {
	executor "local"
	cpus 1

	input:
	tuple path(vcf), path(vcf_index)

	output:
	tuple stdout, path(vcf), path(vcf_index)

	"""
	n_chrom=`bcftools index -s ${vcf} | wc -l`
	if [[ \${n_chrom} -gt 1 ]]; then
		echo "Multiple chromosomes within one reference panel VCF are not allowed." 1>&2
		exit 1
	fi
	chrom=`bcftools index -s ${vcf} | cut -f1`
	printf "\${chrom}"
	"""
}


process HaplotypeCaller {
	//errorStrategy 'retry'
	//maxRetries 3
	cache "lenient"
	cpus 1
	memory "16 GB"
	time "48h"
	//scratch '$SLURM_TMPDIR'
	//stageInMode "copy"

	container "${params.gatkContainer}"
	containerOptions "-B ${params.referenceDir}:/ref"

	input:   
	tuple path(bam), path(bam_index), val(chrom), path(sites_vcf), path(sites_vcf_index)

 	output:
	tuple val(chrom), path("${chrom}.${bam.getSimpleName()}.vcf.gz"), path("${chrom}.${bam.getSimpleName()}.vcf.gz.tbi")
   
	"""
	gatk --java-options -Xmx14G HaplotypeCaller --native-pair-hmm-threads 1 -R /ref/${params.referenceGenome} -L ${sites_vcf} --alleles ${sites_vcf} -I ${bam} -O ${chrom}.${bam.getSimpleName()}.vcf.gz --output-mode EMIT_ALL_ACTIVE_SITES
	"""
}


process join_per_sample_pls {
	cache "lenient"
	cpus 1
	memory "8 GB"
	time "12h"

	input:
	tuple val(chrom), path(vcfs), path(vcfs_indices)

	output:
	tuple val(chrom), path("${chrom}.pls.vcf.gz"), path("${chrom}.pls.vcf.gz.tbi")

	"""
	for f in ${vcfs}; do echo "\${f}"; done | sort  > files.txt
	bcftools merge -i - -m none -l files.txt | bcftools annotate -x INFO | bcftools view -m2 -M2 -Oz -o ${chrom}.pls.vcf.gz
	bcftools index -t ${chrom}.pls.vcf.gz
	"""
}


process study_by_chrom {
	executor "local"
	cpus 1

	input:
	set file(vcf), file(vcf_index) from Channel.fromPath(params.study_vcf).map{ vcf -> [ vcf, vcf + ".tbi" ] }

	output:
	tuple stdout, file(vcf), file(vcf_index) into study

	"""
	tabix -l ${vcf} | tr "\n" ","
	"""
}


process impute_chunks {
	//errorStrategy "retry"
	//maxRetries 3
	cache "lenient"
	cpus 1
	memory "4 GB"
	time "12h"
			
	input:
	tuple val(chrom), path(pls_vcf), path(pls_vcf_index), path(ref_vcf), path(ref_vcf_index), path(chunk)
	
	output:
	tuple val(chrom), path("*.imputed.bcf")
	path("*.imputed.log")

	publishDir "results/logs/impute/", pattern: "*.imputed.log", mode: "copy"
	
	"""
	id=`head -n1 ${chunk} | cut -f1`
	irg=`head -n1 ${chunk} | cut -f3`
	org=`head -n1 ${chunk} | cut -f4`
	chrom_no_prefix=`echo "${chrom}" | sed "s/chr//"`
	
	${params.phase_exec} --input ${pls_vcf} --reference ${ref_vcf} --map ${params.glimpse_maps}/chr\${chrom_no_prefix}.b[0-9]*.gmap.gz --input-region \${irg} --output-region \${org} --output ${chrom}.\${id}.imputed.bcf --log ${chrom}.\${id}.imputed.log
	"""
}


process ligate_chunks {
	//errorStrategy "retry"
	//maxRetries 3
	cache "lenient"
	cpus 1
	memory "4 GB"
	time "12h"
	
	input:
	tuple val(chrom), path(imputed_vcfs)
	
	output:
	tuple path("${chrom}.imputed.bcf"), path("${chrom}.imputed.bcf.csi")
	path("${chrom}.ligate.log")

	publishDir "results", pattern: "${chrom}.imputed.bcf*", mode: "copy"
	publishDir "results/logs/ligate", pattern: "${chrom}.ligate.log", mode: "copy"

	"""
	for f in ${imputed_vcfs}; do bcftools index \${f}; done
	for f in ${imputed_vcfs}; do echo "\${f}"; done | sort -V > files_list.txt
	${params.ligate_exec} --input files_list.txt --output ${chrom}.imputed.bcf --log ${chrom}.ligate.log
	"""
}


workflow {
	chunks = chunk(Channel.fromPath(params.reference_sites_vcfs).map{ vcf -> [vcf, vcf + (vcf.getExtension() == "bcf" ? ".csi" : ".tbi")] })
	
	study_bams = Channel.fromPath(params.study_bams).map { file -> [file, file + (file.getExtension() == "bam" ? ".bai" : ".crai")] }
	reference_vcfs = reference_by_chrom(Channel.fromPath(params.reference_vcfs).map{ file -> [file, file + (file.getExtension() == "bcf" ? ".csi" : ".tbi")] })
	reference_sites_vcfs = reference_sites_by_chrom(Channel.fromPath(params.reference_sites_vcfs).map{ file -> [file, file + (file.getExtension() == "bcf" ? ".csi" : ".tbi")] }) 


	per_sample_pls = HaplotypeCaller(study_bams.combine(reference_sites_vcfs))
	all_pls = join_per_sample_pls(per_sample_pls.groupTuple())

	imputed_chunks = impute_chunks(all_pls.join(reference_vcfs).combine(chunks[0].flatten().map{ file -> [file.getSimpleName(), file] }, by: 0))

	ligate_chunks(imputed_chunks[0].groupTuple())
}
