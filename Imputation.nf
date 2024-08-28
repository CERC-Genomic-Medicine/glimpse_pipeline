#!/usr/bin/env nextflow

/*
* AUTHOR: CERC Genomic Medicine,  Daniel Taliun, PhD <daniel.taliun@mcgill.ca>
* VERSION: 3.0
* YEAR: 2024
*/


process chunk {
	cache "lenient"
	executor "local"
	cpus 1

	input:
	tuple val(chrom), path(sites_vcf), path(sites_vcf_index)

	output:
	path("*.chunk.*.txt")
	path("*.chunks.log")

	publishDir "results/logs/chunk/", pattern: "*.chunks.log", mode: "copy"

	"""
	# Chromosome X in the reference panel is assumed to be split into chrX_nonpar, chrX_par1 and chrX_par2 regions. Each of these three regions are chunked separately.
	${params.chunk_exec} --input ${sites_vcf} --region ${chrom.split("_")[0]} --window-size ${params.window_size} --buffer-size ${params.buffer_size} --output ${chrom}.chunks.txt > ${chrom}.chunks.log
	split -l 1 -d --additional-suffix=.txt ${chrom}.chunks.txt ${chrom}.chunk.
	"""	
}


process HaplotypeCaller {
	errorStrategy 'retry'
	maxRetries 3
	cache "lenient"
	cpus 1
	memory "4 GB"
	time "8h"
	//scratch '$SLURM_TMPDIR'
	//stageInMode "copy"

	container "${params.gatkContainer}"
	containerOptions "-B ${params.referenceDir}:/ref"

	input:   
	tuple path(bam), path(bam_index), val(chrom), path(sites_vcf), path(sites_vcf_index)

 	output:
	tuple val(chrom), path("${chrom}.${bam.getSimpleName()}.vcf.gz"), path("${chrom}.${bam.getSimpleName()}.vcf.gz.tbi")
   
	"""
	gatk --java-options "-Xmx3G" HaplotypeCaller --native-pair-hmm-threads 1 -R /ref/${params.referenceGenome} -L ${sites_vcf} --alleles ${sites_vcf} -I ${bam} -O ${chrom}.${bam.getSimpleName()}.vcf.gz --output-mode EMIT_ALL_ACTIVE_SITES
	"""
}


process join_per_sample_pls {
	cache "lenient"
	cpus 1
	memory "8 GB"
	time "72h"

	input:
	tuple val(chrom), path(vcfs), path(vcfs_indices)

	output:
	tuple val(chrom), path("${chrom}.pls.vcf.gz"), path("${chrom}.pls.vcf.gz.tbi")

        publishDir "results/joined_per_sample_pls", pattern: "${chrom}.pls.vcf*", mode: "copy" 

	"""
	for f in ${vcfs}; do echo "\${f}"; done | sort  > files.txt
	bcftools merge -i - -m none -l files.txt -Ou | bcftools annotate -x INFO -Ou | bcftools view -m2 -M2 -Oz -o ${chrom}.pls.vcf.gz
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

	# chrX_nonpar -> X, chrX_par1 -> X_par1, chrX_par2 -> X_par2, chr1 -> 1, ..., chr22 -> 22
	map_file_name=`echo "${chrom}" | sed "s/^chr//" | sed "s/_nonpar\$//"`

        # GLIMPSE1 imputation can be negatively affected when using GLs computed for indels. We keep only SNPs.
        bcftools view -v snps ${pls_vcf} \${irg} -Ob -o only_snps_input.bcf
	bcftools index only_snps_input.bcf

	# For GLIMPSE1, to impute INDELs which are present only in reference, we need to use `--impute-reference-only-variants` flag.
	if [ "${chrom}" = "chrX_nonpar" ]; then
		${params.phase_exec} --impute-reference-only-variants --samples-file ${params.study_sample_ploidy} --input only_snps_input.bcf --reference ${ref_vcf} --map ${params.glimpse_maps}/chr\${map_file_name}.b[0-9]*.gmap.gz --input-region \${irg} --output-region \${org} --output ${chrom}.\${id}.imputed.bcf --log ${chrom}.\${id}.imputed.log
	else
		${params.phase_exec} --impute-reference-only-variants --input only_snps_input.bcf --reference ${ref_vcf} --map ${params.glimpse_maps}/chr\${map_file_name}.b[0-9]*.gmap.gz --input-region \${irg} --output-region \${org} --output ${chrom}.\${id}.imputed.bcf --log ${chrom}.\${id}.imputed.log
	fi
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
	// We assume that the chromosome name is encoded in the file prefix e.g. chr1.name.vcf.gz, chr2.name.vcf.gz, ..., chrX_nonpar.name.vcf.gz, chrX_par1.name.vcf.gz, chrX_par2.name.vcf.gz
	reference_vcfs = Channel.fromFilePairs("${params.reference_vcfs}", size: -1) { file -> file.getName().replaceAll("(.csi|.tbi)\$", "") }.map {it -> [ (it[0] =~ /chrX{1}[.]{1}|chrX_nonpar|chrX_par1|chrX_par2|chr[1-9][0-9]?/)[0], it[1][0], it[1][1]] }
	reference_sites_vcfs = Channel.fromFilePairs("${params.reference_sites_vcfs}", size: -1) { file -> file.getName().replaceAll("(.csi|.tbi)\$", "") }.map {it -> [ (it[0] =~ /chrX{1}[.]{1}|chrX_nonpar|chrX_par1|chrX_par2|chr[1-9][0-9]?/)[0], it[1][0], it[1][1]] }

	// Pair *.bam and *.bai (or *.cram and *.crai) files by the filename prefix, and make sure that *.bam (or *.cram) file is first in the emitted pair 
	study_bams = Channel.fromFilePairs("${params.study_bams}", flat: true, size: 2).map(it -> ((it[1].getExtension() == "bam") || (it[1].getExtension() == "cram")) ? [it[1], it[2]] : [it[2], it[1]])

	per_sample_pls = HaplotypeCaller(study_bams.combine(reference_sites_vcfs))
	all_pls = join_per_sample_pls(per_sample_pls.groupTuple())

        chunks = chunk(Channel.fromFilePairs("${params.reference_sites_vcfs}", size: -1) { file -> file.getName().replaceAll("(.csi|.tbi)\$", "") }.map {it -> [ (it[0] =~ /chrX{1}[.]{1}|chrX_nonpar|chrX_par1|chrX_par2|chr[1-9][0-9]?/)[0], it[1][0], it[1][1]] })

	imputed_chunks = impute_chunks(all_pls.join(reference_vcfs).combine(chunks[0].flatten().map { file -> [file.getSimpleName(), file] }, by: 0))
	
	ligate_chunks(imputed_chunks[0].groupTuple())
}
