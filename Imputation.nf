

process chunk {
	executor "local"
	cpus 1

	input:
	set file(sites_vcf), file(sites_vcf_index) from Channel.fromPath(params.reference_sites_vcfs).map{ vcf -> [ vcf, vcf + ".csi" ] }

	output:
	tuple stdout, file("*.chunk.*.txt") into chunks
	file "*.chunks.log"

	publishDir "results/logs", pattern: "*.chunks.log", mode: "move"

	"""
	n_chrom=`bcftools index -s ${sites_vcf} | wc -l`
	if [[ \${n_chrom} -gt 1 ]]; then
		echo "Multiple chromosomes within one reference panel VCF are not allowed." 1>&2
		exit 1
	fi
	chrom=`bcftools index -s ${sites_vcf} | cut -f1`

 	${params.chunk_exec} --input ${sites_vcf} --region \${chrom} --window-size ${params.window_size} --buffer-size ${params.buffer_size} --output \${chrom}.chunks.txt > \${chrom}.chunks.log
	split -l 1 -d --additional-suffix=.txt \${chrom}.chunks.txt \${chrom}.chunk.
	printf "\${chrom}"
	"""	
}


process reference_by_chrom {
	executor "local"
	cpus 1

	input:
	set file(vcf), file(vcf_index) from Channel.fromPath(params.reference_vcfs).map{ vcf -> [ vcf, vcf + ".csi" ] }

	output:
	tuple stdout, file(vcf), file(vcf_index) into references

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


chunks = chunks.transpose().combine(references, by: 0)
study = study.flatMap { chroms, vcf, vcf_index -> chroms.split(',').collect { [it, vcf, vcf_index] }} 


process impute_chunks {
	errorStrategy "retry"
	maxRetries 3
	cpus 1
			
	input:
	set val(chrom), file(study_vcf), file(study_vcf_index), file(chunk), file(ref_vcf), file(ref_vcf_index) from study.combine(chunks, by: 0)

	output:
	tuple val(chrom), val("${study_vcf.getBaseName()}"), file("*.imputed.bcf") into imputed
	file "*.imputed.log"

	publishDir "results/logs", pattern: "*.imputed.log", mode: "copy"
	
	"""
	id=`head -n1 ${chunk} | cut -f1`
	irg=`head -n1 ${chunk} | cut -f3`
	org=`head -n1 ${chunk} | cut -f4`
	chrom_no_prefix=`echo "${chrom}" | sed "s/chr//"`
	${params.phase_exec} --input ${study_vcf} --reference ${ref_vcf} --map ${params.glimpse_maps}chr\${chrom_no_prefix}.b37.gmap.gz --input-region \${irg} --output-region \${org} --output ${chrom}.\${id}.${study_vcf.getBaseName()}.imputed.bcf --log ${chrom}.\${id}.${study_vcf.getBaseName()}.imputed.log
	"""
}


process ligate_chunks {
	errorStrategy "retry"
	maxRetries 3
	cpus 1
	
	input:
	set val(chrom), val(base_name), file(imputed_vcf) from imputed.groupTuple(by: [0, 1])

	output:
	tuple val(base_name), file("${chrom}.${base_name}.imputed.bcf"), file("${chrom}.${base_name}.imputed.bcf.csi") into ligated_vcfs
	file "${chrom}.${base_name}.ligate.log"

	publishDir "results", pattern: "${chrom}.${base_name}.imputed.bcf*", enabled: !params.concatenate, mode: "copy"
	publishDir "results/logs", pattern: "*.ligate.log", mode: "copy"

	"""
	for f in ${imputed_vcf}; do bcftools index \${f}; done
	for f in ${imputed_vcf}; do echo "\${f}"; done | sort -V > files_list.txt
	${params.ligate_exec} --input files_list.txt --output ${chrom}.${base_name}.imputed.bcf --log ${chrom}.${base_name}.ligate.log
	bcftools index ${chrom}.${base_name}.imputed.bcf
	"""
}


process concat_chroms {
	errorStrategy "retry"
	maxRetries 3
	cpus 1

	when:
	params.concatenate

	input:
	set val(base_name), file(imputed_vcf), file(imputed_vcf_index) from ligated_vcfs.groupTuple()

	output:
	set file("${base_name}.imputed.bcf"), file("${base_name}.imputed.bcf.csi")

	publishDir "results", mode: "copy"

	"""
	for f in ${imputed_vcf}; do echo \${f}; done | sort -V > files.txt
	bcftools concat -f files.txt -Ob -o ${base_name}.imputed.bcf
	bcftools index ${base_name}.imputed.bcf
	"""
}


