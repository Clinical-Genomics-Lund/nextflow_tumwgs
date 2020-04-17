#!/usr/bin/env nextflow

// GENERAL PATHS //
OUTDIR = params.outdir+'/'+params.subdir
CRONDIR = params.crondir

// SENTIEON CONFIGS //
K_size      = 100000000
bwa_num_shards = params.bwa_shards
bwa_shards = Channel.from( 0..bwa_num_shards-1 )
genomic_num_shards = params.genomic_shards_num

// FASTA //
genome_file = params.genome_file

// VEP REFERENCES AND ANNOTATION DBS //
CADD = params.CADD
VEP_FASTA = params.VEP_FASTA
MAXENTSCAN = params.MAXENTSCAN
VEP_CACHE = params.VEP_CACHE
GNOMAD = params.GNOMAD
GERP = params.GERP
PHYLOP =  params.PHYLOP
PHASTCONS = params.PHASTCONS

// ANNOTATION DBS GENERAL //
CLINVAR = params.CLINVAR


PON = [F: params.GATK_PON_FEMALE, M: params.GATK_PON_MALE]

group_id = "HEJ"

csv = file(params.csv)
mode = csv.countLines() > 2 ? "paired" : "unpaired"
println(mode)

Channel
    .fromPath(params.csv)
    .splitCsv(header:true)
    .map{ row-> tuple(row.group, row.id, file(row.read1), file(row.read2)) }
    .into { fastq_sharded; fastq; vcf_info }

Channel
    .fromPath(params.csv)
    .splitCsv(header:true)
    .map{ row-> tuple(row.id, row.diagnosis, row.read1, row.read2) }
    .set{ qc_extra }

Channel
    .fromPath(params.csv)
    .splitCsv(header:true)
    .map{ row-> tuple(row.group, row.id, row.sex, row.type) }
    .set { meta_gatkcov }

Channel
    .fromPath(params.csv)
    .splitCsv(header:true)
    .map{ row-> tuple(row.group, row.id, row.type) }
    .into { meta_manta ; meta_concatVCF; meta_vardict; meta_freebayes; meta_dnascope; meta_aggregate}



// Split bed file in to smaller parts to be used for parallel variant calling
Channel
    .fromPath("${params.intersect_bed}")
    .ifEmpty { exit 1, "Regions bed file not found: ${params.intersect_bed}" }
    .splitText( by: 20000, file: 'bedpart.bed' )
    .into { beds_freebayes; beds_vardict }



if(genome_file ){
    bwaId = Channel
            .fromPath("${genome_file}.bwt")
            .ifEmpty { exit 1, "BWA index not found: ${genome_file}.bwt" }
}


Channel
    .fromPath(params.genomic_shards_file)
    .splitCsv(header:false)
    .into { shards1; shards2; shards3; shards4; shards5; }

// A channel to pair neighbouring bams and vcfs. 0 and top value removed later
// Needs to be 0..n+1 where n is number of shards in shards.csv
Channel
    .from( 0..(genomic_num_shards+1) )
    .collate( 3,1, false )
    .into{ shardie1; shardie2 }


// Align fractions of fastq files with BWA
process bwa_align_sharded {
	cpus 50
	memory '64 GB'

	input:
		set val(shard), val(group), val(id), r1, r2 from bwa_shards.combine(fastq_sharded)

	output:
		set val(id), file("${id}_${shard}.bwa.sort.bam"), file("${id}_${shard}.bwa.sort.bam.bai") into bwa_shards_ch

	when:
		params.shardbwa

	"""
	sentieon bwa mem -M \\
		-R '@RG\\tID:${id}\\tSM:${id}\\tPL:illumina' \\
		-K $K_size \\
		-t ${task.cpus} \\
		-p $genome_file '<sentieon fqidx extract -F $shard/$bwa_num_shards -K $K_size $r1 $r2' | sentieon util sort \\
		-r $genome_file \\
		-o ${id}_${shard}.bwa.sort.bam \\
		-t ${task.cpus} --sam2bam -i -
	"""
}

// Merge the fractioned bam files
process bwa_merge_shards {
	cpus 50

	input:
		set val(id), file(shard), file(shard_bai) from bwa_shards_ch.groupTuple()

	output:
		set id, file("${id}_merged.bam"), file("${id}_merged.bam.bai") into merged_bam, qc_merged_bam

	when:
		params.shardbwa
    
	script:
		bams = shard.sort(false) { a, b -> a.getBaseName() <=> b.getBaseName() } .join(' ')

	"""
	sentieon util merge -o ${id}_merged.bam ${bams}
	"""
}

// ALTERNATIVE PATH: Unsharded BWA, utilize local scratch space.
process bwa_align {
	cpus 27
	memory '64 GB'
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'

	input:
		set val(group), val(id), file(r1), file(r2) from fastq

	output:
		set id, file("${id}_merged.bam"), file("${id}_merged.bam.bai") into bam, qc_bam

	when:
		!params.shardbwa

	"""
	sentieon bwa mem \\
		-M \\
		-R '@RG\\tID:${id}\\tSM:${id}\\tPL:illumina' \\
		-t ${task.cpus} \\
		$genome_file $r1 $r2 \\
		| sentieon util sort \\
		-r $genome_file \\
		-o ${id}_merged.bam \\
		-t ${task.cpus} --sam2bam -i -
	"""
}



// Collect information that will be used by to remove duplicate reads.
// The output of this step needs to be uncompressed (Sentieon manual uses .gz)
// or the command will occasionally crash in Sentieon 201808.07 (works in earlier)
process locus_collector {
	cpus 16
	errorStrategy 'retry'
	maxErrors 5

	input:
		set id, file(bam), file(bai), val(shard_name), val(shard) from bam.mix(merged_bam).combine(shards1)

	output:
		set val(id), file("${shard_name}_${id}.score"), file("${shard_name}_${id}.score.idx") into locus_collector_scores
		set val(id), file(bam), file(bai) into merged_bam_id

	"""
	sentieon driver \\
		-t ${task.cpus} \\
		-i $bam $shard \\
		--algo LocusCollector \\
		--fun score_info ${shard_name}_${id}.score
	"""
}



locus_collector_scores
    .groupTuple()
    .join(merged_bam_id)
    .combine(shards2)
    .set{ all_scores }


// Remove duplicate reads
process dedup {
	cpus 16
	cache 'deep'
	errorStrategy 'retry'
	maxErrors 5

	input:
		set val(id), file(score), file(idx), file(bam), file(bai), val(shard_name), val(shard) from all_scores

	output:
		set val(id), file("${shard_name}_${id}.bam"), file("${shard_name}_${id}.bam.bai") into shard_dedup_bam
		set val(group_id), file("${shard_name}_${id}.bam"), file("${shard_name}_${id}.bam.bai") into tndnascope_bams
		set id, file("${shard_name}_${id}_dedup_metrics.txt") into dedup_metrics

	script:
		scores = score.sort(false) { a, b -> a.getBaseName().tokenize("_")[0] as Integer <=> b.getBaseName().tokenize("_")[0] as Integer } .join(' --score_info ')

	"""
	sentieon driver \\
		-t ${task.cpus} \\
		-i $bam $shard \\
		--algo Dedup --score_info $scores \\
		--metrics ${shard_name}_${id}_dedup_metrics.txt \\
		--rmdup ${shard_name}_${id}.bam
	"""

}

shard_dedup_bam
    .groupTuple()
    .into{ all_dedup_bams1; all_dedup_bams2; all_dedup_bams4 }


//merge shards with shard combinations
shards3
    .merge(tuple(shardie1))
    .into{ shard_shard; shard_shard2 }

process dedup_metrics_merge {

	input:
		set id, file(dedup) from dedup_metrics.groupTuple()

	output:
		set id, file("dedup_metrics.txt") into merged_dedup_metrics

	"""
	sentieon driver --passthru --algo Dedup --merge dedup_metrics.txt $dedup
	"""
}


//Collect various QC data: TODO MOVE qc_sentieon to container!
process sentieon_qc {
	cpus 54
	memory '64 GB'
	publishDir "${OUTDIR}/qc", mode: 'copy' , overwrite: 'true'

	input:
		set id, file(bam), file(bai), file(dedup) from qc_bam.mix(qc_merged_bam).join(merged_dedup_metrics)

	output:
		set id, file("${id}.QC") into qc_cdm

	"""
	sentieon driver \\
		-r $genome_file -t ${task.cpus} \\
		-i ${bam} \\
		--algo MeanQualityByCycle mq_metrics.txt \\
		--algo QualDistribution qd_metrics.txt \\
		--algo GCBias --summary gc_summary.txt gc_metrics.txt \\
		--algo AlignmentStat aln_metrics.txt \\
		--algo InsertSizeMetricAlgo is_metrics.txt \\
		--algo WgsMetricsAlgo wgs_metrics.txt
	qc_sentieon.pl $id wgs > ${id}.QC
	"""
}


// Load QC data into CDM (via middleman)
process qc_to_cdm {
	cpus 1
	publishDir "${CRONDIR}/qc", mode: 'copy' , overwrite: 'true'
	
	input:
		set id, file(qc) from qc_cdm
		set id, diagnosis, r1, r2 from qc_extra

	output:
		file("${id}.cdm") into cdm_done

	script:
		parts = r1.split('/')
		idx =  parts.findIndexOf {it ==~ /......_......_...._........../}
		rundir = parts[0..idx].join("/")

	"""
	echo "--run-folder $rundir --sample-id $id --subassay $diagnosis --assay tumwgs --qc ${OUTDIR}/qc/${id}.QC" > ${id}.cdm
	"""
}



process bqsr {
	cpus 16
	errorStrategy 'retry'
	maxErrors 5

	input:
		set val(id), file(bams), file(bai), val(shard_name), val(shard), val(one), val(two), val(three) from all_dedup_bams1.combine(shard_shard)

	output:
		set val(id), file("${shard_name}_${id}.bqsr.table") into bqsr_table

	script:
		combo = [one, two, three]
		combo = (combo - 0) //first dummy value
		combo = (combo - (genomic_num_shards+1)) //last dummy value
		commons = combo.collect{ "${it}_${id}.bam" }   //add .bam to each shardie, remove all other bams
		bam_neigh = commons.join(' -i ')

	"""
	sentieon driver \\
		-t ${task.cpus} \\
		-r $genome_file \\
		-i $bam_neigh $shard \\
		--algo QualCal -k $params.KNOWN1 -k $params.KNOWN2 ${shard_name}_${id}.bqsr.table
	"""
}

// Merge the bqrs shards
process merge_bqsr {
	input:
		set id, file(tables) from bqsr_table.groupTuple()

	output:
		set val(id), file("${id}_merged.bqsr.table") into bqsr_merged

	"""
	sentieon driver \\
		--passthru \\
		--algo QualCal \\
		--merge ${id}_merged.bqsr.table $tables
	"""
}

process merge_dedup_bam {
	cpus 1
	publishDir "${OUTDIR}/bam", mode: 'copy', overwrite: 'true'

	input:
		set val(id), file(bams), file(bais) from all_dedup_bams4

	output:
		set group, id, file("${id}_merged_dedup.bam"), file("${id}_merged_dedup.bam.bai") into cov_bam, freebayes_bam, vardict_bam, manta_bam, dnascope_bam

	script:
		bams_sorted_str = bams.sort(false) { a, b -> a.getBaseName().tokenize("_")[0] as Integer <=> b.getBaseName().tokenize("_")[0] as Integer } .join(' -i ')
		group = "bams"

	"""
	sentieon util merge -i ${bams_sorted_str} -o ${id}_merged_dedup.bam --mergemode 10
	"""
}

bqsr_merged
    .groupTuple()
    .into{ bqsr_merged1; bqsr_merged2;}

all_dedup_bams2
    .join(bqsr_merged1)
    .set{ all_dedup_bams3 }


tndnascope_bams.groupTuple().set { allbams }

all_dedup_bams3
    .combine(shard_shard2).groupTuple(by:5).combine(allbams)
    .set{ tnscope_bam_shards }

// Do somatic SNV calling using TNscope, sharded
process tnscope {
	cpus 16
	errorStrategy 'retry'
	maxErrors 5

	input:
		set id, bams_dummy, bai_dummy, bqsr, val(shard_name), val(shard), val(one), val(two), val(three), val(grid), file(bams), file(bai) from tnscope_bam_shards

	output:
		set grid, file("${shard_name[0]}.vcf"), file("${shard_name[0]}.vcf.idx") into vcf_shard

	script:
		combo = [one[0], two[0], three[0]] // one two three take on values 0 1 2, 1 2 3...30 31 32
		combo = (combo - 0) //first dummy value removed (0)
		combo = (combo - (genomic_num_shards+1)) //last dummy value removed (32)
		commonsT = (combo.collect{ "${it}_${id[0]}.bam" })   //add .bam to each combo to match bam files from input channel
		commonsN = (combo.collect{ "${it}_${id[1]}.bam" })   //add .bam to each combo to match bam files from input channel
		bam_neighT = commonsT.join(' -i ') 
		bam_neighN = commonsN.join(' -i ') 

	"""
	sentieon driver \\
		-t ${task.cpus} \\
		-r $genome_file \\
		-i $bam_neighT -i $bam_neighN $shard \\
		-q ${bqsr[0][0]} -q ${bqsr[1][0]} \\
		--algo TNscope --disable_detector sv --tumor_sample ${id[0]} --normal_sample ${id[1]}  ${shard_name[0]}.vcf
	"""
}


// Merge vcf shards from TNscope
process merge_vcf {
	cpus 16

	input:
		set id, file(vcfs), file(idx) from vcf_shard.groupTuple()
        
	output:
		set group, file("${id}.tnscope.vcf"), file("${id}.tnscope.vcf.idx") into complete_vcf

	script:
		group = "vcfs"
		vcfs_sorted = vcfs.sort(false) { a, b -> a.getBaseName().tokenize("_")[0] as Integer <=> b.getBaseName().tokenize("_")[0] as Integer } .join(' ')

	"""
	sentieon driver \\
		-t ${task.cpus} \\
		--passthru \\
		--algo DNAscope \\
		--merge ${id}.tnscope.vcf $vcfs_sorted
	"""
}



// Do germline SNV calling using DNAscope, sharded
process dnascope {
	cpus 50
	errorStrategy 'retry'
	maxErrors 5

	input:
		set gr, id, file(bam), file(bai) from dnascope_bam.groupTuple()
		set group, smpl_id, type from meta_dnascope.groupTuple()

	output:
		set ID_Tumor, file("${ID_Tumor}_dnascope.vcf.gz") into gvcf_gens

	script:
		Tumor_index = type.findIndexOf{ it == 'tumor' }
		ID_Tumor = smpl_id[Tumor_index]
		tumor_index= id.findIndexOf{it == "$ID_Tumor" }
		bam_tumor = bam[tumor_index]

	"""
	sentieon driver \\
		-t ${task.cpus} \\
		-r $genome_file \\
		-i $bam_tumor \\
		--algo DNAscope --emit_mode GVCF ${ID_Tumor}_dnascope.vcf.gz
	"""
}



// Variant calling with freebayes
process freebayes {
	cpus 1
	errorStrategy 'retry'
	maxErrors 5

	input:
		set val(gr), id, file(bam), file(bai) from freebayes_bam.groupTuple()
		set val(group), smpl_id , val(type) from meta_freebayes.groupTuple()
		each file(bed) from beds_freebayes

	output:
		set val("freebayes"), group , file("freebayes_${bed}.vcf") into vcfparts_freebayes
		
	script:
		if( mode == "paired" ) {
			Tumor_index = type.findIndexOf{ it == 'tumor' }
			ID_Tumor = smpl_id[Tumor_index]
			tumor_index= id.findIndexOf{it == "$ID_Tumor" }
			bam_tumor = bam[tumor_index]

			Normal_index = type.findIndexOf{ it == 'normal' }
			ID_normal = smpl_id[Normal_index]
			normal_index = id.findIndexOf{it == "$ID_normal" }
			bam_normal = bam[normal_index]

			"""
			freebayes -f $genome_file -t $bed --pooled-continuous --pooled-discrete --min-repeat-entropy 1 -F 0.03 $bam_tumor $bam_normal  > freebayes_${bed}.vcf.raw
			#vcffilter -F LowCov -f "DP > 30" -f "QA > 150" freebayes_${bed}.vcf.raw | vcffilter -F LowFrq -o -f "AB > 0.05" -f "AB = 0" | vcfglxgt > freebayes_${bed}.filt1.vcf
			filter_freebayes_somatic_wgs.pl freebayes_${bed}.vcf.raw $ID_Tumor $ID_normal | grep -v 'FAIL_' > freebayes_${bed}.vcf
			"""
		}
		else if( mode == "unpaired" ) {
			"""
			freebayes -f $genome_file -t $bed --pooled-continuous --pooled-discrete --min-repeat-entropy 1 -F 0.03 $bam > freebayes_${bed}.vcf
			"""
		}
}



process vardict {
	cpus 4
	errorStrategy 'retry'
	maxErrors 5

	input:
		set gr, id, file(bam), file(bai) from vardict_bam.groupTuple()
		set val(group), smpl_id , val(type) from meta_vardict.groupTuple()
		each file(bed) from beds_vardict
		//.splitText( by: 150, file: 'minibedpart.bed' )

	output:
		set val("vardict"), group , file("vardict_${bed}.vcf") into vcfparts_vardict
		
	script:
		if( mode == "paired" ) {

			Tumor_index = type.findIndexOf{ it == 'tumor' }
			ID_Tumor = smpl_id[Tumor_index]
			tumor_index= id.findIndexOf{it == "$ID_Tumor" }
			bam_tumor = bam[tumor_index]

			Normal_index = type.findIndexOf{ it == 'normal' }
			ID_normal = smpl_id[Normal_index] 
			normal_index = id.findIndexOf{it == "$ID_normal" }
			bam_normal = bam[normal_index]

			"""
			export JAVA_HOME=/opt/conda/envs/CMD-TUMWGS
			vardict-java -U -th 4 -G $genome_file -f 0.03 -N ${ID_Tumor} -b "${bam_tumor}|${bam_normal}" -c 1 -S 2 -E 3 -g 4 ${bed} | testsomatic.R | var2vcf_paired.pl -N "${ID_Tumor}|${ID_normal}" -f 0.03 > vardict_${bed}.raw.vcf
			filter_vardict_somatic_wgs.pl vardict_${bed}.raw.vcf $ID_Tumor $ID_normal | grep -v 'FAIL_' > vardict_${bed}.vcf
			"""
		}
		else if( mode == "unpaired" ) {
			"""
			export JAVA_HOME=/opt/conda/envs/CMD-TUMWGS
			vardict-java -U -G $genome_file -f 0.03 -N ${id} -b $bam -c 1 -S 2 -E 3 -g 4 $bed | teststrandbias.R | var2vcf_valid.pl -N  ${id} -E -f 0.03 > vardict_${bed}.vcf
			"""
		}
}



// Prepare vcf parts for concatenation
vcfparts_freebayes = vcfparts_freebayes.groupTuple(by:[0,1])
vcfparts_vardict = vcfparts_vardict.groupTuple(by:[0,1])
vcfs_to_concat = vcfparts_freebayes.mix(vcfparts_vardict)




process concatenate_vcfs {
	publishDir "${OUTDIR}/vcf", mode: 'copy', overwrite: true

	input:
		set vc, gr, file(vcfs) from vcfs_to_concat
		

	output:
		set gr, vc, file("${gr}_${vc}.vcf.gz") into concatenated_vcfs

	"""
	vcf-concat $vcfs | vcf-sort -c | gzip -c > ${vc}.concat.vcf.gz
	vt decompose ${vc}.concat.vcf.gz -o ${vc}.decomposed.vcf.gz
	vt index ${vc}.decomposed.vcf.gz
	vt sort -m chrom  ${vc}.decomposed.vcf.gz -o  ${vc}.decomposed.sorted.vcf.gz
	vt normalize ${vc}.decomposed.sorted.vcf.gz -r $genome_file | vt uniq - -o ${gr}_${vc}.vcf.gz
	"""
}


process aggregate_vcfs {
	cpus 1
	publishDir "${OUTDIR}/vcf", mode: 'copy', overwrite: true
	time '40m'

	input:
		set group, vc, file(vcfs) from concatenated_vcfs.groupTuple()
		set g, id, type from meta_aggregate.groupTuple()

	output:
		set group, val("${id[tumor_idx]}"), file("${group}.agg.vcf") into vcf_pon

	script:
		sample_order = id[0]
		if( mode == "paired" ) {
			tumor_idx = type.findIndexOf{ it == 'tumor' || it == 'T' }
			normal_idx = type.findIndexOf{ it == 'normal' || it == 'N' }
			sample_order = id[tumor_idx]+","+id[normal_idx]
		}

	"""
	aggregate_vcf.pl --vcf ${vcfs.sort(false) { a, b -> a.getBaseName() <=> b.getBaseName() }.join(",")} --sample-order ${sample_order} |vcf-sort -c > ${group}.agg.unsorted.vcf
	vcf-sort -c ${group}.agg.unsorted.vcf > ${group}.agg.vcf
	"""
}

process pon_filter {
	publishDir "${OUTDIR}/vcf", mode: 'copy', overwrite: true
	cpus 1
	time '1h'

	input:
		set group, tumor_id, file(vcf) from vcf_pon

	output:
		set group, file("${group}.agg.pon.vcf") into vcf_vep

	script:
		def pons = []
		if( params.PON_freebayes ) { pons.push("freebayes="+params.PON_freebayes) }
		if( params.PON_vardict )   { pons.push("vardict="+params.PON_vardict) }
		def pons_str = pons.join(",")

	"""
	filter_with_pon.pl --vcf $vcf --pons $pons_str --tumor-id $tumor_id > ${group}.agg.pon.vcf
	"""
}


complete_vcf
    .groupTuple()
    .set{ gvcfs }

process gvcf_combine {
	cpus 16

	input:
		set id, file(vcf), file(idx) from gvcfs
		set val(group), val(id), r1, r2 from vcf_info

	output:
		set group, file("${group}.combined.vcf"), file("${group}.combined.vcf.idx") into combined_vcf

	script:
		// Om fler än en vcf, GVCF combine annars döp om och skickade vidare
		if (mode == "family" ) {
			ggvcfs = vcf.join(' -v ')

			"""
			sentieon driver \\
				-t ${task.cpus} \\
				-r $genome_file \\
				--algo GVCFtyper \\
				-v $ggvcfs ${group}.combined.vcf
			"""
		}
		// annars ensam vcf, skicka vidare
		else {
			ggvcf = vcf.join('')
			gidx = idx.join('')

			"""
			mv ${ggvcf} ${group}.combined.vcf
			mv ${gidx} ${group}.combined.vcf.idx
			"""
		}
}

// Splitting & normalizing variants:
process split_normalize {
	cpus 1
	publishDir "${OUTDIR}/vcf", mode: 'copy', overwrite: 'true'

	input:
		set group, file(vcf), file(idx) from combined_vcf

	output:
		set group, file("${group}.norm.uniq.DPAF.vcf") into split_norm, vcf_gnomad

	"""
	vcfbreakmulti ${vcf} > ${group}.multibreak.vcf
	bcftools norm -m-both -c w -O v -f $genome_file -o ${group}.norm.vcf ${group}.multibreak.vcf
	vcfstreamsort ${group}.norm.vcf | vcfuniq > ${group}.norm.uniq.vcf
	wgs_DPAF_filter.pl ${group}.norm.uniq.vcf > ${group}.norm.uniq.DPAF.vcf
	"""

}

// Intersect VCF, exome/clinvar introns
process intersect {

	input:
		set group, file(vcf) from split_norm

	output:
		set group, file("${group}.intersected.vcf") into split_vep, split_cadd, vcf_loqus

	"""
	bedtools intersect -a $vcf -b $params.intersect_bed -u -header > ${group}.intersected.vcf
	"""

}


process annotate_vep {
	container = '/fs1/resources/containers/ensembl-vep_latest.sif'
	publishDir "${OUTDIR}/vcf", mode: 'copy', overwrite: true
	cpus 20
	time '1h'

	input:
		set group, file(vcf) from vcf_vep

	output:
		set group, file("${group}.agg.pon.vep.vcf") into vcf_panel

	"""
	vep -i ${vcf} -o ${group}.agg.pon.vep.vcf \\
		--offline --merged --everything --vcf --no_stats \\
		--fork ${task.cpus} \\
		--force_overwrite \\
		--plugin CADD $params.CADD --plugin LoFtool \\
		--fasta $params.VEP_FASTA \\
		--dir_cache $params.VEP_CACHE --dir_plugins $params.VEP_CACHE/Plugins \\
		--distance 200 \\
		-cache -custom $params.GNOMAD \\
	"""
}

process filter_with_panel_snv {
	publishDir "${OUTDIR}/vcf", mode: 'copy', overwrite: true
	cpus 1
	time '1h'

	input:
		set group, file(vcf) from vcf_panel

	output:
		set group, file("${group}.agg.pon.vep.panel.vcf")

	"""
	filter_with_panel_snv.pl $vcf $params.PANEL_SNV > ${group}.agg.pon.vep.panel.vcf
	"""
}



// Create coverage profile using GATK
process gatkcov {
	publishDir "${OUTDIR}/cov", mode: 'copy' , overwrite: 'true'    
    
	cpus 2
	memory '64 GB'

	input:
		set id, group, file(bam), file(bai), gr, sex, type from cov_bam.join(meta_gatkcov, by:1).groupTuple(by:1)

	output:
		set val("${id[tumor_idx]}"), file("${id[tumor_idx]}.standardizedCR.tsv"), file("${id[tumor_idx]}.denoisedCR.tsv") into cov_gens
		file("${id[tumor_idx]}.modeled.png")

	script:
		tumor_idx = type.findIndexOf{ it == 'tumor' }
		normal_idx = type.findIndexOf{ it == 'normal' }


	"""
	source activate gatk4-env

	gatk CollectReadCounts \\
		-I ${bam[tumor_idx]} -L $params.COV_INTERVAL_LIST \\
		--interval-merging-rule OVERLAPPING_ONLY -O ${bam[tumor_idx]}.hdf5

	gatk --java-options "-Xmx50g" DenoiseReadCounts \\
		-I ${bam[tumor_idx]}.hdf5 --count-panel-of-normals ${PON[sex[tumor_idx]]} \\
		--standardized-copy-ratios ${id[tumor_idx]}.standardizedCR.tsv \\
		--denoised-copy-ratios ${id[tumor_idx]}.denoisedCR.tsv

	gatk PlotDenoisedCopyRatios \\
		--standardized-copy-ratios ${id[tumor_idx]}.standardizedCR.tsv \\
		--denoised-copy-ratios ${id[tumor_idx]}.denoisedCR.tsv \\
		--sequence-dictionary $params.GENOMEDICT \\
		--minimum-contig-length 46709983 --output . --output-prefix ${id[tumor_idx]}

	gatk --java-options "-Xmx50g" CollectAllelicCounts \\
		-L $params.GATK_GNOMAD \\
		-I ${bam[tumor_idx]} \\
		-R $genome_file \\
		-O ${id[tumor_idx]}.allelicCounts.tsv

	gatk --java-options "-Xmx50g" CollectAllelicCounts \\
		-L $params.GATK_GNOMAD \\
		-I ${bam[normal_idx]} \\
		-R $genome_file \\
		-O ${id[normal_idx]}.allelicCounts.tsv

	gatk --java-options "-Xmx40g" ModelSegments \\
		--denoised-copy-ratios ${id[tumor_idx]}.denoisedCR.tsv \\
		--allelic-counts ${id[tumor_idx]}.allelicCounts.tsv \\
		--normal-allelic-counts ${id[normal_idx]}.allelicCounts.tsv \\
		--minimum-total-allele-count-normal 20 \\
		--output . \\
		--output-prefix ${id[tumor_idx]}

	gatk CallCopyRatioSegments \\
		--input ${id[tumor_idx]}.cr.seg \\
		--output ${id[tumor_idx]}.called.seg

	gatk PlotModeledSegments \\
		--denoised-copy-ratios ${id[tumor_idx]}.denoisedCR.tsv \\
		--allelic-counts ${id[tumor_idx]}.hets.tsv \\
		--segments ${id[tumor_idx]}.modelFinal.seg \\
		--sequence-dictionary $params.GENOMEDICT \\
		--minimum-contig-length 46709983 \\
		--output . \\
		--output-prefix ${id[tumor_idx]}
	"""
}

process generate_gens_data {
	publishDir "${OUTDIR}/plot_data", mode: 'copy' , overwrite: 'true'
	tag "$group"
	cpus 1

	input:
		set id, file(gvcf), file(cov_stand), file(cov_denoise) from gvcf_gens.join(cov_gens)

	output:
		set file("${id}.cov.bed.gz"), file("${id}.baf.bed.gz"), file("${id}.cov.bed.gz.tbi"), file("${id}.baf.bed.gz.tbi")

	"""
	generate_gens_data.pl $cov_stand $gvcf $id $params.GENS_GNOMAD
	"""
}


//Somatic Variant Calling - Manat 
process manta{
	publishDir "$OUTDIR/manta" , mode:'copy'
	cpus 20
	memory '64 GB'

	input: 
		set val(gr), id, file(bam), file(bai) from manta_bam.groupTuple()
		set val(group), smpl_id , val(type) from meta_manta.groupTuple()

	output:
		set val(group), file("${group}_manta.vcf") into manta_vcf

	script:
	
	if( mode == "paired" ) {

		Tumor_index = type.findIndexOf{ it == 'tumor' }
		ID_Tumor = smpl_id[Tumor_index]
		tumor_index= id.findIndexOf{it == "$ID_Tumor" }
		bam_tumor = bam[tumor_index]

		Normal_index = type.findIndexOf{ it == 'normal' }
		ID_normal = smpl_id[Normal_index] 
		normal_index = id.findIndexOf{it == "$ID_normal" }
		bam_normal = bam[normal_index]

		"""
		configManta.py \\
		--tumorBam ${bam_tumor} \\
		--normalBam ${bam_normal} \\
		--reference ${genome_file} \\
		--runDir .

		python runWorkflow.py -m local -j ${task.cpus}
		mv ./results/variants/somaticSV.vcf.gz ${group}_manta.vcf.gz
		gunzip ${group}_manta.vcf.gz
		
		"""
		}
	else if( mode == "unpaired" ) {
		"""
		configManta.py \\
			--tumorBam ${bam} \\
			--reference ${genome_file} \\
			--generateEvidenceBam \\
			--region \\
			--runDir .
		python runWorkflow.py -m local -j ${task.cpus}
		
		mv ./results/variants/tumorSV.vcf.gz ${group}_manta.vcf.gz
		gunzip ${group}_manta.vcf.gz 
		"""
		}	
	}


process annotate_manta {
	publishDir "$OUTDIR/manta" , mode:'copy'
	cpus 2
	memory '8 GB'

	input:
		set group, file(vcf) from manta_vcf

	output:
		set group, file("${group}.manta.snpeff.vcf")

	"""
	snpEff -Xmx4g -configOption data.dir=$params.SNPEFF_DIR GRCh37.75 \\
		$vcf > ${group}.manta.snpeff.vcf
	"""
}
