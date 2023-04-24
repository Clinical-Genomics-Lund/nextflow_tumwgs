#!/usr/bin/env nextflow

/*
* 'TumWgs' - A tumor whole genome sequencing pipeline that takes use senteion genomic distributed workflow  mode.   
*/ 

nextflow.enable.dsl = 2

/*
* General parameters 
*/



log.info """\
======================================================================
TUMWGS -NF v2
======================================================================
outdir                  :       $params.outdir
subdir                  :       $params.subdir
crondir                 :       $params.crondir
genome                  :       $params.genome_file
csv                     :       $params.csv
panel                   :       $params.panel
bwa_num_shards          :       $params.bwa_shards
genome_num              :       $params.genomic_shards_num
PON_freebayes           :       $params.PON_freebayes
PON_vardict             :       $params.PON_vardict
Intersect_bed           :       $params.intersect_bed
CADD                    :       $params.CADD 
VEP_FASTA               :       $params.VEP_FASTA
VEP_CACHE               :       $params.VEP_CACHE
GNOMAD                  :       $params.GNOMAD
SNV_HARD_FILTER         :       $params.SNV_HARD_FILTER
=====================================================================
"""

/*
* Import modules
*/

include {   COPY_FASTQ;
            BWA_ALIGN_SHARDED;
            DELETE_FASTQ;
            BWA_MERGE_SHARDS;
            BAM_CRAM_ALL;
            LOCUS_COLLECTOR;
            DEDUP;
            DEDUP_METRICS_MERGE;
            SENTIEON_QC;
            BQSR;
            MERGE_BQSR;
            MERGE_DEDUP_CRAM;
            QC_TO_CDM;
            TNSCOPE;
            MERGE_VCF;
            DNASCOPE_TUM;
            DNASCOPE_NOR;
            FREEBAYES;
            VARDICT;
            CONCATENATE_VCFS;
            AGGREGATE_VCFS;
            PON_FILTER;
            GVCF_COMBINE;
            SPLIT_NORMALIZE;
            INTERSECT;
            ANNOTATE_VEP;
            FILTER_WITH_PANEL_SNV;
            MANTA;
            ANNOTATE_MANTA;
            FILTER_WITH_PANEL_FUSIONS;
            GATKCOV_BAF;
            GATKCOV_COUNT_TUM;
            GATKCOV_CALL_TUM;
            GATKCOV_COUNT_NOR;
            GATKCOV_CALL_NOR;
            CNVS_ANNOTATE;
            FILTER_WITH_PANEL_CNVS;
            GENERATE_GENS_DATA;
            GENERATE_GENS_DATA_NOR;
            COYOTE
        } from './modules/sentieon/main.nf' params(params) 
    

/*  Subworkflows for the pipelines could be made modular  */

workflow sentieon_workflow {            
    take:
        K_size
        genomeShards
        bwa_num_shards
        bwa_shards
        fastq_sharded

    main:
     
    if (params.copy) {
         COPY_FASTQ (    fastq_sharded   )
         BWA_ALIGN_SHARDED ( K_size, 
                            bwa_num_shards,
                            bwa_shards.combine(COPY_FASTQ.out)   )
         BWA_MERGE_SHARDS (  BWA_ALIGN_SHARDED.out.groupTuple()  )
         DELETE_FASTQ ( BWA_MERGE_SHARDS.out, COPY_FASTQ.out )
        } else {
         BWA_ALIGN_SHARDED ( K_size, 
                            bwa_num_shards,
                            bwa_shards.combine(fastq_sharded)   )
         BWA_MERGE_SHARDS (  BWA_ALIGN_SHARDED.out.groupTuple()  )
        }
        BAM_CRAM_ALL (  params.genome_file,
                        BWA_MERGE_SHARDS.out.groupTuple()    )        
        
        LOCUS_COLLECTOR (   params.genome_file,
                            BAM_CRAM_ALL.out.combine(shards)    )

        dedupInput  =   LOCUS_COLLECTOR.out.groupTuple().join(BAM_CRAM_ALL.out).combine(shards)
        DEDUP ( params.genome_file, dedupInput )
        
        dedupMetricesInput  =   DEDUP.out[1].groupTuple()
        DEDUP_METRICS_MERGE ( dedupMetricesInput  )

        SENTIEON_QC (   params.genome_file, 
                        BAM_CRAM_ALL.out.join(DEDUP_METRICS_MERGE.out))

        bqsrInput = DEDUP.out[0].groupTuple().combine(shards)
        BQSR (  genomeShards,
                params.genome_file,
                params.KNOWN,
                bqsrInput   )

        MERGE_BQSR (    BQSR.out.groupTuple()   )
        MERGE_DEDUP_CRAM (  params.genome_file,
                            DEDUP.out[0].groupTuple()   )

    emit:
        bam = BWA_MERGE_SHARDS.out
        cram = MERGE_DEDUP_CRAM.out
        dedup = DEDUP.out[0]
        bqsr = MERGE_BQSR.out
        sentieonqc = SENTIEON_QC.out

}

workflow cdm_qc_workflow {
    take:
        cron
        sentieonqc
    main:
        QC_TO_CDM ( sentieonqc.join(cron)  )
}

workflow tnscope_workflow {
    take:
        genomeShards
        metaid
        bqsr
        dedup
        shards
    main:
        bqsr_1      =   bqsr.groupTuple().first()
        bqsr_2      =   bqsr.groupTuple().last()
        sample1     =   dedup.groupTuple().first()
        sample2     =   dedup.groupTuple().last()
        metaShards  =   shards.combine(metaId.groupTuple())
        TNSCOPE (   genomeShards,
                    metaShards,
                    sample1,
                    sample2,
                    bqsr_1,
                    bqsr_2,
                    params.genome_file)
        MERGE_VCF ( TNSCOPE.out.groupTuple()    )
    
    emit:
        vcf  = MERGE_VCF.out
        

}

workflow dnascope_tum_workflow {
    take:
        metaid
        dedup
    main:
        sample1     =   dedup.groupTuple().first()
        sample2     =   dedup.groupTuple().last()
        DNASCOPE_TUM (  params.genome_file, 
                        sample1,
                        sample2,
                        metaId.groupTuple() )
    emit:
        vcf = DNASCOPE_TUM.out
}

workflow dnascope_nor_workflow {
    take:
        metaid
        dedup
    main:
        sample1     =   dedup.groupTuple().first()
        sample2     =   dedup.groupTuple().last()
        DNASCOPE_NOR (  params.genome_file, 
                        sample1,
                        sample2,
                        metaId.groupTuple() )
    emit:
        vcf = DNASCOPE_NOR.out

}

workflow snv_calling_workflow {
    take:
        mode
        beds
        bam
        tnVcf
        fastq_sharded
        metaId

    main:
        bed             =    beds.combine(metaId.groupTuple())
        bamMerged       =    metaId.flatten().first().combine(bam).groupTuple()
        all             =   bamMerged.combine(bed)

        FREEBAYES ( mode, 
                    all, 
                    params.genome_file )
        VARDICT ( mode, 
                    all, 
                    params.genome_file )

        vcfPartFreebayes    =   FREEBAYES.out.groupTuple(by:[0,1])
        vcfPartVardict      =   VARDICT.out.groupTuple(by:[0,1])    
        vcfConcat           =   vcfPartFreebayes.mix(vcfPartVardict)
        CONCATENATE_VCFS (  vcfConcat,
                            params.genome_file  )
        conVcfs             =   CONCATENATE_VCFS.out.groupTuple()
        idMeta              =   metaId.groupTuple()
        AGGREGATE_VCFS (    mode,
                            conVcfs,
                            idMeta  )
        inputPonVcf         =   AGGREGATE_VCFS.out
        PON_FILTER (    params.PON_freebayes,
                        params.PON_vardict,
                        inputPonVcf )
        groupFasta          =   fastq_sharded.groupTuple() 
        mergVcf             =   tnVcf.groupTuple()
        
        GVCF_COMBINE (  mode,
                        groupFasta,
                        mergVcf )

        SPLIT_NORMALIZE (   params.genome_file, 
                            GVCF_COMBINE.out )

        INTERSECT ( params.intersect_bed, 
                    SPLIT_NORMALIZE.out )

        ANNOTATE_VEP (  params.CADD,
                        params.VEP_FASTA,
                        params.VEP_CACHE,
                        params.GNOMAD,
                        PON_FILTER.out
                    )
        def snv_panel = params.panels[params.panel]['PANEL_SNV']
        FILTER_WITH_PANEL_SNV ( snv_panel,
                                ANNOTATE_VEP.out )
    emit:
        ann_vcf = ANNOTATE_VEP.out
        vcf     = FILTER_WITH_PANEL_SNV.out
        
}

workflow sv_calling_workflow {
    take:
        mode
        cram
        metaId
    main:
        idManta                 =   metaId.groupTuple()
        fileManta               =   metaId.flatten().first().combine(cram).groupTuple()
        MANTA ( mode,
                params.genome_file, 
                idManta,
                fileManta   )
        ANNOTATE_MANTA (    params.SNPEFF_DIR,
                            MANTA.out )

        def fus_panel = params.panels[params.panel]['PANEL_FUS']    
        FILTER_WITH_PANEL_FUSIONS ( fus_panel,
                                    ANNOTATE_MANTA.out  )
    emit:
        vcf = FILTER_WITH_PANEL_FUSIONS.out

}

workflow cnv_calling_workflow {
    take:
        gatkid
        cram
     
    main:
        gatkBaf = metaId.flatten().first().combine(cram).join(gatkid, by:1)

        GATKCOV_BAF (   params.GATK_GNOMAD ,
                        params.genome_file ,
                        gatkBaf )
        gatkCovCount =  gatkBaf.groupTuple(by:1)

        GATKCOV_COUNT_TUM ( params.COV_INTERVAL_LIST,
                            params.GATK_PON_FEMALE,
                            params.GATK_PON_MALE,
                            params.GENOMEDICT,
                            params.genome_file,
                            gatkCovCount    )
        gatkcovcall =  GATKCOV_BAF.out.join(GATKCOV_COUNT_TUM.out[0],  by:1, remainder:true).groupTuple(by:1)
        GATKCOV_CALL_TUM (  params.GENOMEDICT,
                            gatkcovcall  )

        GATKCOV_COUNT_NOR ( params.COV_INTERVAL_LIST,
                            params.GATK_PON_FEMALE,
                            params.GATK_PON_MALE,
                            params.GENOMEDICT,
                            params.genome_file,
                            gatkCovCount    )
        gatkCovCallN =  GATKCOV_BAF.out.join(GATKCOV_COUNT_NOR.out[0],  by:1, remainder:true).groupTuple(by:1)
        GATKCOV_CALL_NOR (  params.GENOMEDICT,
                            gatkCovCallN   )
                            
        CNVS_ANNOTATE ( params.GENE_BED_PC,
                        GATKCOV_CALL_TUM.out[1]  )

        def cnv_panel = params.panels[params.panel]['PANEL_CNV']
        FILTER_WITH_PANEL_CNVS (    cnv_panel, 
                                    CNVS_ANNOTATE.out   )
    
    emit:
        tumplot =   GATKCOV_CALL_TUM.out[0]
        norplot =   GATKCOV_CALL_NOR.out[0]
        tumCov  =   GATKCOV_COUNT_TUM.out[1]
        norCov  =   GATKCOV_COUNT_NOR.out[1]
        vcf =   FILTER_WITH_PANEL_CNVS.out
}

workflow gens_tum_workflow {
    take:
        gatkcovTum
        dnascopeTum

    main:
        gens = dnascopeTum.join(gatkcovTum)
        GENERATE_GENS_DATA ( params.GENS_GNOMAD, gens)
}

workflow gens_nor_workflow {
    take:
        gatkcovNor
        dnascopeNor

    main:
        gensN = dnascopeNor.join(gatkcovNor)
        GENERATE_GENS_DATA_NOR ( params.GENS_GNOMAD, gensN)
}

workflow coyote_workflow {
    take:
        metaCoyote
        snvVcf
        svVcf
        cnvVcf
        tumplot
    main:
        filteredData = snvVcf.join(cnvVcf).join(svVcf)

    COYOTE (    filteredData,
                metaCoyote.groupTuple(),
                tumplot)
}

/*  Main Channel and workflow    */

/* Channels */

bwa_num_shards = params.bwa_shards
genomeShards   = params.genomic_shards_num
K_size = 100000000
mode =  file(params.csv).countLines() > 2 ? "paired" : "unpaired"

Channel
    .from( 0..bwa_num_shards-1)
    .set { bwa_shards }

Channel
    .fromPath(params.csv)
    .splitCsv(header:true)
    .map{ row-> tuple(row.group, row.id, file(row.read1), file(row.read2)) }
    .set { fastq_sharded }

Channel
    .fromPath(params.genomic_shards_file)
    .splitCsv(header:false)
    .set { shards }
    
Channel
    .fromPath(params.csv)
    .splitCsv(header:true)
    .map{ row-> tuple(row.id, row.diagnosis, file(row.read1)) }
    .set{ cron } 

Channel
    .fromPath(params.csv)
    .splitCsv(header:true)
    .map{ row-> tuple(row.group, row.id, row.type) }
    .set { metaId }

Channel
    .fromPath("${params.intersect_bed}")
    .ifEmpty { exit 1, "Regions bed file not found: ${params.intersect_bed}" }
    .splitText( by: 20000, file: 'bedpart.bed' )
    .map { it -> [it.getBaseName().tokenize(".")[1],it] }
    .set { beds }

Channel
    .fromPath(params.csv)
    .splitCsv(header:true)
    .map{ row-> tuple(row.group, row.id, row.sex, row.type) }
    .set { gatkId }

Channel
    .fromPath(params.csv)
    .splitCsv(header:true)
    .map{ row-> tuple(row.group, row.id, row.type, row.clarity_sample_id, row.clarity_pool_id) }
    .set { metaCoyote }


/*  Main workflow  */
workflow {

    sentieon_workflow (  K_size,
                        genomeShards,
                        bwa_num_shards, 
                        bwa_shards, 
                        fastq_sharded)
    cdm_qc_workflow (   cron, 
                        sentieon_workflow.out.sentieonqc )
    tnscope_workflow (  genomeShards,
                        metaId,
                        sentieon_workflow.out.bqsr,
                        sentieon_workflow.out.dedup,
                        shards )
    dnascope_tum_workflow ( metaId,
                            sentieon_workflow.out.dedup )
    dnascope_nor_workflow ( metaId,
                            sentieon_workflow.out.dedup )
    snv_calling_workflow (  mode,
                            beds,
                            sentieon_workflow.out.bam,
                            tnscope_workflow.out.vcf,
                            fastq_sharded,
                            metaId  )
     sv_calling_workflow (   mode,
                            sentieon_workflow.out.cram,
                            metaId)
    cnv_calling_workflow (  gatkId,
                            sentieon_workflow.out.cram    )

    gens_tum_workflow (cnv_calling_workflow.out.tumCov,        dnascope_tum_workflow.out.vcf )                   

    gens_nor_workflow ( cnv_calling_workflow.out.norCov,            dnascope_nor_workflow.out.vcf  ) 

    coyote_workflow (   metaCoyote,
                        snv_calling_workflow.out.vcf,
                        sv_calling_workflow.out.vcf,
                        cnv_calling_workflow.out.vcf,
                        cnv_calling_workflow.out.tumplot )    

    
}

/* Report the workflow */
workflow.onComplete {

	def msg = """\
		Pipeline execution summary
		---------------------------
		Completed at: ${workflow.complete}
		Duration    : ${workflow.duration}
		Success     : ${workflow.success}
		scriptFile  : ${workflow.scriptFile}
		workDir     : ${workflow.workDir}
		exit status : ${workflow.exitStatus}
		errorMessage: ${workflow.errorMessage}
		errorReport :
		"""
		.stripIndent()
	def error = """\
		${workflow.errorReport}
		"""
		.stripIndent()

	base = file(params.csv).getBaseName()
	logFile = file("$params.outdir/cron/logs/" + base + ".complete")
	logFile.text = msg
	logFile.append(error)
}