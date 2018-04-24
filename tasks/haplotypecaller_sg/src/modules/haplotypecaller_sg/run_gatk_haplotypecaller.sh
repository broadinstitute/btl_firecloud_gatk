#!/bin/bash -l

#. /broad/tools/scripts/useuse

#unuse Java-1.6
#use Java-1.7


libdir=$1
bamfile=$2
sampleID=$3
jobspecmemory=$4
dbsnp=$5
reference=$6
interval_file=$7
counter=$8


#command: sh <libdir>run_gatk_haplotypecaller.sh <libdir> ${text bam.file} ${text sampleID} ${text job.spec.memory} ${text dbsnp} ${text reference} ${text interval_file} ${scatter}

#Get scatter_interval:
scatter_interval=$(head -n ${counter} ${interval_file} | tail -n 1)

java -Xmx${jobspecmemory}g -jar $libdir/GenomeAnalysisTK-3.3-0-g37228af/GenomeAnalysisTK.jar \
	-T HaplotypeCaller \
	-R $reference \
	-I $bamfile \
	--max_alternate_alleles 3 \
	--dbsnp $dbsnp \
	--genotyping_mode DISCOVERY \
	--variant_index_type LINEAR \
	-ERC GVCF \
	--variant_index_parameter 128000 \
	--minPruning 2 \
	-stand_call_conf 30.0 \
	-stand_emit_conf 30.0 \
	-A DepthPerSampleHC \
	-A StrandBiasBySample \
	-pairHMM VECTOR_LOGLESS_CACHING \
	-o ${sampleID}_scatter${counter}.gvcf.gz \
    --disable_auto_index_creation_and_locking_when_reading_rods \
	${scatter_interval}

