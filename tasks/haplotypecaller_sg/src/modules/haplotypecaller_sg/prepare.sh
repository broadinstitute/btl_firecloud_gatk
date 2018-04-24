#!/bin/bash -l

#. /broad/tools/scripts/useuse

#unuse Java-1.6
#use Java-1.8 > /dev/null

#export JAVA_HOME=/opt/java8/jdk1.8.0_31

libdir=$1
reference_sequence=$2
number_scatter_jobs=$3
jobspecmemory=$4
bamfile=$5
sampleID=$6
dbsnp_file=$7

#<libdir> ${text reference_sequence} ${text number_scatter_jobs} ${text jobspecmemory} ${bamfile} ${sampleID} ${dbsnp_file}

#Create list of even intervals for Scatter-Gather jobs
$JAVA_HOME/jre/bin/java -jar ${libdir}/hellbender-all-1.0.jar CreateEvenIntervals -R $reference_sequence -n $number_scatter_jobs --output_intervals even_intervals.txt > /dev/null


#Print parameters for Scatter-Gather jobs
current_dir=$(pwd)
for scatter in $(eval echo {1..${number_scatter_jobs}})
do
echo "$libdir $bamfile $sampleID $jobspecmemory $dbsnp_file $reference_sequence ${current_dir}/even_intervals.txt $scatter"
#command: sh <libdir>run_gatk_haplotypecaller.sh <libdir> ${text bam.file} ${text sampleID} ${text job.spec.memory} ${text dbsnp} ${text reference} ${text interval_file} ${scatter}
done


#Print parameters for Gather step
echo "$libdir $sampleID"
