#!/bin/bash -l

#. /broad/tools/scripts/useuse

#unuse Java-1.6
#use Java-1.7
#use Tabix
#use .vcftools-0.1.10


libdir=$1
sampleID=$2


#dir=$(pwd)
dir=.
find -L $dir -name '*.gvcf.gz' | sort -t "_" -k 2 -g > ${dir}/gvcfs_final.list
files=$(cat ${dir}/gvcfs_final.list)
/opt/src/modules/haplotypecaller_sg/vcftools/bin/vcf-concat $files | /opt/src/modules/tabix/bgzip -c > ${dir}/${sampleID}.gvcf.gz  &&   /opt/src/modules/tabix/tabix -p vcf ${dir}/${sampleID}.gvcf.gz

