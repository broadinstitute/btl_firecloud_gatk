#!/bin/bash
LIBDIR=$1;
shift;
#chmod u+x $LIBDIR/samtools;

out=`echo $1 | tr '/' '\n' | tail -1`

#$LIBDIR/samtools view -H $1 > $out.header.txt;
$LIBDIR/samtools view -H $1 | grep -P "RG\t" | sed 's/\t/\n/g' | grep SM: | sort | uniq | head -1 | sed 's/SM://' > upload.txt
