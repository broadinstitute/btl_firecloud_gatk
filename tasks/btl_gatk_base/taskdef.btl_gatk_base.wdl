# GATK WDL
# import "hc_scatter.wdl" as sub

task VersionCheck {
    String gatk
    Float version
    command {
        set -euo pipefail
        /opt/src/algutil/monitor_start.py
        python /opt/src/version_check.py -i "java -jar ${gatk} --version" -v ${version}
        /opt/src/algutil/monitor_stop.py
    } output {
        String out = "done"
        File monitor_start="monitor_start.log"
        File monitor_stop="monitor_stop.log"
        File dstat="dstat.log"

    } runtime {
        task_name: "VersionCheck"
        docker : "gcr.io/btl-dockers/btl_gatk:1"

    }
} 

# TODO remove this
#task **MakeOutputDir** {
#    String go
#    String output_dir
#    command {
#        mkdir -p ${output_dir}
#    }
#    output {
#        String out = output_dir
#    } runtime {
#        task_name: "**MakeOutputDir**"
#    }
#    parameter_meta {
#        output_dir: "The root directory for where all sample directories and ref index files will be deposited."
#    }
#}


#TODO not called anywhere, do we remove it? Or, use this code to conditionally generate an index.
# note uses only a single index file to determine if all files are good
#task CheckIndex {
#    File ref_file
#    String out_file = sub(ref_file, "\\.fa*$", ".dict")
#    String seq_dict = "${out_file}"
#    command {
#        set -euo pipefail
#        /opt/src/algutil/monitor_start.py
#        python -c "import os; print ( '0' if os.path.isfile('${seq_dict}') else '1')"
#        /opt/src/algutil/monitor_stop.py
#    }
#    output {
#        Int retcode = read_int(stdout())
#        File out = "${ref_file}"
#        File sd = seq_dict
#        File monitor_start="monitor_start.log"
#        File monitor_stop="monitor_stop.log"
#        File dstat="dstat.log"
#    } runtime {
#        task_name: "CheckIndex"
#        docker : "gcr.io/btl-dockers/btl_gatk:1"
#
#    }
#    parameter_meta {
#        output_dir: "The root directory for where all sample directories and ref index files will be deposited."
#        ref_file: "The reference fasta file name (without the path)."
#        out_file: "The seq dict filename must be generated and is used to construct seq_dict"
#        seq_dict: "The sequence dictionary created by Picard CreateSequenceDictionary"
#    }
#}

task copyFile{
    File source
    String dest
    command {
        echo "copy ${source} to the output"
    }
    output {
        File out = source
    }
}

#task copyFile {
#    File source
#    String dest
#    command {
#        cp ${source} ${dest}
#    }
#}

task IndexReference {
    File ref_fasta
    String picard

    #String out_file = sub(ref_fasta, "\\.fasta*$", ".dict")
    String out_file = "ref.dict"
    command {
        set -euo pipefail
        ln -sT `pwd` /opt/execution
        ln -sT `pwd`/../inputs /opt/inputs

        /opt/src/algutil/monitor_start.py
        ln ${ref_fasta} ref.fasta

        bwa index ref.fasta
        samtools faidx ref.fasta
        java -jar ${picard} CreateSequenceDictionary REFERENCE=ref.fasta O=ref.dict
        ls 
        /opt/src/algutil/monitor_stop.py

    }
    output {
        File out = "ref.fasta"
        File dict = "ref.dict"
        File amb = "ref.fasta.amb"
        File ann = "ref.fasta.ann"
        File bwt = "ref.fasta.bwt"
        File fai = "ref.fasta.fai"
        File pac = "ref.fasta.pac"
        File sa = "ref.fasta.sa"
        File monitor_start="monitor_start.log"
        File monitor_stop="monitor_stop.log"
        File dstat="dstat.log"
    } runtime {
        task_name: "IndexReference"
        docker : "gcr.io/btl-dockers/btl_gatk:1"
        
    }
    parameter_meta {
        ref_fasta: "The name of the reference file (without the path)."
        picard: "The path to the picard executable jar file."
        output_dir: "The root directory for where all sample directories and ref index files will be deposited."
        out_file: "The seq dict filename must be generated and is used to construct seq_dict"
        seq_dict: "The sequence dictionary created by Picard CreateSequenceDictionary"
        old_ref: "The full path to the original reference file location."
        ref: "The absolute path of the reference file to be used by the workflow."
    }
}



task removeFile {
    File file
    String ready
    command {
        echo this task has been disabled
    } runtime {
        task_name: "removeFile"
    }
}

task CreateIntervalsList {
    File ref
    Float interval_size
    command {
        set -euo pipefail
        ln -sT `pwd` /opt/execution
        ln -sT `pwd`/../inputs /opt/inputs
        /opt/src/algutil/monitor_start.py
        python /opt/src/intervals_creator.py -r ${ref} \
        -i ${interval_size} > intervals.list
        /opt/src/algutil/monitor_stop.py

    }
    output{
        File out = "intervals.list"
        File monitor_start="monitor_start.log"
        File monitor_stop="monitor_stop.log"
        File dstat="dstat.log"
    } runtime {
        task_name: "CreateIntervalList"
        docker : "gcr.io/btl-dockers/btl_gatk:1"

    }
    parameter_meta {
        ref: "The absolute path of the reference file to be used by the workflow."
        interval_size: "The size in gigabases that each interval should approximate."
        output_dir: "The root directory for where all sample directories and ref index files will be deposited."
    }
}

# TODO delete this, it is no longer called.
#task **SymlinkCromwellExecutionDir** {
#    String output_dir
#
#    command {
#        ln -s `dirname \`dirname $PWD\`` ${output_dir}
#    }
#}

# TODO remove??
#task GenerateFastqNames {
#    String sample_name
#    String ready
#    command {
#        set -euo pipefail
#        /opt/src/algutil/monitor_start.py
#        echo "GenerateFastqNames: ${sample_name}.1.fq, ${sample_name}.2.fq"
#        /opt/src/algutil/monitor_stop.py
#
#    }
#    output {
#        Array[File] fastq_out = ["${sample_name}.1.fq", "${sample_name}.2.fq"]
#        File monitor_start="monitor_start.log"
#        File monitor_stop="monitor_stop.log"
#        File dstat="dstat.log"
#    } runtime {
#        task_name: "GenerateFastqNames"
#        docker : "gcr.io/btl-dockers/btl_gatk:1"
#
#    }
#    parameter_meta {
#        sample_name: "The name of the sample as indicated by the 1st column of the gatk.samples_file json input."
#    }
#}

task SamToFastq {
    String picard
    File in_bam
    String sample_name
    String out_fq1 = "${sample_name}.1.fq"
    String out_fq2 = "${sample_name}.2.fq"
    command {
        set -euo pipefail
        /opt/src/algutil/monitor_start.py
        java -Xmx12G -jar ${picard} SamToFastq INPUT=${in_bam} FASTQ=${out_fq1} SECOND_END_FASTQ=${out_fq2} VALIDATION_STRINGENCY=LENIENT
        /opt/src/algutil/monitor_stop.py

    }
    output {
        String done = "Done"
        File fq1 = out_fq1
        File fq2 = out_fq2
        File monitor_start="monitor_start.log"
        File monitor_stop="monitor_stop.log"
        File dstat="dstat.log"
    } runtime {
        task_name: "SamToFastq"
        docker : "gcr.io/btl-dockers/btl_gatk:1"
    }
    parameter_meta {
        picard: "The absolute path to the picard jar to execute."
        in_bam: "The bam file to convert to fastq."
        sample_dir: "The sample-specific directory inside output_dir for each sample."
        sample_name: "The name of the sample as indicated by the 1st column of the gatk.samples_file json input."
        out_fq1: "The fastq file containing the first read of each pair."
        out_fq2: "The fastq file containing the second read of each pair"
    }
}

# TODO - change this to be a passthrough for a fastq file

#task **CopyFastq** {
#    File fq
#    Int pair
#    String sample_name
#    File out_fq = "${sample_name}.${pair}.fq"
#    command {
#        cp ${fq} ${out_fq}
#    } output {
#    String out = out_fq
#    } runtime {
#        task_name: "**CopyFastq**"
#    }
#    parameter_meta {
#        fq: "The fastq file to copy."
#        sample_dir: "The sample-specific directory inside output_dir for each sample."
#        sample_name: "The name of the sample as indicated by the 1st column of the gatk.samples_file json input."
#        out_fq: "the absolute path of the copied fastq file."
#    }
#}

task AlignBAM {
    String sample_name
    File ref
    File dict
    File amb
    File ann
    File bwt
    File fai
    File pac
    File sa
    File fq1
    File fq2
    String read_group = "'@RG\\tID:FLOWCELL_${sample_name}\\tSM:${sample_name}\\tPL:ILLUMINA\\tLB:LIB_${sample_name}'"
    command {
        set -euo pipefail
        /opt/src/algutil/monitor_start.py
        bwa mem -t 8 -R ${read_group} ${ref} ${fq1} ${fq2} | samtools view -bS -> ${sample_name}.aligned.bam
        /opt/src/algutil/monitor_stop.py
    }
    output {
        File bam = "${sample_name}.aligned.bam"
        String done = "done"
        File monitor_start="monitor_start.log"
        File monitor_stop="monitor_stop.log"
        File dstat="dstat.log"
    } runtime {
        task_name: "AlignBAM"
        docker : "gcr.io/btl-dockers/btl_gatk:1"
    }
    parameter_meta {
        ref: "fasta file of reference genome"
        sample_dir: "The sample-specific directory inside output_dir for each sample."
        sample_name: "The name of the sample as indicated by the 1st column of the gatk.samples_file json input."
        fq_array: "An array containing the paths to the first and second fastq files."
        read_group: "The read group string that will be included in the bam header."
    }
}

task SortBAM {
    String picard
    String sample_name
    File aligned_bam
    command {
        set -euo pipefail
        /opt/src/algutil/monitor_start.py
        java -Xmx8G -jar ${picard} SortSam I=${aligned_bam} O=${sample_name}.sorted.bam SO=coordinate
        /opt/src/algutil/monitor_stop.py

    }
    output {
        File bam = "${sample_name}.sorted.bam"
        String done = "Done"
        File monitor_start="monitor_start.log"
        File monitor_stop="monitor_stop.log"
        File dstat="dstat.log"
    } runtime {
        task_name: "SortBAM"
        docker : "gcr.io/btl-dockers/btl_gatk:1"
    }
    parameter_meta {
        picard: "The absolute path to the picard jar to execute."
        aligned_sam: "The sam file to be sorted."
        sample_dir: "The sample-specific directory inside output_dir for each sample."
        bam: "The sorted bam file."
    }
}

task MarkDuplicates {
    String picard
    String sample_name
    File sorted_bam
    command {
        set -euo pipefail
        /opt/src/algutil/monitor_start.py
        java -Xmx8G -jar ${picard} MarkDuplicates I=${sorted_bam} O=${sample_name}.marked_duplicates.bam M=${sample_name}.marked_duplicates.metrics
        /opt/src/algutil/monitor_stop.py
    }
    output {
        File bam = "${sample_name}.marked_duplicates.bam"
        String done = "done"
        File monitor_start="monitor_start.log"
        File monitor_stop="monitor_stop.log"
        File dstat="dstat.log"
    } runtime {
        task_name: "MarkDuplicates"
        docker : "gcr.io/btl-dockers/btl_gatk:1"
    }
    parameter_meta {
        picard: "The absolute path to the picard jar to execute."
        sorted_bam: "The sorted bam file to mark duplicates on."
        out_bam: "The bam file where duplicates are marked."
        metrics: "The marked duplicates metrics file."
    }
}

task ReorderBAM {
    String picard
    String sample_name
    File marked_bam
    String out_bam = "${sample_name}.reordered.bam"
    File ref
    File dict
    command {
        set -euo pipefail
        /opt/src/algutil/monitor_start.py
        java -Xmx8G -jar ${picard} ReorderSam I=${marked_bam} O=${out_bam} R=${ref}
        samtools index ${out_bam}
        /opt/src/algutil/monitor_stop.py

    }
    output {
        File bam = out_bam
        File index = "${sample_name}.reordered.bam.bai"
        String done = "done"
        File monitor_start="monitor_start.log"
        File monitor_stop="monitor_stop.log"
        File dstat="dstat.log"
    } runtime {
        task_name: "ReorderBAM"
        docker : "gcr.io/btl-dockers/btl_gatk:1"
    }
    parameter_meta {
        picard: "The absolute path to the picard jar to execute."
        marked_bam: "The bam file where duplicates are marked."
        out_bam: "The reordered bam file."
        ref: "fasta file of reference genome"
    }
}

task IndexBAM {
    File in_bam
    File index
    command {
        set -euo pipefail
        /opt/src/algutil/monitor_start.py
        #TODO what is the "rm index" doing?
        rm index
        samtools index ${in_bam}
        /opt/src/algutil/monitor_stop.py

    } 
    output {
        File monitor_start="monitor_start.log"
        File monitor_stop="monitor_stop.log"
        File dstat="dstat.log"
    } 
    
    runtime {
        task_name: "IndexBAM"
        docker : "gcr.io/btl-dockers/btl_gatk:1"
    }
    parameter_meta {
        in_bam: "the reordered bam file to be indexed."
    }
}

task RealignerTargetCreator {
    String gatk
    File ref
    File dict
    File amb
    File ann
    File bwt
    File fai
    File pac
    File sa
    File in_bam
    File index
    String sample_name
    String out = "${sample_name}.interval_list"
    command {
        set -euo pipefail
        /opt/src/algutil/monitor_start.py
        java -Xmx8G -jar ${gatk} -T RealignerTargetCreator -nct 1 -nt 24 -R ${ref} -I ${in_bam} -o ${out}
        /opt/src/algutil/monitor_stop.py

    }
    output {
        File intervals = out
        File monitor_start="monitor_start.log"
        File monitor_stop="monitor_stop.log"
        File dstat="dstat.log"
    } runtime {
        task_name: "RealignerTargetCreator"
        docker : "gcr.io/btl-dockers/btl_gatk:1"
    }
    parameter_meta {
        gatk: "The absolute path to the gatk executable jar."
        ref: "fasta file of reference genome"
        in_bam: "The input bam for the gatk task"
        out: "The intervals list to be used by IndelRealigner"
    }
}

task IndelRealigner {
    String gatk
    File ref
    File dict
    File amb
    File ann
    File bwt
    File fai
    File pac
    File sa
    File in_bam
    File intervals
    String sample_name
    String out = "${sample_name}.indels_realigned.bam"

    command {
        set -euo pipefail
        /opt/src/algutil/monitor_start.py
        samtools index ${in_bam}
        java -Xmx4G -jar ${gatk} -T IndelRealigner -nct 1 -nt 1 -R ${ref} -I ${in_bam} -targetIntervals ${intervals} -o ${out}
        samtools index ${out}
        /opt/src/algutil/monitor_stop.py

    }
    output {
        File bam = "${sample_name}.indels_realigned.bam"
        File index = "${sample_name}.indels_realigned.bam.bai"
        String done = "done"
        File monitor_start="monitor_start.log"
        File monitor_stop="monitor_stop.log"
        File dstat="dstat.log"
    } runtime {
        task_name: "IndelRealigner"
        docker : "gcr.io/btl-dockers/btl_gatk:1"
    }
    parameter_meta {
        gatk: "The absolute path to the gatk executable jar."
        ref: "fasta file of reference genome"
        in_bam: "The input bam for the gatk task"
        intervals: "The intervals list to be used by IndelRealigner"
        out: "the bam including realigned indels."
    }
}

#task CreateBamList {
#    File samples_file
#    String output_dir
#    File out_file = "${output_dir}bqsr_bams_list.txt"
#    command {
#        set -euo pipefail
#        /opt/src/algutil/monitor_start.py
#        cat ${samples_file} | cut -f 2 > ${out_file}
#        /opt/src/algutil/monitor_stop.py
#
#    }
#    output {
#        File out = out_file
#        File monitor_start="monitor_start.log"
#        File monitor_stop="monitor_stop.log"
#        File dstat="dstat.log"
#    } runtime {
#        task_name: "CreateBamList"
#        docker : "gcr.io/btl-dockers/btl_gatk:1"
#    }
#    parameter_meta {
#        samples_file: "The text file containing sample names and bam file paths."
#        output_dir: "The root directory for where all sample directories and ref index files will be deposited."
#        out_file: "The text file containing absolute paths to all bam files, one per line."
#    }
#}

task BaseRecalibrator {
    String gatk
    File ref
    File dict
    File amb
    File ann
    File bwt
    File fai
    File pac
    File sa
    File bam
    File index
    Array[String] known_sites
    String ? bqsr
    String out_file
    command {
        set -euo pipefail
        /opt/src/algutil/monitor_start.py
        samtools index ${bam}
        java -Xmx4G  -jar ${gatk} -T BaseRecalibrator -nct 8 -nt 1 -R ${ref} -I  ${bam} -knownSites ${sep=" -knownSites " known_sites} -o ${out_file} ${" -BQSR " + bqsr}
        /opt/src/algutil/monitor_stop.py
    }
    output {
        File out = "${out_file}"
        File monitor_start="monitor_start.log"
        File monitor_stop="monitor_stop.log"
        File dstat="dstat.log"
    } runtime {
        task_name: "BaseRecalibrator"
        docker : "gcr.io/btl-dockers/btl_gatk:1"
    }
    parameter_meta {
        gatk: "The absolute path to the gatk executable jar."
        ref: "fasta file of reference genome"
        bam_list: "The text file containing absolute paths to all bam files, one per line."
        known_sites: "An array of databases of known polymorphic sites."
        bqsr: "Full path to bqsr file."
    }
}

task AnalyzeCovariates {
    String gatk
    File ref
    File dict
    File amb
    File ann
    File bwt
    File fai
    File pac
    File sa
    String sample_name
    String before
    String after
    command {
        set -euo pipefail
        /opt/src/algutil/monitor_start.py

        #use R-3.1
        java -Xmx8G -jar ${gatk} \
        -T AnalyzeCovariates \
        -R ${ref} \
        -before ${before} \
        -after ${after} \
        -plots ${sample_name}.recalibration_plots.pdf
        /opt/src/algutil/monitor_start.py

    } output {
        String done = "done"
        File out = "${sample_name}.recalibration_plots.pdf"
        File monitor_start="monitor_start.log"
        File monitor_stop="monitor_stop.log"
        File dstat="dstat.log"
    } runtime {
        task_name: "AnalyzeCovariates"
        docker : "gcr.io/btl-dockers/btl_gatk:1"
    }
    parameter_meta {
        gatk: "The absolute path to the gatk executable jar."
        ref: "fasta file of reference genome"
        sample_dir: "The sample-specific directory inside output_dir for each sample."
        sample_name: "The name of the sample as indicated by the 1st column of the gatk.samples_file json input."
        known_sites: "An array of databases of known polymorphic sites."
        before: "The table output from the first BaseRecalibrator step."
        after: "The table output from the second BaseRecalibrator step."
    }
}

task PrintReads {
    String gatk
    String sample_name
    String ref
    File in_bam
    File dict
    File amb
    File ann
    File bwt
    File fai
    File pac
    File sa
    String bqsr
    String out =  "${sample_name}.bqsr.bam"
    command {
        set -euo pipefail
        /opt/src/algutil/monitor_start.py
        samtools index ${in_bam}
        java -Xmx4G -jar ${gatk} -T PrintReads -nct 8 -nt 1 -R ${ref} -I ${in_bam} -BQSR ${bqsr} -o ${out}
        samtools index ${out}
        /opt/src/algutil/monitor_stop.py

    }
    output {
        File bam = out
        # File index = "$[sample_name}.bqsr.bam.bai"
        String done = "done"
        File monitor_start="monitor_start.log"
        File monitor_stop="monitor_stop.log"
        File dstat="dstat.log"
    } runtime {
        task_name: "PrintReads"
        docker : "gcr.io/btl-dockers/btl_gatk:1"
    }
    parameter_meta {
        gatk: "The absolute path to the gatk executable jar."
        ref: "fasta file of reference genome"
        in_bam: "The input bam for PrintReads."
        bqsr: "Full path to bqsr file."
        out: "The bam file with bqsr applied to it."
    }
}
# NOTE: HaplotypeCaller complains if I don't provide --variant_index_type and variant_index_parameter
task HaplotypeCaller {
    String gatk
	File ref
    File dict
    File amb
    File ann
    File bwt
    File fai
    File pac
    File sa
	String sample_name
	File in_bam
	#File index
	File ? intervals
	File ? bqsr_file
	Int ? ploidy
	String ? erc
	String ? extra_hc_params
    String out = "${sample_name}.g.vcf"
	command {
        set -euo pipefail
        /opt/src/algutil/monitor_start.py
	    samtools index ${in_bam}
		java -Xmx8G -jar ${gatk} \
			-T HaplotypeCaller \
			-nt 1 \
			-R ${ref} \
			--input_file ${in_bam} \
			${"--intervals " + intervals} \
			${"-BQSR " + bqsr_file} \
			-ERC ${default="GVCF" erc} \
			-ploidy ${default="2" ploidy} \
			--interval_padding 100 \
			-o ${out} \
			-variant_index_type LINEAR \
			-variant_index_parameter 128000 \
            ${default="\n" extra_hc_params}
        /opt/src/algutil/monitor_stop.py

	}
	output {
		#To track additional outputs from your task, please manually add them below
		File vcf = out
		String done = "done"
        File monitor_start="monitor_start.log"
        File monitor_stop="monitor_stop.log"
        File dstat="dstat.log"
    } runtime {
        task_name: "HaplotypeCaller"
        docker : "gcr.io/btl-dockers/btl_gatk:1"
    }
	parameter_meta {
		gatk: "Executable jar for the GenomeAnalysisTK"
		ref: "fasta file of reference genome"
        sample_name: "The name of the sample as indicated by the 1st column of the gatk.samples_file json input."
        sample_dir: "The sample-specific directory inside output_dir for each sample."
        in_bam: "The bam file to call HaplotypeCaller on."
        intervals: "An array of intervals to restrict processing to."
        bqsr_file: "The full path to the BQSR file."
        erc: "Mode for emitting reference confidence scores."
        extra_hc_params: "A parameter that allows users to pass any additional paramters to the task."
        out: "VCF file produced by haplotype caller."
	}
}

task GenotypeGVCFs {
	String gatk
	File ref
    File dict
    File amb
    File ann
    File bwt
    File fai
    File pac
    File sa
	String ? extra_gg_params # If a parameter you'd like to use is missing from this task, use this term to add your own string
	File ? intervals
	Array[File] variant_files
    Boolean ? all_sites
    String gcvf_out = "genotypeGVCFs.vcf"
	command {
        set -euo pipefail
        /opt/src/algutil/monitor_start.py
		java -Xmx8G -jar ${gatk} \
			-T GenotypeGVCFs \
			-R ${ref} \
			${sep=" --intervals " "--intervals " + intervals} \
			-o ${gcvf_out} \
			-V ${sep=" -V " variant_files} \
			${true="-allSites" false="" all_sites} \
			${default="\n" extra_gg_params}
        /opt/src/algutil/monitor_stop.py

	}
	output {
		#To track additional outputs from your task, please manually add them below
		File out = gcvf_out
		String done = "done"
        File monitor_start="monitor_start.log"
        File monitor_stop="monitor_stop.log"
        File dstat="dstat.log"
    } runtime {
        task_name: "GenotypeGVCFs"
        docker : "gcr.io/btl-dockers/btl_gatk:1"
    }
	parameter_meta {
		gatk: "Executable jar for the GenomeAnalysisTK"
		ref: "fasta file of reference genome"
		extra_gg_params: "An optional parameter which allows the user to specify additions to the command line at run time"
		out: "File to which variants should be written"
		ploidy: "Ploidy (number of chromosomes) per sample. For pooled data, set to (Number of samples in each pool * Sample Ploidy)."
		variant_files: "One or more input gVCF files"
		intervals: "One or more genomic intervals over which to operate"
		all_sites: "Include loci found to be non-variant after genotyping"
		gcvf_out: "The output vcf of GenotypeGVCFs"
	}
}

# https://software.broadinstitute.org/gatk/documentation/article.php?id=1259
task VariantRecalibrator {
	String gatk
	File ref
    File dict
    File amb
    File ann
    File bwt
    File fai
    File pac
    File sa
	String mode
	File ? intervals
	File task_input
    Array[String] resource
    Array[String] annotation
    Int ? max_gaussians
    Int ? mq_cap
	String tranches_file = "${mode}.tranches"
    String recal_file = "${mode}.recal"
    String rscript_file = "${mode}.plots.R"
	String ? extra_vr_params # If a parameter you'd like to use is missing from this task, use this term to add your own string
	command {
        set -euo pipefail
        /opt/src/algutil/monitor_start.py
		java -Xmx8G -jar ${gatk} \
			-T VariantRecalibrator \
			-R ${ref} \
			${default="" "--intervals " + intervals} \
			-input ${task_input} \
			-mode ${mode} \
            -resource:${sep=" -resource:" resource} \
			-recalFile ${recal_file} \
			-tranchesFile ${tranches_file} \
			-rscriptFile ${rscript_file} \
			-an ${sep=" -an " annotation} \
			--maxGaussians ${max_gaussians} \
			--MQCapForLogitJitterTransform ${mq_cap}
			${default="\n" extra_vr_params}
        /opt/src/algutil/monitor_stop.py

	}

	output {
		#To track additional outputs from your task, please manually add them below
		File out = task_input
		File tranches = tranches_file
		File recal = recal_file
		File rscript = rscript_file
		String done = "done"
        File monitor_start="monitor_start.log"
        File monitor_stop="monitor_stop.log"
        File dstat="dstat.log"
    } runtime {
        task_name: "VariantRecalibrator"
        docker : "gcr.io/btl-dockers/btl_gatk:1"
    }
	parameter_meta {
		gatk: "Executable jar for the GenomeAnalysisTK"
		ref: "fasta file of reference genome"
		extra_vr_params: "An optional parameter which allows the user to specify additions to the command line at run time"
		aggregate: "Additional raw input variants to be used in building the model"
		task_input: "One or more VCFs of raw input variants to be recalibrated"
		recal_file: "The output recal file used by ApplyRecalibration"
		resource: "A list of sites for which to apply a prior probability of being correct but which aren't used by the algorithm (training and truth sets are required to run)"
		tranches_file: "The output tranches file used by ApplyRecalibration"
		intervals: "One or more genomic intervals over which to operate"
		mode: "The mode for recalibration (indel or snp)."
		annotation: "An array of annotations to use for calculations."
		max_gaussians: "Max number of Gaussians for the positive model"
		mq_cap: "Apply logit transform and jitter to MQ values"
	}
}

task ApplyRecalibration {
    String gatk
    File ref
    File dict
    File amb
    File ann
    File bwt
    File fai
    File pac
    File sa
    File vcf_in
    Float ts_filter
    File recal_file
    File tranches
    String mode
    String prefix
    String vcf_out = "${prefix}.recalibrated.filtered.vcf"
    command {
        set -euo pipefail
        /opt/src/algutil/monitor_start.py
		java -Xmx8G -jar ${gatk} \
            -T ApplyRecalibration \
            -R ${ref} \
            -input ${vcf_in} \
            --ts_filter_level ${ts_filter} \
            -tranchesFile ${tranches} \
            -recalFile ${recal_file} \
            -mode ${mode} \
            -o ${prefix}.recalibrated.filtered.vcf
        /opt/src/algutil/monitor_stop.py

		}
    output {
        File out = "${prefix}.recalibrated.filtered.vcf"
        String done = "done"
        File monitor_start="monitor_start.log"
        File monitor_stop="monitor_stop.log"
        File dstat="dstat.log"
    } runtime {
        task_name: "ApplyRecalibration"
        docker : "gcr.io/btl-dockers/btl_gatk:1"
    }
    parameter_meta {
		gatk: "Executable jar for the GenomeAnalysisTK"
		ref: "fasta file of reference genome"
		vcf_in: "The raw input variants to be recalibrated."
		ts_filter: "The truth sensitivity level at which to start filtering"
		recal_file: "The output recal file used by ApplyRecalibration"
		mode: "Recalibration mode to employ: 1.) SNP for recalibrating only SNPs (emitting indels untouched in the output VCF); 2.) INDEL for indels; and 3.) BOTH for recalibrating both SNPs and indels simultaneously."
        vcf_out: "The output filtered and recalibrated VCF file in which each variant is annotated with its VQSLOD value"
    }
}

# http://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set
task SelectVariants {
    String gatk
    File ref
    File dict
    File amb
    File ann
    File bwt
    File fai
    File pac
    File sa
    File vcf_in
    String mode
    command {
        set -euo pipefail
        /opt/src/algutil/monitor_start.py
        java -Xmx8G -jar ${gatk} \
            -T SelectVariants \
            -R ${ref} \
            -V ${vcf_in} \
            -selectType ${mode} \
            -o select${mode}.vcf
        /opt/src/algutil/monitor_stop.py

    }
    output {
        File out = "select${mode}.vcf"
        String done = "Done"
        File monitor_start="monitor_start.log"
        File monitor_stop="monitor_stop.log"
        File dstat="dstat.log"
    } runtime {
        task_name: "SelectVariants"
        docker : "gcr.io/btl-dockers/btl_gatk:1"
    }
    parameter_meta {
        gatk: "Executable jar for the GenomeAnalysisTK"
        ref: "fasta file of reference genome"
        vcf_in: "The input variants file."
        vcf_out: "The output variants file."
        }
}

task HardFiltration {
    String gatk
    File ref
    File dict
    File amb
    File ann
    File bwt
    File fai
    File pac
    File sa
    File vcf_in
    String variant_type
    String vcf_out = "filtered_${variant_type}.vcf"
    String filter_expression
    command {
        set -euo pipefail
        /opt/src/algutil/monitor_start.py
        java -Xmx8G -jar ${gatk} \
            -T VariantFiltration \
            -R ${ref} \
            -V ${vcf_in} \
            --filterExpression "${filter_expression}" \
            --filterName "my_variant_filter" \
            -o ${vcf_out}
        /opt/src/algutil/monitor_stop.py

    }
    output {
        File out = vcf_out
        String done = "done"
        File monitor_start="monitor_start.log"
        File monitor_stop="monitor_stop.log"
        File dstat="dstat.log"
    } runtime {
        task_name: "HardFiltration"
        docker : "gcr.io/btl-dockers/btl_gatk:1"
    }
    parameter_meta {
        gatk: "Executable jar for the GenomeAnalysisTK"
        ref: "fasta file of reference genome"
        vcf_in: "The input variants file."
        vcf_out: "The output variants file."
        filter_expression: "The user-defined expressions to indicate which variants to filter."
    }
}

task CombineVariants {
    String gatk
    File ref
    File dict
    File amb
    File ann
    File bwt
    File fai
    File pac
    File sa
    File vcf1
    File vcf2
    String outfile = "filtered.combined.vcf"
    command {
        set -euo pipefail
        /opt/src/algutil/monitor_start.py
        java -jar -Xmx8G ${gatk} \
            -T CombineVariants \
            -R ${ref} \
            --variant ${vcf1} \
            --variant ${vcf2} \
            -o ${outfile} \
            -genotypeMergeOptions UNIQUIFY
        /opt/src/algutil/monitor_stop.py

    } runtime {
        task_name: "CombineVariants"
        docker : "gcr.io/btl-dockers/btl_gatk:1"
    }
    output {
        File out = outfile
        File monitor_start="monitor_start.log"
        File monitor_stop="monitor_stop.log"
        File dstat="dstat.log"
    }
}

# Based on http://gatkforums.broadinstitute.org/gatk/discussion/50/adding-genomic-annotations-using-snpeff-and-variantannotator
task SnpEff {
    File vcf_in
    String vcf_out = "SnpEff.annotated.vcf"
    String snpeff
    String snpeff_db
    String ? snpeff_extra_params
    command {
        set -euo pipefail
        /opt/src/algutil/monitor_start.py
        cp /cil/shed/apps/external/snpEff/snpEff-4.1g/snpEff.config .
        java -Xmx4G -jar ${snpeff} -formatEff -no-downstream -no-intergenic -no-upstream -no-utr \
            ${snpeff_db} ${vcf_in} ${snpeff_extra_params} > ${vcf_out}
        /opt/src/algutil/monitor_stop.py

    }
    runtime {
        task_name: "CombineVariants"
        docker : "gcr.io/btl-dockers/btl_gatk:1"
    }
    output {
        File out = vcf_out
        File monitor_start="monitor_start.log"
        File monitor_stop="monitor_stop.log"
        File dstat="dstat.log"
    }
    parameter_meta {
        vcf_in: "The input variants file."
        vcf_out: "The output variants file."
        snp_eff: "The path to the snpEff executable."
        snpeff_db: "The snpeff database to use."
        snpeff_extra_params: "A string for passing any additional parameters to snpeff."
    }
}

workflow gatk {
    # Initialize workflow
    # Global parameters
    File samples_file
    File ref_fasta
    String picard = "/cil/shed/apps/external/picard/current/bin/picard.jar"
    String gatk = "/humgen/gsa-hpprojects/GATK/bin/GenomeAnalysisTK-3.7-93-ge9d8068/GenomeAnalysisTK.jar"
    Float gatk_version
    # cleaning will not work on Firecloud, so this has been set to always be false.
    Boolean clean = false
    Boolean extreme_deletion
    # TCIR Selection
    Boolean tcir
    # BQSR selection and relevant parameters
    Boolean bqsr
    Array[String] known_sites
    # HaplotypeCaller parameters
    Float interval_size
    String ? erc
    Int ? ploidy
    String ? extra_hc_params
    String ? bqsr_recal_report
    # GenotypeGVCF parameters
    String ? extra_gg_params
    Boolean ? all_sites
    # VQSR selection and relevant parameters
    Boolean vqsr
    Array[String] snp_resource
    Array[String] indel_resource
    Array[String] snp_annotation
    Array[String] indel_annotation
    Int ? snp_max_gaussians
    Int ? indel_max_gaussians
    Int ? mq_cap_snp
    Int ? mq_cap_indel
    Float ts_filter_snp
    Float ts_filter_indel
    String ? extra_vr_params
    # Hard filtration selection and relevant parameters
    Boolean combined_filtration
    Boolean variant_filtration
    String filter_expression
    String indel_filter_expression
    # SNPEff selection and relevant parameters
    Boolean use_snpeff
    String snpeff = "/cil/shed/apps/external/snpEff/snpEff-4.1g/snpEff.jar"
    String snpeff_db
    String ? snpeff_extra_params

    call VersionCheck{
        input:
        gatk = gatk,
        version = gatk_version
    }


#TODO remove
#    call **MakeOutputDir** {
#        input:
#        go = VersionCheck.out,
#        output_dir = output_dir
#    }

#TODO remove
    #call **SymlinkCromwellExecutionDir** {
    #    input:
    #    output_dir = ""
    #}

    # Check's if index files exist(using .dict file as marker). Will always output path to localized reference
    # With assumption that it either exists or will be created by IndexReference.
    call IndexReference {
        input:
        ref_fasta = ref_fasta,
        picard = picard,
    }

#TODO can use this output to scatter within a single bam
    call CreateIntervalsList {
        input:
        ref = IndexReference.out,
        interval_size = interval_size,
    }

# scatter by input sample.  
# each sample can have either 1 bam file or 2 fastq's
    scatter(sample in read_tsv(samples_file)) {
        # Within specified output_dir, create a subdir using sample_name specified in 1st column of tsv file.

        # one bam file provided
        if (length(sample) == 2) {
            call SamToFastq {
                input:
                picard = picard,
                in_bam = sample[1],
                sample_name = sample[0],
            }
        }
        # two fastq's provided.
        # TODO - this case is not implemented yet
#        if (length(sample) == 3) {


#            call copyFile as CopyFastq1{
#                input:
#                source = sample[1],
#                dest = ""

#            call copyFile as CopyFastq2{
#                input:
#                source = sample[2],
#                dest = ""


#            call **CopyFastq** as **CopyFastq**1 {
#                input:
#                fq = sample[1],
#                pair = 1,
#                sample_name = sample[0],
#            }

#            call **CopyFastq** as **CopyFastq**2 {
#                input:
#                fq = sample[2],
#                pair = 2,
#                sample_name = sample[0],
#            }
 
#        }
        # TODO later ought to rewire this use a 'select' statement to pick source of fastq's
        call AlignBAM {
            input:
            ref = IndexReference.out,
            dict = IndexReference.dict,
            amb = IndexReference.amb,
            ann = IndexReference.ann,
            bwt = IndexReference.bwt,
            fai = IndexReference.fai,
            pac = IndexReference.pac,
            sa = IndexReference.sa,
            sample_name = sample[0],
            fq1 = SamToFastq.fq1,
            fq2 = SamToFastq.fq2
        }
        call SortBAM {
            input:
            picard = picard,
            sample_name = sample[0],
            aligned_bam = AlignBAM.bam
        }
        # TODO delete cleaning step
        if (clean == true) {
            call removeFile as removeAlignedBam{input: file = AlignBAM.bam, ready = SortBAM.done}
        }
        call MarkDuplicates {
            input:
            picard = picard,
            sample_name = sample[0],
            sorted_bam = SortBAM.bam
        }
        # TODO delete cleaning step
        if (clean == true) {
            call removeFile as removeSortedSAM{ input: file = SortBAM.bam, ready = MarkDuplicates.done}
        }
        # TODO why reorderbam rather than sortbam??
        call ReorderBAM {
            input:
            picard = picard,
            sample_name = sample[0],
            marked_bam = MarkDuplicates.bam,
            ref = IndexReference.out,
            dict = IndexReference.dict
        }
        # TODO delete cleaning step
        if (clean == true) {
            call removeFile as removeMarkedDuplicates{input: file = MarkDuplicates.bam, ready = ReorderBAM.done}
        }
        # TODO why is a bam index being passed into index bam?
        # It is not consumed explicitly by anything downstream, and does not seem to be output anywhere
        #call IndexBAM {
        #    input:
        #    in_bam = ReorderBAM.bam,
        #    index = ReorderBAM.index
        #}

        # TODO doesn't indel-realigner require the bam to be re-indexed afterwards
        if (tcir == true) {
            call RealignerTargetCreator {
                input:
                gatk = gatk,
                sample_name = sample[0],
                ref = IndexReference.out,
                dict = IndexReference.dict,
                amb = IndexReference.amb,
                ann = IndexReference.ann,
                bwt = IndexReference.bwt,
                fai = IndexReference.fai,
                pac = IndexReference.pac,
                sa = IndexReference.sa,
                in_bam = ReorderBAM.bam,
                index = ReorderBAM.index
            }
            call IndelRealigner {
                input:
                gatk = gatk,
                ref = IndexReference.out,
                dict = IndexReference.dict,
                amb = IndexReference.amb,
                ann = IndexReference.ann,
                bwt = IndexReference.bwt,
                fai = IndexReference.fai,
                pac = IndexReference.pac,
                sa = IndexReference.sa,
                in_bam = ReorderBAM.bam,
                sample_name = sample[0],
                intervals = RealignerTargetCreator.intervals
            }
        # TODO delete cleaning step

            if (clean == true) {
                call removeFile as removeRTCI {input: file = RealignerTargetCreator.intervals, ready = IndelRealigner.done}
            }
            # copying out one of the final outputs...
            call copyFile as cpRealignedBam{
                input:
                source = IndelRealigner.bam,
                dest = ""
            }
        }
        if (bqsr == true) {
            # TODO why is bqsr run twice? it is expensive... Just for analyze covariates??
            File bqsr_bam = select_first([IndelRealigner.bam,ReorderBAM.bam])
            File bqsr_index = select_first([IndelRealigner.index, ReorderBAM.index])
            call BaseRecalibrator as BaseRecal_1 {
                input:
                gatk = gatk,
                ref = IndexReference.out,
                dict = IndexReference.dict,
                amb = IndexReference.amb,
                ann = IndexReference.ann,
                bwt = IndexReference.bwt,
                fai = IndexReference.fai,
                pac = IndexReference.pac,
                sa = IndexReference.sa,
                bam = bqsr_bam,
                index = bqsr_index,
                known_sites = known_sites,
                out_file = "recal_data.table"
            }

            call BaseRecalibrator as BaseRecal_2 {
                input:
                gatk = gatk,
                ref = IndexReference.out,
                dict = IndexReference.dict,
                amb = IndexReference.amb,
                ann = IndexReference.ann,
                bwt = IndexReference.bwt,
                fai = IndexReference.fai,
                pac = IndexReference.pac,
                sa = IndexReference.sa,
                bam = bqsr_bam,
                index = bqsr_index,
                known_sites = known_sites,
                bqsr = BaseRecal_1.out,
                out_file = "post_recal_data.table"
            }

            call AnalyzeCovariates {
                input:
                gatk = gatk,
                ref = IndexReference.out,
                dict = IndexReference.dict,
                amb = IndexReference.amb,
                ann = IndexReference.ann,
                bwt = IndexReference.bwt,
                fai = IndexReference.fai,
                pac = IndexReference.pac,
                sa = IndexReference.sa,
                sample_name = sample[0],
                before = BaseRecal_1.out,
                after = BaseRecal_2.out
            }
            # copy files to final output area 
            call copyFile as cpAnalyzeCovariatesOut{
                input:
                source = AnalyzeCovariates.out,
                dest = ""
            }
            call PrintReads {
                input:
                gatk = gatk,
                sample_name = sample[0],
                ref = IndexReference.out,
                dict = IndexReference.dict,
                amb = IndexReference.amb,
                ann = IndexReference.ann,
                bwt = IndexReference.bwt,
                fai = IndexReference.fai,
                pac = IndexReference.pac,
                sa = IndexReference.sa,
                in_bam = bqsr_bam,
                bqsr = BaseRecal_2.out
            }
            # TODO remove the cleaning
            if (clean == true) {
                call removeFile as removeBR1{input: file = BaseRecal_1.out, ready = AnalyzeCovariates.done}
                call removeFile as removeBR2{input: file = BaseRecal_2.out, ready = PrintReads.done}
            }
            # copy files to final output area 
            call copyFile as cpPrintReadsBam{
                input:
                source = PrintReads.bam,
                dest = ""
            }
        }
        File hc_bam = select_first([PrintReads.bam, IndelRealigner.bam, ReorderBAM.bam])
        call HaplotypeCaller {
            input:
            gatk = gatk,
            ref = IndexReference.out,
            dict = IndexReference.dict,
            amb = IndexReference.amb,
            ann = IndexReference.ann,
            bwt = IndexReference.bwt,
            fai = IndexReference.fai,
            pac = IndexReference.pac,
            sa = IndexReference.sa,
            sample_name = sample[0],
            in_bam = hc_bam,
            intervals = CreateIntervalsList.out,
            bqsr_file = bqsr_recal_report,
            ploidy = ploidy,
            erc = erc,
            extra_hc_params = extra_hc_params
        }
        # TODO remove cleaning steps
        if (clean == true) {
            if (hc_bam == ReorderBAM.bam) {
                call removeFile as removeReorderedBAM{input: file = ReorderBAM.bam, ready = HaplotypeCaller.done}
            }
            if (extreme_deletion == true) {
                call removeFile as removeFinalBam{input: file = hc_bam, ready = HaplotypeCaller.done}
            }
        }
        # copy files to final output area 
        call copyFile as cpHaplotypeCallerVcf{
            input:
            source = HaplotypeCaller.vcf,
            dest = ""
        }
    }
    # Scatter block ends

    call GenotypeGVCFs {
	    input:
        gatk = gatk,
        ref = IndexReference.out,
        dict = IndexReference.dict,
        amb = IndexReference.amb,
        ann = IndexReference.ann,
        bwt = IndexReference.bwt,
        fai = IndexReference.fai,
        pac = IndexReference.pac,
        sa = IndexReference.sa,
        extra_gg_params = extra_gg_params,
        all_sites = all_sites,
        intervals = CreateIntervalsList.out,
        variant_files = HaplotypeCaller.vcf # array of files
    }
    if (vqsr == true) {
        call VariantRecalibrator as SnpRecalibration {
            input:
            gatk = gatk,
            ref = IndexReference.out,
            dict = IndexReference.dict,
            amb = IndexReference.amb,
            ann = IndexReference.ann,
            bwt = IndexReference.bwt,
            fai = IndexReference.fai,
            pac = IndexReference.pac,
            sa = IndexReference.sa,
            intervals = CreateIntervalsList.out,
            task_input = GenotypeGVCFs.out,
            resource = snp_resource,
            annotation = snp_annotation,
            mode = "snp",
            max_gaussians = snp_max_gaussians,
            mq_cap = mq_cap_snp,
            extra_vr_params = extra_vr_params
        }

        call ApplyRecalibration as ApplySnpRecalibration {
            input:
            gatk = gatk,
            ref = IndexReference.out,
            dict = IndexReference.dict,
            amb = IndexReference.amb,
            ann = IndexReference.ann,
            bwt = IndexReference.bwt,
            fai = IndexReference.fai,
            pac = IndexReference.pac,
            sa = IndexReference.sa,
            vcf_in = GenotypeGVCFs.out,
            ts_filter = ts_filter_snp,
            recal_file = SnpRecalibration.recal,
            tranches = SnpRecalibration.tranches,
            mode = "snp",
            prefix = "snp"
        }
        call VariantRecalibrator as IndelRecalibration {
            input:
            gatk = gatk,
            ref = IndexReference.out,
            dict = IndexReference.dict,
            amb = IndexReference.amb,
            ann = IndexReference.ann,
            bwt = IndexReference.bwt,
            fai = IndexReference.fai,
            pac = IndexReference.pac,
            sa = IndexReference.sa,
            intervals = CreateIntervalsList.out,
            task_input = ApplySnpRecalibration.out,
            resource = indel_resource,
            mode = "indel",
            max_gaussians = indel_max_gaussians,
            mq_cap = mq_cap_indel,
            annotation = indel_annotation,
            extra_vr_params = extra_vr_params
        }
        call ApplyRecalibration as ApplyIndelRecalibration {
            input:
            gatk = gatk,
            ref = IndexReference.out,
            dict = IndexReference.dict,
            amb = IndexReference.amb,
            ann = IndexReference.ann,
            bwt = IndexReference.bwt,
            fai = IndexReference.fai,
            pac = IndexReference.pac,
            sa = IndexReference.sa,
            vcf_in = ApplySnpRecalibration.out,
            ts_filter = ts_filter_indel,
            recal_file = IndelRecalibration.recal,
            tranches = IndelRecalibration.tranches,
            mode = "indel",
            prefix = "snp_indel"
        }
        # TODO remove cleaning step
        if (clean == true) {
            call removeFile as RemoveGenoTypeGCVF{input: file =  GenotypeGVCFs.out, ready=ApplySnpRecalibration.done}
            call removeFile as removeIntervalsVQSR{input: file = CreateIntervalsList.out, ready = ApplyIndelRecalibration.done}
        }
    }
    # TODO remove cleaning step
    if (clean == true) {
        if (vqsr == false) {
        call removeFile as removeIntervals{input: file = CreateIntervalsList.out, ready = GenotypeGVCFs.done}
        }
    }
    File sv_vcf = select_first([ApplyIndelRecalibration.out, GenotypeGVCFs.out])
    if (combined_filtration == true) {
        call HardFiltration as CombinedFiltration {
            input:
            gatk = gatk,
            ref = IndexReference.out,
            dict = IndexReference.dict,
            amb = IndexReference.amb,
            ann = IndexReference.ann,
            bwt = IndexReference.bwt,
            fai = IndexReference.fai,
            pac = IndexReference.pac,
            sa = IndexReference.sa,
            vcf_in = sv_vcf,
            variant_type = "ALL",
            filter_expression = filter_expression
        }
    }
    if (variant_filtration == true) {
        call SelectVariants as SelectSnps{
            input:
            gatk = gatk,
            dict = IndexReference.dict,
            amb = IndexReference.amb,
            ann = IndexReference.ann,
            bwt = IndexReference.bwt,
            fai = IndexReference.fai,
            pac = IndexReference.pac,
            sa = IndexReference.sa,
            ref = IndexReference.out,
            vcf_in = sv_vcf,
            mode = "SNP"
        }
        call HardFiltration as FilterSnps{
            input:
            gatk = gatk,
            ref = IndexReference.out,
            dict = IndexReference.dict,
            amb = IndexReference.amb,
            ann = IndexReference.ann,
            bwt = IndexReference.bwt,
            fai = IndexReference.fai,
            pac = IndexReference.pac,
            sa = IndexReference.sa,
            vcf_in = SelectSnps.out,
            variant_type = "SNPs",
            filter_expression = filter_expression
        }
        call SelectVariants as SelectIndels{
            input:
            gatk = gatk,
            dict = IndexReference.dict,
            amb = IndexReference.amb,
            ann = IndexReference.ann,
            bwt = IndexReference.bwt,
            fai = IndexReference.fai,
            pac = IndexReference.pac,
            sa = IndexReference.sa,
            ref = IndexReference.out,
            vcf_in = sv_vcf,
            mode = "INDEL"
        }
        call HardFiltration as FilterIndels{
            input:
            gatk = gatk,
            ref = IndexReference.out,
            dict = IndexReference.dict,
            amb = IndexReference.amb,
            ann = IndexReference.ann,
            bwt = IndexReference.bwt,
            fai = IndexReference.fai,
            pac = IndexReference.pac,
            sa = IndexReference.sa,
            vcf_in = SelectIndels.out,
            variant_type = "INDELS",
            filter_expression = indel_filter_expression
        }
        call CombineVariants {
            input:
            gatk = gatk,
            dict = IndexReference.dict,
            amb = IndexReference.amb,
            ann = IndexReference.ann,
            bwt = IndexReference.bwt,
            fai = IndexReference.fai,
            pac = IndexReference.pac,
            sa = IndexReference.sa,
            ref = IndexReference.out,
            vcf1 = FilterSnps.out,
            vcf2 = FilterIndels.out
        }
        # TODO remove cleaning
        if (clean == true) {
            call removeFile as RemoveSV_VCF{input: file = sv_vcf, ready=SelectSnps.done}
            call removeFile as RemoveSelectedSNPs{input: file = SelectSnps.out, ready= FilterSnps.done}
            call removeFile as RemoveSelectedIndels{input: file = SelectIndels.out, ready = FilterIndels.done}
        }
    }
    #TODO either there should be enforcement that the caller does not run both variant_filtration and combined_filtration, or it ought to handle the case where both are selected.
    File snpeff_vcf = select_first([CombinedFiltration.out, CombineVariants.out])
    if (use_snpeff == true) {
        call SnpEff {
            input:
            vcf_in = snpeff_vcf,
            snpeff_db = snpeff_db,
            snpeff = snpeff,
            snpeff_extra_params = snpeff_extra_params
        }
    }
    # Final handoff steps
    Array[File] handoff_files = select_all([SnpRecalibration.tranches, SnpRecalibration.recal, SnpRecalibration.rscript, IndelRecalibration.tranches, IndelRecalibration.recal, IndelRecalibration.rscript, ApplyIndelRecalibration.out, FilterSnps.out, FilterIndels.out, snpeff_vcf, SnpEff.out])
    scatter(file in handoff_files) {
        call copyFile as copyFinal {
            input:
            source = file,
            dest = ""
        }
    }
}
