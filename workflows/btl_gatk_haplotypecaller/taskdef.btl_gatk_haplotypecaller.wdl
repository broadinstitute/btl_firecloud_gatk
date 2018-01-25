workflow gatk_haplotypecaller {


    File ? bqsr_bam
    File ? indelrealigner_bam
    File ? uncleaned_bam
    File ? bqsr_bam_index
    File ? indelrealigner_bam_index
    File ? uncleaned_bam_index
    File index_reference_out
    File index_reference_dict
    File index_reference_amb
    File index_reference_ann
    File index_reference_bwt
    File index_reference_fai
    File index_reference_pac
    File index_reference_sa
    String sample_name
    File ? bqsr_recal_report
    String ploidy
    String erc
    String extra_hc_params
    String gatk = "/humgen/gsa-hpprojects/GATK/bin/GenomeAnalysisTK-3.7-93-ge9d8068/GenomeAnalysisTK.jar"



#TODO can use this output to scatter within a single bam
    call CreateIntervalsList {
        input:
        ref = index_reference_out,
        dict = index_reference_dict,
        amb = index_reference_amb,
        ann = index_reference_ann,
        bwt = index_reference_bwt,
        fai = index_reference_fai,
        pac = index_reference_pac,
        sa = index_reference_sa,        interval_size = interval_size,
    }

    File hc_bam = select_first([bqsr_bam, indelrealigner_bam, uncleaned_bam])
    File hc_bam_index = select_first([bqsr_bam_index, indelrealigner_bam_index, uncleaned_bam_index])
    call HaplotypeCaller {
        input:
        gatk = gatk,
        ref = index_reference_out,
        dict = index_reference_dict,
        amb = index_reference_amb,
        ann = index_reference_ann,
        bwt = index_reference_bwt,
        fai = index_reference_fai,
        pac = index_reference_pac,
        sa = index_reference_sa,
        sample_name = sample_name,
        in_bam = hc_bam,
        in_bam_index = hc_bam_index,
        intervals = CreateIntervalsList.out,
        bqsr_file = bqsr_recal_report,
        ploidy = ploidy,
        erc = erc,
        extra_hc_params = extra_hc_params
    }


}


task CreateIntervalsList {
	File ref
    File dict
    File amb
    File ann
    File bwt
    File fai
    File pac
    File sa
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
    File in_bam_index
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
	    #samtools index ${in_bam}
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
		File gvcf = out
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