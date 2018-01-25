workflow gatk_vqsr{

    String gatk = "/humgen/gsa-hpprojects/GATK/bin/GenomeAnalysisTK-3.7-93-ge9d8068/GenomeAnalysisTK.jar"

    File index_reference_out
    File index_reference_dict
    File index_reference_amb
    File index_reference_ann
    File index_reference_bwt
    File index_reference_fai
    File index_reference_pac
    File index_reference_sa
    String sample_name
    #TODO how should intervals list be wired?
    File CreateIntervalsList_out 
    Array[File] HaplotypeCaller_vcfs
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



    call VariantRecalibrator as SnpRecalibration {
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
        intervals = CreateIntervalsList_out,
        task_input = HaplotypeCaller_vcfs,
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
        ref = index_reference_out,
        dict = index_reference_dict,
        amb = index_reference_amb,
        ann = index_reference_ann,
        bwt = index_reference_bwt,
        fai = index_reference_fai,
        pac = index_reference_pac,
        sa = index_reference_sa,
        vcf_in = HaplotypeCaller_vcfs,
        ts_filter = ts_filter_snp,
        recal_file = SnpRecalibration.recal,
        tranches = SnpRecalibration.tranches,
        mode = "snp",
        prefix = "snp"
    }
    call VariantRecalibrator as IndelRecalibration {
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
        intervals = CreateIntervalsList_out,
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
        ref = index_reference_out,
        dict = index_reference_dict,
        amb = index_reference_amb,
        ann = index_reference_ann,
        bwt = index_reference_bwt,
        fai = index_reference_fai,
        pac = index_reference_pac,
        sa = index_reference_sa,
        vcf_in = ApplySnpRecalibration.out,
        ts_filter = ts_filter_indel,
        recal_file = IndelRecalibration.recal,
        tranches = IndelRecalibration.tranches,
        mode = "indel",
        prefix = "snp_indel"
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
