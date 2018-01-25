workflow gatk_aggregate_calls {
    File index_reference_out
    File index_reference_dict
    File index_reference_amb
    File index_reference_ann
    File index_reference_bwt
    File index_reference_fai
    File index_reference_pac
    File index_reference_sa
    String ? extra_gg_params
    #TODO how should intervals be wired?
    File CreateIntervalsList_out 
    Array[File] HaplotypeCaller_vcfs
    Boolean ? all_sites
    String gatk = "/humgen/gsa-hpprojects/GATK/bin/GenomeAnalysisTK-3.7-93-ge9d8068/GenomeAnalysisTK.jar"

    call GenotypeGVCFs {
	    input:

        ref = index_reference_out,
        dict = index_reference_dict,
        amb = index_reference_amb,
        ann = index_reference_ann,
        bwt = index_reference_bwt,
        fai = index_reference_fai,
        pac = index_reference_pac,
        sa = index_reference_sa,
        extra_gg_params = extra_gg_params,
        intervals = CreateIntervalsList_out,
        HaplotypeCaller_vcfs = HaplotypeCaller_vcfs, # array of files
        all_sites = all_sites,
        gatk = gatk,

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
	Array[File] HaplotypeCaller_vcfs
    Boolean ? all_sites
    String gvcf_out = "genotypeGVCFs.vcf"
	command {
        set -euo pipefail
        /opt/src/algutil/monitor_start.py
		java -Xmx8G -jar ${gatk} \
			-T GenotypeGVCFs \
			-R ${ref} \
			${sep=" --intervals " "--intervals " + intervals} \
			-o ${gvcf_out} \
			-V ${sep=" -V " HaplotypeCaller_vcfs} \
			${true="-allSites" false="" all_sites} \
			${default="\n" extra_gg_params}
        /opt/src/algutil/monitor_stop.py

	}
	output {
		#To track additional outputs from your task, please manually add them below
		File out = gvcf_out
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