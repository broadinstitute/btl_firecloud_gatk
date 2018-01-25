
workflow gatk_snpeff {
    File ? vqsr_vcf
    File ? genotype_vcf
    File ? filtration_vcf

    File snpeff_vcf = select_first([filtration_vcf, vqsr_vcf, genotype_vcf])

    String snpeff = "/cil/shed/apps/external/snpEff/snpEff-4.1g/snpEff.jar"
    String snpeff_db
    String ? snpeff_extra_params

    call SnpEff {
        input:
        vcf_in = snpeff_vcf,
        snpeff_db = snpeff_db,
        snpeff = snpeff,
        snpeff_extra_params = snpeff_extra_params
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
