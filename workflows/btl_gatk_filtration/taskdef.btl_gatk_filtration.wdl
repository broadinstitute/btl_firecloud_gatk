workflow gatk_filtration {

    File index_reference_out
    File index_reference_dict
    File index_reference_amb
    File index_reference_ann
    File index_reference_bwt
    File index_reference_fai
    File index_reference_pac
    File index_reference_sa
    String sample_name
    String gatk = "/humgen/gsa-hpprojects/GATK/bin/GenomeAnalysisTK-3.7-93-ge9d8068/GenomeAnalysisTK.jar"

    File ? vqsr_vcf
    File ? genotype_vcf

    File sv_vcf = select_first([vqsr_vcf, genotype_vcf])
    String filtration_type #combined or variant
    String filter_expression
    String indel_filter_expression

       # File snpeff_vcf = select_first([CombinedFiltration.out, CombineVariants.out])


    if (filtration_type == "combined") {
        call HardFiltration as CombinedFiltration {
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
            vcf_in = sv_vcf,
            variant_type = "ALL",
            filter_expression = filter_expression
        }
    }
    if (filtration_type == "variant") {
        call SelectVariants as SelectSnps{
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
            vcf_in = sv_vcf,
            mode = "SNP"
        }
        call HardFiltration as FilterSnps{
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
            variant_type = "SNPs",
            filter_expression = filter_expression
        }
        call SelectVariants as SelectIndels{
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
            vcf_in = sv_vcf,
            mode = "INDEL"
        }
        call HardFiltration as FilterIndels{
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
            vcf_in = SelectIndels.out,
            variant_type = "INDELS",
            filter_expression = indel_filter_expression
        }
        call CombineVariants {
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
            vcf1 = FilterSnps.out,
            vcf2 = FilterIndels.out

        }
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