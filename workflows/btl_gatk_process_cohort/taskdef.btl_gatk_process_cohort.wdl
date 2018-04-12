import "taskdef.btl_gatk_joint_genotype.wdl" as btl_gatk_joint_genotype
import "taskdef.btl_gatk_vqsr.wdl" as btl_gatk_vqsr
import "taskdef.btl_gatk_variant_filtration.wdl" as btl_gatk_variant_filtration
import "taskdef.btl_gatk_filter_genotypes.wdl" as btl_gatk_filter_genotypes
import "taskdef.btl_gatk_snpeff.wdl" as btl_gatk_snpeff

workflow gatk_process_cohort {
    String? onprem_download_path
    Map[String, String]? handoff_files

    File gvcf_fofn
    # I've made gvcf_list 'optional' just so widdler validate doesn't complain when it isn't in json file. 
    Array[File] ? gvcf_list = read_lines(gvcf_fofn)
    String cohort_name
    File reference_tgz
    File snpeff_db_tgz
    String snpeff_db_name
    Boolean run_gatk_filter_genotypes

    call btl_gatk_joint_genotype.gatk_joint_genotype_task as gatk_joint_genotype_task{
        input:
            HaplotypeCaller_gvcfs = gvcf_list,
            cohort_name = cohort_name,
            reference_tgz = reference_tgz,
            output_disk_gb = "10",
            debug_dump_flag = "onfail"
    }


    call btl_gatk_variant_filtration.gatk_variant_filtration_task as gatk_variant_filtration_task{
        input:
            sv_vcf = gatk_joint_genotype_task.vcf_out, 
            cohort_name = cohort_name,
            reference_tgz = reference_tgz,
            output_disk_gb = "10",
            debug_dump_flag = "onfail",
            snp_filter_expression = "VQSLOD <= 0.0",
            indel_filter_expression = "VQSLOD <= 0.0"

    }
    if (run_gatk_filter_genotypes == true) {
        call btl_gatk_filter_genotypes.gatk_filter_genotypes_task as gatk_filter_genotypes_task {
            input:
                vcf_in = gatk_variant_filtration_task.vcf_out,
                cohort_name = cohort_name,
                output_disk_gb = "10",
                debug_dump_flag = "onfail"
        }
    }
    File snpeff_vcf_in = select_first([gatk_filter_genotypes_task.vcf_out, gatk_variant_filtration_task.vcf_out])
    call btl_gatk_snpeff.gatk_snpeff_task as gatk_snpeff_task {
        input:
            vcf_in = snpeff_vcf_in,
            snpeff_db_tgz = snpeff_db_tgz,
            snpeff_db_name = snpeff_db_name
    }

    output {
        File out_vcf = gatk_variant_filtration_task.vcf_out
        }
}



