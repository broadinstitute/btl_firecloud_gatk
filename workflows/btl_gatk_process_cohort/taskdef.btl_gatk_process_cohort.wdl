import "taskdef.btl_gatk_joint_genotype.wdl" as btl_gatk_joint_genotype
import "taskdef.btl_gatk_vqsr.wdl" as btl_gatk_vqsr
import "taskdef.btl_gatk_variant_filtration.wdl" as btl_gatk_variant_filtration
import "taskdef.btl_gatk_snpeff.wdl" as btl_gatk_snpeff


workflow gatk_process_cohort {
    String? onprem_download_path
    Map[String, String]? handoff_files


    Array[File] gvcf_list
    String cohort_name
    File reference_tgz
    

    call btl_gatk_joint_genotype.gatk_joint_genotype_task as gatk_joint_genotype_task{
        input:
            HaplotypeCaller_gvcfs = gvcf_list,
            cohort_name = cohort_name,
            reference_tgz = reference_tgz,
            output_disk_gb = "10",
            debug_dump_flag = "onfail",
    }


    call btl_gatk_variant_filtration.gatk_variant_filtration_task as gatk_variant_filtration_task{
        input:
            sv_vcf = gatk_joint_genotype_task.vcf_out, 
            cohort_name = cohort_name,
            reference_tgz = reference_tgz,
            output_disk_gb = "10",
            debug_dump_flag = "onfail",
            snp_filter_expression = "VQSLOD <= 0.0",
            indel_filter_expression = "VQSLOD <= 0.0",

    }



    output {
        File out_vcf = gatk_variant_filtration_task.vcf_out
        }

}



