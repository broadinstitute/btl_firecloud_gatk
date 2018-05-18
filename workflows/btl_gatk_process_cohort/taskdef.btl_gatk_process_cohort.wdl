import "taskdef.btl_gatk_joint_genotype.wdl" as btl_gatk_joint_genotype
import "taskdef.btl_gatk_vqsr.wdl" as btl_gatk_vqsr
import "taskdef.btl_gatk_variant_filtration.wdl" as btl_gatk_variant_filtration
import "taskdef.btl_gatk_filter_genotypes.wdl" as btl_gatk_filter_genotypes
import "taskdef.btl_gatk_snpeff.wdl" as btl_gatk_snpeff

workflow gatk_process_cohort {
    String? onprem_download_path
    Map[String, String]? handoff_files

    File gvcf_fofn
    # I've made gvcf_list 'optional' just so widdler validate doesn't complain when it isn't in json file, since
    # this is extrapolated from the fofn.
    Array[File] gvcf_list = read_lines(gvcf_fofn)
    String cohort_name
    File reference_tgz
    File snpeff_db_tgz
    String snpeff_db_name

    Array[String] snp_annotation
    Array[String] indel_annotation
    Array[File] known_sites_vcfs
    Array[File] known_sites_vcf_tbis
    Array[String] snp_resource_params
    Array[String] indel_resource_params

    Int snp_max_gaussians
    Int indel_max_gaussians
    Int mq_cap_snp
    Int mq_cap_indel
    Float ts_filter_snp
    Float ts_filter_indel
    String ? extra_vr_params
    String output_disk_gb = "10"

    Boolean run_gatk_vqsr = true
    Boolean run_gatk_filter_genotypes

    call btl_gatk_joint_genotype.gatk_joint_genotype_task as gatk_joint_genotype_task{
        input:
            HaplotypeCaller_gvcfs = gvcf_list,
            cohort_name = cohort_name,
            reference_tgz = reference_tgz,
            output_disk_gb = "10",
            debug_dump_flag = "onfail"
    }

    if (run_gatk_vqsr == true) {
        call btl_gatk_vqsr.gatk_vqsr_task as gatk_vqsr_task {
            input:
                cohort_name = cohort_name,
                genotype_caller_vcf = gatk_joint_genotype_task.vcf_out,
                reference_tgz = reference_tgz,
                snp_annotation = snp_annotation,
                indel_annotation = indel_annotation,
                snp_max_gaussians = snp_max_gaussians,
                indel_max_gaussians = indel_max_gaussians,
                snp_resource_params= snp_resource_params,
                indel_resource_params = indel_resource_params,
                mq_cap_snp = mq_cap_snp,
                mq_cap_indel = mq_cap_indel,
                ts_filter_snp = ts_filter_snp,
                ts_filter_indel = ts_filter_indel,
                known_sites_vcf_tbis = known_sites_vcf_tbis,
                known_sites_vcfs =known_sites_vcfs,
                extra_vr_params = extra_vr_params,
                output_disk_gb = output_disk_gb,
                debug_dump_flag = "onfail"
        }
    }
    Array[File] filtration_vcf_in = select_all([gatk_vqsr_task.vcf_out, gatk_joint_genotype_task.vcf_out])
    call btl_gatk_variant_filtration.gatk_variant_filtration_task as gatk_variant_filtration_task{
        input:
            sv_vcf = filtration_vcf_in[0],
            cohort_name = cohort_name,
            reference_tgz = reference_tgz,
            output_disk_gb = output_disk_gb,
            debug_dump_flag = "onfail",
            snp_filter_expression = "VQSLOD <= 0.0",
            indel_filter_expression = "VQSLOD <= 0.0"

    }

    if (run_gatk_filter_genotypes == true) {
        call btl_gatk_filter_genotypes.gatk_filter_genotypes_task as gatk_filter_genotypes_task {
            input:
                vcf_in = gatk_variant_filtration_task.vcf_out,
                cohort_name = cohort_name,
                output_disk_gb = output_disk_gb,
                debug_dump_flag = "onfail"
        }
    }
    Array[File] snpeff_vcf_in = select_all([gatk_filter_genotypes_task.vcf_out, gatk_variant_filtration_task.vcf_out])
    call btl_gatk_snpeff.gatk_snpeff_task as gatk_snpeff_task {
        input:
            cohort_name = cohort_name,
            vcf_in = snpeff_vcf_in[0],
            snpeff_db_tgz = snpeff_db_tgz,
            snpeff_db_name = snpeff_db_name,
            output_disk_gb = output_disk_gb,
            debug_dump_flag = "onfail"
    }

    output {
        File out_vcf = gatk_variant_filtration_task.vcf_out
        }
}



