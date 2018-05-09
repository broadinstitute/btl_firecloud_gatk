import "taskdef.btl_gatk_alignbam.wdl" as btl_gatk_alignbam
import "taskdef.btl_gatk_tcir.wdl" as btl_gatk_tcir
import "taskdef.btl_gatk_bqsr.wdl" as btl_gatk_bqsr
import "taskdef.btl_gatk_haplotypecaller.wdl" as btl_gatk_haplotypecaller
import "taskdef.btl_gatk_indexref.wdl" as btl_gatk_indexref


workflow gatk_process_samples {
    String? onprem_download_path
    File samples_tsv_fofn

    Boolean prealigned = false
    Boolean use_tcir
    Boolean use_bqsr

    File reference_fasta
    Array[File] known_sites_vcfs
    Array[File] known_sites_vcf_tbis
    String ref_base_name = basename(reference_fasta)
    # This regex should cover all three species of fasta files seen in the wild: .fasta, .fa, and .fna. Tested at regex101.com
    String ref_name = sub(ref_base_name, ".f[nasta]", "")


    call btl_gatk_indexref.gatk_indexref_task as gatk_indexref_task {
        input:
            ref_fasta = reference_fasta,
            ref_name = ref_name,
            output_disk_gb = "10",
            debug_dump_flag = "onfail"
    }


    scatter (sample_entry in read_tsv(samples_tsv_fofn)) {
        File bam_entry = sample_entry[1]
        String sample_name = sample_entry[0]

        if (prealigned) {
            File prealigned_bam_idx = sample_entry[2]
            File prealigned_bam = bam_entry
        }
        
        if (!prealigned) {
            call btl_gatk_alignbam.gatk_alignbam_task as gatk_alignbam_task{
                input:
                    in_bam = bam_entry,
                    sample_name = sample_name,
                    reference_tgz = gatk_indexref_task.reference_tgz,
                    output_disk_gb = "10",
                    debug_dump_flag = "onfail"
            }
        }

        File in_bam = select_first([prealigned_bam, gatk_alignbam_task.out_bam])
        File in_bam_index = select_first([prealigned_bam_idx, gatk_alignbam_task.out_bam_index])

        if (use_tcir) {
            call btl_gatk_tcir.gatk_tcir_task as gatk_tcir_task{
                input:
                    in_bam = in_bam,
                    in_bam_index = in_bam_index,
                    sample_name = sample_name,
                    reference_tgz = gatk_indexref_task.reference_tgz,
                    output_disk_gb = "10",
                    debug_dump_flag = "onfail"
            }
        }

        # note the conditional possibility (that tcir runs) has to be first option in select_first or it will never be assigned to bqsr_in_bam.
        File in_bam_two = select_first([gatk_tcir_task.out_bam, in_bam])
        File in_bam_two_index = select_first([gatk_tcir_task.out_bam_index, in_bam_index])

        if (use_bqsr) {
            call btl_gatk_bqsr.gatk_bqsr_task as gatk_bqsr_task{
                input:
                    in_bam = in_bam_two,
                    in_bam_index = in_bam_two_index,
                    known_sites_vcfs = known_sites_vcfs,
                    known_sites_vcf_tbis = known_sites_vcf_tbis,
                    sample_name = sample_name,
                    reference_tgz = gatk_indexref_task.reference_tgz,
                    output_disk_gb = "10",
                    debug_dump_flag = "onfail"
            }
        }

        File hc_in_bam = select_first([gatk_bqsr_task.out_bam, in_bam_two])
        File hc_in_bam_index = select_first([gatk_bqsr_task.out_bam_index, in_bam_two_index])
        call btl_gatk_haplotypecaller.gatk_haplotypecaller_task as gatk_haplotypecaller_task {
	        input:
                in_bam = hc_in_bam,
                in_bam_index = hc_in_bam_index,
                bqsr_table = gatk_bqsr_task.out_bqsr_table,
                sample_name = sample_name,
                reference_tgz = gatk_indexref_task.reference_tgz,
                output_disk_gb = "10",
                debug_dump_flag = "onfail"
	    }


        output {
            Array[File] out_gvcf = gatk_haplotypecaller_task.out_gvcf
        }
    }
}