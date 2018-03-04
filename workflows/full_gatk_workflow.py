from runcromwellremote import run_wdl


def check_run_wdl(wdl_path, inputs_dict)
    passed, outputs = run_wdl(wdl_path, inputs_dict)
    if not passed:
        sys.exit(1)
    return outputs


def process_sample(inputs)




    tcir_wdl_path = 'btl_gatk_tcir/taskdef.btl_gatk_tcir.wdl'
    tcir_inputs = {
        'gatk_tcir.gatk_tcir_task.reference_tgz':inputs['reference_tgz'],
        'gatk_tcir.gatk_tcir_task.in_bam':inputs['in_bam'],
        'gatk_tcir.gatk_tcir_task.in_bam_index':inputs['in_bam_index'],
        'gatk_tcir.gatk_tcir_task.output_disk_gb':'10',
        'gatk_tcir.gatk_tcir_task.sample_name':inputs['sample_name'],
        'gatk_tcir.gatk_tcir_task.debug_dump_flag':'onfail',

        }
    tcir_outputs = check_run_wdl(tcir_wdl_path, tcir_inputs)




    bqsr_wdl_path = 'btl_gatk_bqsr/taskdef.btl_gatk_bqsr.wdl'
    bqsr_inputs = {
        'gatk_bqsr.indelrealigner_bam':tcir_outputs['out_bam'],
        'gatk_bqsr.indelrealigner_bam_index':tcir_outputs['out_bam_index'],
        'gatk_bqsr.gatk_bqsr_task.debug_dump_flag':'onfail',
        'gatk_bqsr.gatk_bqsr_task.output_disk_gb':'10',
        'gatk_bqsr.gatk_bqsr_task.reference_tgz':inputs['reference_tgz'],
        'gatk_bqsr.gatk_bqsr_task.known_sites_vcfs':["gs://broad-cil-devel-bucket/input_data/7g8_gb4.combined.final.vcf.gz","gs://broad-cil-devel-bucket/input_data/hb3_dd2.combined.final.vcf.gz","gs://broad-cil-devel-bucket/input_data/3d7_hb3.combined.final.vcf.gz"],
        'gatk_bqsr.gatk_bqsr_task.known_sites_vcf_tbis':["gs://broad-cil-devel-bucket/input_data/7g8_gb4.combined.final.vcf.gz.tbi","gs://broad-cil-devel-bucket/input_data/hb3_dd2.combined.final.vcf.gz.tbi","gs://broad-cil-devel-bucket/input_data/3d7_hb3.combined.final.vcf.gz.tbi"],
        'gatk_bqsr.gatk_bqsr_task.sample_name':inputs['sample_name'],
    }
    bqsr_outputs = check_run_wdl(bqsr_wdl_path, bqsr_inputs)



    haplotypecaller_wdl = 'btl_gatk_haplotypecaller/taskdef.btl_gatk_haplotypecaller.wdl'
    haplotypecaller_inputs = {
        'gatk_haplotypecaller.gatk_haplotypecaller_task.reference_tgz':inputs['reference_tgz'],
        'gatk_haplotypecaller.indelrealigner_bam':bqsr_outputs['out_bam']
        'gatk_haplotypecaller.indelrealigner_bam_index':bqsr_outputs['out_bam_index']
        'gatk_haplotypecaller.gatk_haplotypecaller_task.output_disk_gb':'10',
        'gatk_haplotypecaller.gatk_haplotypecaller_task.sample_name':inputs['sample_name'],
        'gatk_haplotypecaller.gatk_haplotypecaller_task.debug_dump_flag':'onfail',

        'gatk_haplotypecaller.gatk_haplotypecaller_task.known_sites_vcf_tbis':["gs://broad-cil-devel-bucket/input_data/7g8_gb4.combined.final.vcf.gz.tbi","gs://broad-cil-devel-bucket/input_data/hb3_dd2.combined.final.vcf.gz.tbi","gs://broad-cil-devel-bucket/input_data/3d7_hb3.combined.final.vcf.gz.tbi"],
        'gatk_haplotypecaller.gatk_haplotypecaller_task.known_sites_vcfs':["gs://broad-cil-devel-bucket/input_data/7g8_gb4.combined.final.vcf.gz","gs://broad-cil-devel-bucket/input_data/hb3_dd2.combined.final.vcf.gz","gs://broad-cil-devel-bucket/input_data/3d7_hb3.combined.final.vcf.gz"],

        }
    haplotypecaller_outputs = check_run_wdl(haplotypecaller_wdl_path, haplotypecaller_inputs)




