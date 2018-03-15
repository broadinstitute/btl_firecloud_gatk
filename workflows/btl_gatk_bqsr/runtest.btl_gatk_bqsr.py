from runcromwellremote import run_wdl

wdl = 'taskdef.btl_gatk_bqsr.wdl'


inputs = {
    'gatk_bqsr.gatk_bqsr_task.reference_tgz':'gs://broad-cil-devel-bucket/input_data/minion_illumina_hybrid_clean_MT.tgz',
    'gatk_bqsr.gatk_bqsr_task.in_bam':"gs://broad-cil-devel-bucket/gatk_tcir/235faed4-cc6b-4d79-a213-843af2f3a58f/call-gatk_tcir_task/Candida_Auris.tcir.bam",
    'gatk_bqsr.gatk_bqsr_task.in_bam_index':'gs://broad-cil-devel-bucket/gatk_tcir/235faed4-cc6b-4d79-a213-843af2f3a58f/call-gatk_tcir_task/Candida_Auris.tcir.bam.bai',
    'gatk_bqsr.gatk_bqsr_task.output_disk_gb':'10',
    'gatk_bqsr.gatk_bqsr_task.sample_name':'Candida_Auris',
    'gatk_bqsr.gatk_bqsr_task.debug_dump_flag':'onfail',
    'gatk_bqsr.gatk_bqsr_task.known_sites_vcfs':["gs://broad-cil-devel-bucket/input_data/7g8_gb4.combined.final.vcf.gz","gs://broad-cil-devel-bucket/input_data/hb3_dd2.combined.final.vcf.gz","gs://broad-cil-devel-bucket/input_data/3d7_hb3.combined.final.vcf.gz"],
    'gatk_bqsr.gatk_bqsr_task.known_sites_vcf_tbis':["gs://broad-cil-devel-bucket/input_data/7g8_gb4.combined.final.vcf.gz.tbi","gs://broad-cil-devel-bucket/input_data/hb3_dd2.combined.final.vcf.gz.tbi","gs://broad-cil-devel-bucket/input_data/3d7_hb3.combined.final.vcf.gz.tbi"],

    }

(status,outputs) = run_wdl(wdl,inputs)

print(outputs)