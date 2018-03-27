from runcromwellremote import run_wdl

wdl = 'taskdef.btl_gatk_process_sample.wdl'


inputs = {
    'gatk_process_sample.reference_tgz':'gs://broad-cil-devel-bucket/input_data/minion_illumina_hybrid_clean_MT.tgz',
    'gatk_process_sample.in_bam':"gs://broad-cil-devel-bucket/input_data/Candida_Auris.bam",
    'gatk_process_sample.sample_name':'Candida_Auris',
    'gatk_process_sample.known_sites_vcf_tbis':["gs://broad-cil-devel-bucket/input_data/7g8_gb4.combined.final.vcf.gz.tbi","gs://broad-cil-devel-bucket/input_data/hb3_dd2.combined.final.vcf.gz.tbi","gs://broad-cil-devel-bucket/input_data/3d7_hb3.combined.final.vcf.gz.tbi"],
    'gatk_process_sample.known_sites_vcfs':["gs://broad-cil-devel-bucket/input_data/7g8_gb4.combined.final.vcf.gz","gs://broad-cil-devel-bucket/input_data/hb3_dd2.combined.final.vcf.gz","gs://broad-cil-devel-bucket/input_data/3d7_hb3.combined.final.vcf.gz"],
    }

(status,outputs) = run_wdl(wdl,inputs, workflow_dependencies="wdl_bundle.zip")

print(outputs)
