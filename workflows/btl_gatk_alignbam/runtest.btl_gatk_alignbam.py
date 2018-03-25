from runcromwellremote import run_wdl

wdl = 'taskdef.latecompress.btl_gatk_alignbam.wdl'

inputs = {
    'gatk_alignbam.gatk_alignbam_task.reference_tgz':'gs://broad-cil-devel-bucket/gatk_indexref/0e409344-5f52-43dc-b811-4b00f517226c/call-gatk_indexref_task/AgPEST_v4.tgz',
    'gatk_alignbam.gatk_alignbam_task.in_bam':"gs://broad-cil-devel-bucket/input_data/171226_Anopheles_Longranger_topoff__Mali1_wgs__phased_possorted_bam.bam",
    'gatk_alignbam.gatk_alignbam_task.output_disk_gb':'10000',
    'gatk_alignbam.gatk_alignbam_task.sample_name':'Mali1_wgs',
    'gatk_alignbam.gatk_alignbam_task.debug_dump_flag':'onfail',

    }

(status,outputs) = run_wdl(wdl,inputs)

print(outputs)