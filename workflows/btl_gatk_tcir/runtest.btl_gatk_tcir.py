from runcromwellremote import run_wdl

wdl = 'taskdef.btl_gatk_tcir.wdl'

inputs = {
    'gatk_tcir.gatk_tcir_task.reference_tgz':'gs://broad-cil-devel-bucket/input_data/minion_illumina_hybrid_clean_MT.tgz',
    'gatk_tcir.gatk_tcir_task.in_bam':'gs://broad-cil-devel-bucket/input_data/Candida_Auris.tcir.bam',
    'gatk_tcir.gatk_tcir_task.in_bam_index':'gs://broad-cil-devel-bucket/input_data/Candida_Auris.tcir.bam.bai',
    'gatk_tcir.gatk_tcir_task.output_disk_gb':'10',
    'gatk_tcir.gatk_tcir_task.sample_name':'Candida_Auris',
    'gatk_tcir.gatk_tcir_task.debug_dump_flag':'onfail',

    }

(status,outputs) = run_wdl(wdl,inputs)

print(outputs)

