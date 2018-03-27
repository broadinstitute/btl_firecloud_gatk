from runcromwellremote import run_wdl

wdl = 'taskdef.btl_gatk_process_cohort.wdl'


inputs = {
    'gatk_process_cohort.reference_tgz':'gs://broad-cil-devel-bucket/input_data/minion_illumina_hybrid_clean_MT.tgz',
    'gatk_process_cohort.gvcf_list':["gs://broad-cil-devel-bucket/gatk_haplotypecaller/dbf79433-58e5-4931-b77c-3fcf801f3db2/call-gatk_haplotypecaller_task/Candida_Auris.gvcf"],
    'gatk_process_cohort.cohort_name':'Candida_Auris_cohort',
    }

(status,outputs) = run_wdl(wdl,inputs, workflow_dependencies="wdl_bundle.zip")

print(outputs)
