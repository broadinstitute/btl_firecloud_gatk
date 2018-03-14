from runcromwellremote import run_wdl

wdl = 'taskdef.btl_gatk_filtration.wdl'

inputs = {
    'gatk_filtration.gatk_filtration_task.reference_tgz':'gs://broad-cil-devel-bucket/input_data/minion_illumina_hybrid_clean_MT.tgz',
    'gatk_filtration.gatk_filtration_task.output_disk_gb':'10',
    'gatk_filtration.gatk_filtration_task.debug_dump_flag':'onfail',
    'gatk_filtration.gatk_filtration_task.cohort_name':"test_cohort",
    'gatk_filtration.genotype_caller_vcf':"gs://broad-cil-devel-bucket/gatk_joint_genotype/5929ccf3-d281-45ad-bf63-20d7b558184a/call-gatk_joint_genotype_task/test_cohort.vcf",
    'gatk_filtration.gatk_filtration_task.filter_expression':"'VQSLOD <= 0.0'",
    'gatk_filtration.filtration_type':'variant',

    
    }

(status,outputs) = run_wdl(wdl,inputs)

print(outputs)

               
