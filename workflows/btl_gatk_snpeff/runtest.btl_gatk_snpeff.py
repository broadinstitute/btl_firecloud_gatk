from runcromwellremote import run_wdl

wdl = 'taskdef.btl_gatk_snpeff.wdl'

inputs = {
    'gatk_snpeff.gatk_snpeff_task.reference_tgz':'gs://broad-cil-devel-bucket/input_data/minion_illumina_hybrid_clean_MT.tgz',
    'gatk_snpeff.gatk_snpeff_task.output_disk_gb':'10',
    'gatk_snpeff.gatk_snpeff_task.debug_dump_flag':'onfail',
    'gatk_snpeff.gatk_snpeff_task.vcf_in':"gs://broad-cil-devel-bucket/gatk_joint_genotype/5929ccf3-d281-45ad-bf63-20d7b558184a/call-gatk_joint_genotype_task/test_cohort.vcf",
    'gatk_snpeff.gatk_snpeff_task.cohort_name':"test_cohort",
    
    }

(status,outputs) = run_wdl(wdl,inputs)

print(outputs)
