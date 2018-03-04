from runcromwellremote import run_wdl

wdl = 'taskdef.btl_gatk_joint_genotype.wdl'

inputs = {
    'gatk_joint_genotype.gatk_joint_genotype_task.reference_tgz':'gs://broad-cil-devel-bucket/input_data/minion_illumina_hybrid_clean_MT.tgz',
    'gatk_joint_genotype.gatk_joint_genotype_task.output_disk_gb':'10',
    'gatk_joint_genotype.gatk_joint_genotype_task.debug_dump_flag':'onfail',

    'gatk_joint_genotype.gatk_joint_genotype_task.known_sites_vcf_tbis':["gs://broad-cil-devel-bucket/input_data/7g8_gb4.combined.final.vcf.gz.tbi","gs://broad-cil-devel-bucket/input_data/hb3_dd2.combined.final.vcf.gz.tbi","gs://broad-cil-devel-bucket/input_data/3d7_hb3.combined.final.vcf.gz.tbi"],
    'gatk_joint_genotype.gatk_joint_genotype_task.known_sites_vcfs':["gs://broad-cil-devel-bucket/input_data/7g8_gb4.combined.final.vcf.gz","gs://broad-cil-devel-bucket/input_data/hb3_dd2.combined.final.vcf.gz","gs://broad-cil-devel-bucket/input_data/3d7_hb3.combined.final.vcf.gz"],

    'gatk_joint_genotype.gatk_joint_genotype_task.cohort_name':"test_cohort",
    'gatk_joint_genotype.gatk_joint_genotype_task.HaplotypeCaller_gvcfs':["gs://broad-cil-devel-bucket/gatk_haplotypecaller/9ad1d390-f507-44ae-9c70-72cb2ed9acdb/call-gatk_haplotypecaller_task/Candida_Auris.gvcf"],

    }

(status,outputs) = run_wdl(wdl,inputs)

print(outputs)

