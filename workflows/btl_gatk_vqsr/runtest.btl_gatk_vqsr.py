from runcromwellremote import run_wdl

wdl = 'taskdef.btl_gatk_vqsr.wdl'

inputs = {
    'gatk_vqsr.gatk_vqsr_task.reference_tgz':'gs://broad-cil-devel-bucket/input_data/minion_illumina_hybrid_clean_MT.tgz',
    'gatk_vqsr.gatk_vqsr_task.output_disk_gb':'10',
    'gatk_vqsr.gatk_vqsr_task.debug_dump_flag':'onfail',


    'gatk_vqsr.gatk_vqsr_task.cohort_name':"test_cohort",
    'gatk_vqsr.gatk_vqsr_task.indel_resource':["7g8_gb4,known=false,training=true,truth=true,prior=12.0 /gsap/garage-protistvector/U19_Aim4/Pf3K/7g8_gb4.combined.final.vcf.gz", "hb3_dd2,known=false,training=true,truth=true,prior=12.0 /gsap/garage-protistvector/U19_Aim4/Pf3K/hb3_dd2.combined.final.vcf.gz", "3d7_hb3,known=false,training=true,truth=true,prior=12.0 /gsap/garage-protistvector/U19_Aim4/Pf3K/3d7_hb3.combined.final.vcf.gz"],
    'gatk_vqsr.gatk_vqsr_task.snp_resource': ["7g8_gb4,known=false,training=true,truth=true,prior=15.0 /gsap/garage-protistvector/U19_Aim4/Pf3K/7g8_gb4.combined.final.vcf.gz", "hb3_dd2,known=false,training=true,truth=true,prior=15.0 /gsap/garage-protistvector/U19_Aim4/Pf3K/hb3_dd2.combined.final.vcf.gz", "3d7_hb3,known=false,training=true,truth=true,prior=15.0 /gsap/garage-protistvector/U19_Aim4/Pf3K/3d7_hb3.combined.final.vcf.gz"],
    'gatk_vqsr.gatk_vqsr_task.ts_filter_indel': 99.0,
    'gatk_vqsr.gatk_vqsr_task.snp_annotation': ["QD", "FS", "SOR", "DP", "MQ"],
    'gatk_vqsr.gatk_vqsr_task.genotype_caller_vcf':"gs://broad-cil-devel-bucket/gatk_joint_genotype/5929ccf3-d281-45ad-bf63-20d7b558184a/call-gatk_joint_genotype_task/test_cohort.vcf",
    'gatk_vqsr.gatk_vqsr_task.ts_filter_snp': 99.5,
    'gatk_vqsr.gatk_vqsr_task.indel_annotation': ["QD", "FS", "MQ"],
    'gatk_vqsr.gatk_vqsr_task.snp_max_gaussians': 2,
    'gatk_vqsr.gatk_vqsr_task.indel_max_gaussians': 2,
    'gatk_vqsr.gatk_vqsr_task.mq_cap_snp': 70,
    'gatk_vqsr.gatk_vqsr_task.mq_cap_indel': 70,
    

    'gatk_vqsr.gatk_vqsr_task.known_sites_vcf_tbis':["gs://broad-cil-devel-bucket/input_data/7g8_gb4.combined.final.vcf.gz.tbi","gs://broad-cil-devel-bucket/input_data/hb3_dd2.combined.final.vcf.gz.tbi","gs://broad-cil-devel-bucket/input_data/3d7_hb3.combined.final.vcf.gz.tbi"],
    'gatk_vqsr.gatk_vqsr_task.known_sites_vcfs':["gs://broad-cil-devel-bucket/input_data/7g8_gb4.combined.final.vcf.gz","gs://broad-cil-devel-bucket/input_data/hb3_dd2.combined.final.vcf.gz","gs://broad-cil-devel-bucket/input_data/3d7_hb3.combined.final.vcf.gz"],
    'gatk_vqsr.gatk_vqsr_task.snp_resource_params':["7g8_gb4,known=false,training=true,truth=true,prior=15.0", "hb3_dd2,known=false,training=true,truth=true,prior=15.0", "3d7_hb3,known=false,training=true,truth=true,prior=15.0"],

    'gatk_vqsr.gatk_vqsr_task.indel_resource_params':["7g8_gb4,known=false,training=true,truth=true,prior=12.0", "hb3_dd2,known=false,training=true,truth=true,prior=12.0", "3d7_hb3,known=false,training=true,truth=true,prior=12.0"],


    }

(status,outputs) = run_wdl(wdl,inputs)

print(outputs)

               
