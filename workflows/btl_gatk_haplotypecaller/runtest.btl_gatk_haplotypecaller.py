from runcromwellremote import run_wdl

wdl = 'taskdef.btl_gatk_haplotypecaller.wdl'

inputs = {
    'gatk_haplotypecaller.gatk_haplotypecaller_task.reference_tgz':'gs://broad-cil-devel-bucket/input_data/minion_illumina_hybrid_clean_MT.tgz',
    'gatk_haplotypecaller.indelrealigner_bam':'gs://broad-cil-devel-bucket/input_data/Candida_Auris.tcir.bam',
    'gatk_haplotypecaller.indelrealigner_bam_index':'gs://broad-cil-devel-bucket/input_data/Candida_Auris.tcir.bam.bai',
    'gatk_haplotypecaller.gatk_haplotypecaller_task.output_disk_gb':'10',
    'gatk_haplotypecaller.gatk_haplotypecaller_task.sample_name':'Candida_Auris',
    'gatk_haplotypecaller.gatk_haplotypecaller_task.debug_dump_flag':'onfail',

    'gatk_haplotypecaller.gatk_haplotypecaller_task.known_sites_vcf_tbis':["gs://broad-cil-devel-bucket/input_data/7g8_gb4.combined.final.vcf.gz.tbi","gs://broad-cil-devel-bucket/input_data/hb3_dd2.combined.final.vcf.gz.tbi","gs://broad-cil-devel-bucket/input_data/3d7_hb3.combined.final.vcf.gz.tbi"],
    'gatk_haplotypecaller.gatk_haplotypecaller_task.known_sites_vcfs':["gs://broad-cil-devel-bucket/input_data/7g8_gb4.combined.final.vcf.gz","gs://broad-cil-devel-bucket/input_data/hb3_dd2.combined.final.vcf.gz","gs://broad-cil-devel-bucket/input_data/3d7_hb3.combined.final.vcf.gz"],

    }

(status,outputs) = run_wdl(wdl,inputs)

print(outputs)

