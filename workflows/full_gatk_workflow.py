from runcromwellremote import run_wdl


def check_run_wdl(wdl_path, inputs_dict)
    passed, outputs = run_wdl(wdl_path, inputs_dict)
    if not passed:
        sys.exit(1)
    return outputs

def index_reference(inputs)

    indexref_wdl_path = 'btl_gatk_indexref/taskdef.btl_gatk_indexref.wdl'
    indexref_inputs = {
        "gatk_indexref.IndexReference.ref_name": "minion_illumina_hybrid_clean_MT",
        "gatk_indexref.IndexReference.output_disk_gb": "10",
        "gatk_indexref.IndexReference.ref_fasta": "gs://broad-cil-devel-bucket/input_data/minion_illumina_hybrid_clean_MT.fasta",
        "gatk_indexref.IndexReference.debug_dump_flag": "onfail"
        }
    indexref_outputs = check_run_wdl(indexref_wdl_path, indexref_inputs)

    outputs = {'reference_tgz':indexref_outputs['reference_tgz']}
    return outputs


def process_sample(inputs)


    alignbam_wdl_path = 'btl_gatk_alignbam/taskdef.btl_gatk_alignbam.wdl'
    alignbam_inputs = {
        'gatk_alignbam.gatk_alignbam_task.reference_tgz':inputs['reference_tgz'],
        'gatk_alignbam.gatk_alignbam_task.in_bam':inputs['in_bam'],
        'gatk_alignbam.gatk_alignbam_task.in_bam_index':inputs['in_bam_index'],
        'gatk_alignbam.gatk_alignbam_task.output_disk_gb':'10',
        'gatk_alignbam.gatk_alignbam_task.sample_name':inputs['sample_name'],
        'gatk_alignbam.gatk_alignbam_task.debug_dump_flag':'onfail',

        }
    alignbam_outputs = check_run_wdl(alignbam_wdl_path, alignbam_inputs)



    tcir_wdl_path = 'btl_gatk_tcir/taskdef.btl_gatk_tcir.wdl'
    tcir_inputs = {
        'gatk_tcir.gatk_tcir_task.reference_tgz':inputs['reference_tgz'],
        'gatk_tcir.gatk_tcir_task.in_bam':alignbam_outputs['out_bam'],
        'gatk_tcir.gatk_tcir_task.in_bam_index':alignbam_outputs['out_bam_index'],
        'gatk_tcir.gatk_tcir_task.output_disk_gb':'10',
        'gatk_tcir.gatk_tcir_task.sample_name':inputs['sample_name'],
        'gatk_tcir.gatk_tcir_task.debug_dump_flag':'onfail',

        }
    tcir_outputs = check_run_wdl(tcir_wdl_path, tcir_inputs)




    bqsr_wdl_path = 'btl_gatk_bqsr/taskdef.btl_gatk_bqsr.wdl'
    bqsr_inputs = {
        'gatk_bqsr.indelrealigner_bam':tcir_outputs['out_bam'],
        'gatk_bqsr.indelrealigner_bam_index':tcir_outputs['out_bam_index'],
        'gatk_bqsr.gatk_bqsr_task.debug_dump_flag':'onfail',
        'gatk_bqsr.gatk_bqsr_task.output_disk_gb':'10',
        'gatk_bqsr.gatk_bqsr_task.reference_tgz':inputs['reference_tgz'],
        'gatk_bqsr.gatk_bqsr_task.known_sites_vcfs':["gs://broad-cil-devel-bucket/input_data/7g8_gb4.combined.final.vcf.gz","gs://broad-cil-devel-bucket/input_data/hb3_dd2.combined.final.vcf.gz","gs://broad-cil-devel-bucket/input_data/3d7_hb3.combined.final.vcf.gz"],
        'gatk_bqsr.gatk_bqsr_task.known_sites_vcf_tbis':["gs://broad-cil-devel-bucket/input_data/7g8_gb4.combined.final.vcf.gz.tbi","gs://broad-cil-devel-bucket/input_data/hb3_dd2.combined.final.vcf.gz.tbi","gs://broad-cil-devel-bucket/input_data/3d7_hb3.combined.final.vcf.gz.tbi"],
        'gatk_bqsr.gatk_bqsr_task.sample_name':inputs['sample_name'],
    }
    bqsr_outputs = check_run_wdl(bqsr_wdl_path, bqsr_inputs)



    haplotypecaller_wdl_path = 'btl_gatk_haplotypecaller/taskdef.btl_gatk_haplotypecaller.wdl'
    haplotypecaller_inputs = {
        'gatk_haplotypecaller.gatk_haplotypecaller_task.reference_tgz':inputs['reference_tgz'],
        'gatk_haplotypecaller.indelrealigner_bam':bqsr_outputs['out_bam']
        'gatk_haplotypecaller.indelrealigner_bam_index':bqsr_outputs['out_bam_index']
        'gatk_haplotypecaller.gatk_haplotypecaller_task.output_disk_gb':'10',
        'gatk_haplotypecaller.gatk_haplotypecaller_task.sample_name':inputs['sample_name'],
        'gatk_haplotypecaller.gatk_haplotypecaller_task.debug_dump_flag':'onfail',

        'gatk_haplotypecaller.gatk_haplotypecaller_task.known_sites_vcf_tbis':["gs://broad-cil-devel-bucket/input_data/7g8_gb4.combined.final.vcf.gz.tbi","gs://broad-cil-devel-bucket/input_data/hb3_dd2.combined.final.vcf.gz.tbi","gs://broad-cil-devel-bucket/input_data/3d7_hb3.combined.final.vcf.gz.tbi"],
        'gatk_haplotypecaller.gatk_haplotypecaller_task.known_sites_vcfs':["gs://broad-cil-devel-bucket/input_data/7g8_gb4.combined.final.vcf.gz","gs://broad-cil-devel-bucket/input_data/hb3_dd2.combined.final.vcf.gz","gs://broad-cil-devel-bucket/input_data/3d7_hb3.combined.final.vcf.gz"],

        }
    haplotypecaller_outputs = check_run_wdl(haplotypecaller_wdl_path, haplotypecaller_inputs)


    outputs = {'out_gvcf':haplotypecaller_outputs['out_gvcf']}
    return outputs

def process_cohort(inputs):

    joint_genotype_wdl_path = 'btl_gatk_joint_genotype/taskdef.btl_gatk_joint_genotype.wdl'
    joint_genotype_inputs = {
        'gatk_joint_genotype.gatk_joint_genotype_task.reference_tgz':inputs['reference_tgz'],
        'gatk_joint_genotype.gatk_joint_genotype_task.output_disk_gb':'10',
        'gatk_joint_genotype.gatk_joint_genotype_task.debug_dump_flag':'onfail',

        'gatk_joint_genotype.gatk_joint_genotype_task.known_sites_vcf_tbis':["gs://broad-cil-devel-bucket/input_data/7g8_gb4.combined.final.vcf.gz.tbi","gs://broad-cil-devel-bucket/input_data/hb3_dd2.combined.final.vcf.gz.tbi","gs://broad-cil-devel-bucket/input_data/3d7_hb3.combined.final.vcf.gz.tbi"],
        'gatk_joint_genotype.gatk_joint_genotype_task.known_sites_vcfs':["gs://broad-cil-devel-bucket/input_data/7g8_gb4.combined.final.vcf.gz","gs://broad-cil-devel-bucket/input_data/hb3_dd2.combined.final.vcf.gz","gs://broad-cil-devel-bucket/input_data/3d7_hb3.combined.final.vcf.gz"],

        'gatk_joint_genotype.gatk_joint_genotype_task.cohort_name':inputs['cohort_name'],
        'gatk_joint_genotype.gatk_joint_genotype_task.HaplotypeCaller_gvcfs':inputs['gvcf_list'],

        }
    joint_genotype_outputs = check_run_wdl(joint_genotype_wdl_path, joint_genotype_inputs)


    vqsr_wdl_path = 'btl_gatk_vqsr/taskdef.btl_gatk_vqsr.wdl'

    vqsr_inputs = {
        'gatk_vqsr.gatk_vqsr_task.reference_tgz':inputs['reference_tgz'],
        'gatk_vqsr.gatk_vqsr_task.output_disk_gb':'10',
        'gatk_vqsr.gatk_vqsr_task.debug_dump_flag':'onfail',
        'gatk_vqsr.gatk_vqsr_task.genotype_caller_vcf':joint_genotype_outputs['vcf_out'],
        'gatk_vqsr.gatk_vqsr_task.cohort_name':inputs['cohort_name'],

        'gatk_vqsr.gatk_vqsr_task.ts_filter_indel': 99.0,
        'gatk_vqsr.gatk_vqsr_task.snp_annotation': ["QD", "FS", "SOR", "DP", "MQ"],
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
    vqsr_outputs = check_run_wdl(vqsr_outputs_wdl_path, vqsr_inputs)

    variant_filtration_wdl = 'btl_gatk_variant_filtration/taskdef.btl_gatk_variant_filtration.wdl'
#        'gatk_variant_filtration.genotype_caller_vcf':vqsr_outputs['out_vcf'],
    variant_filtration_inputs = {
        'gatk_variant_filtration.gatk_variant_filtration_task.reference_tgz':inputs['reference_tgz'],
        'gatk_variant_filtration.sv_vcf':joint_genotype_outputs['out_vcf'],
        'gatk_variant_filtration.gatk_variant_filtration_task.output_disk_gb':'10',
        'gatk_variant_filtration.gatk_variant_filtration_task.debug_dump_flag':'onfail',
        'gatk_variant_filtration.gatk_variant_filtration_task.cohort_name':inputs['cohort_name'],
        'gatk_variant_filtration.gatk_variant_filtration_task.snp_filter_expression':"VQSLOD <= 0.0",
        'gatk_variant_filtration.gatk_variant_filtration_task.indel_filter_expression':"VQSLOD <= 0.0",
        }

    variant_filtration_outputs = check_run_wdl(variant_filtration_outputs_wdl_path, variant_filtration_inputs)

#need to add snpeff
