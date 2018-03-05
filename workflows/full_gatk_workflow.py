from runcromwellremote import run_wdl


def check_run_wdl(wdl_path, inputs_dict):
    passed, outputs = run_wdl(wdl_path, inputs_dict)
    if not passed:
        sys.exit(1)
    return outputs

def run_index_reference(inputs):

    indexref_wdl_path = 'btl_gatk_indexref/taskdef.btl_gatk_indexref.wdl'
    indexref_inputs = {
        "gatk_indexref.IndexReference.ref_name": inputs['ref_name'],
        "gatk_indexref.IndexReference.output_disk_gb": "10",
        "gatk_indexref.IndexReference.ref_fasta": inputs['ref_fasta'],
        "gatk_indexref.IndexReference.debug_dump_flag": "onfail"
        }

    if False:
        indexref_outputs = check_run_wdl(indexref_wdl_path, indexref_inputs)
    else:
        indexref_outputs = {
            'gatk_indexref.IndexReference.reference_tgz':'gs://broad-cil-devel-bucket/gatk_indexref/5e5363b4-4ff5-4e49-a8eb-7db0218c4125/call-IndexReference/minion_illumina_hybrid_clean_MT.tgz'
        }

    outputs = {'reference_tgz':indexref_outputs['gatk_indexref.IndexReference.reference_tgz']}
    return outputs


def run_process_sample(inputs):


    alignbam_wdl_path = 'btl_gatk_alignbam/taskdef.btl_gatk_alignbam.wdl'
    alignbam_inputs = {
        'gatk_alignbam.gatk_alignbam_task.reference_tgz':inputs['reference_tgz'],
        'gatk_alignbam.gatk_alignbam_task.in_bam':inputs['in_bam'],
        'gatk_alignbam.gatk_alignbam_task.output_disk_gb':'10',
        'gatk_alignbam.gatk_alignbam_task.sample_name':inputs['sample_name'],
        'gatk_alignbam.gatk_alignbam_task.debug_dump_flag':'onfail',

        }
    if True:
        alignbam_outputs = check_run_wdl(alignbam_wdl_path, alignbam_inputs)
    else:
        alignbam_outputs = {
            'gatk_alignbam.gatk_alignbam_task.out_bam':'gs://broad-cil-devel-bucket/gatk_alignbam/1cd6a7fa-080f-4008-8e06-bc7eaad72b60/call-gatk_alignbam_task/Candida_Auris.bam',
            'gatk_alignbam.gatk_alignbam_task.out_bam_index':'gs://broad-cil-devel-bucket/gatk_alignbam/1cd6a7fa-080f-4008-8e06-bc7eaad72b60/call-gatk_alignbam_task/Candida_Auris.bam.bai'
        }



    tcir_wdl_path = 'btl_gatk_tcir/taskdef.btl_gatk_tcir.wdl'
    tcir_inputs = {
        'gatk_tcir.gatk_tcir_task.reference_tgz':inputs['reference_tgz'],
        'gatk_tcir.gatk_tcir_task.in_bam':alignbam_outputs['gatk_alignbam.gatk_alignbam_task.out_bam'],
        'gatk_tcir.gatk_tcir_task.in_bam_index':alignbam_outputs['gatk_alignbam.gatk_alignbam_task.out_bam_index'],
        'gatk_tcir.gatk_tcir_task.output_disk_gb':'10',
        'gatk_tcir.gatk_tcir_task.sample_name':inputs['sample_name'],
        'gatk_tcir.gatk_tcir_task.debug_dump_flag':'onfail',

        }

    if True:
        tcir_outputs = check_run_wdl(tcir_wdl_path, tcir_inputs)
    else:
        tcir_outputs = {
            'gatk_tcir.gatk_tcir_task.out_bam':"gs://broad-cil-devel-bucket/gatk_tcir/235faed4-cc6b-4d79-a213-843af2f3a58f/call-gatk_tcir_task/Candida_Auris.tcir.bam",
            'gatk_tcir.gatk_tcir_task.out_bam_index':'gs://broad-cil-devel-bucket/gatk_tcir/235faed4-cc6b-4d79-a213-843af2f3a58f/call-gatk_tcir_task/Candida_Auris.tcir.bam.bai'
        }




    bqsr_wdl_path = 'btl_gatk_bqsr/taskdef.btl_gatk_bqsr.wdl'
    bqsr_inputs = {
        'gatk_bqsr.gatk_bqsr_task.in_bam':tcir_outputs['gatk_tcir.gatk_tcir_task.out_bam'],
        'gatk_bqsr.gatk_bqsr_task.in_bam_index':tcir_outputs['gatk_tcir.gatk_tcir_task.out_bam_index'],
        'gatk_bqsr.gatk_bqsr_task.debug_dump_flag':'onfail',
        'gatk_bqsr.gatk_bqsr_task.output_disk_gb':'10',
        'gatk_bqsr.gatk_bqsr_task.reference_tgz':inputs['reference_tgz'],
        'gatk_bqsr.gatk_bqsr_task.known_sites_vcfs':["gs://broad-cil-devel-bucket/input_data/7g8_gb4.combined.final.vcf.gz","gs://broad-cil-devel-bucket/input_data/hb3_dd2.combined.final.vcf.gz","gs://broad-cil-devel-bucket/input_data/3d7_hb3.combined.final.vcf.gz"],
        'gatk_bqsr.gatk_bqsr_task.known_sites_vcf_tbis':["gs://broad-cil-devel-bucket/input_data/7g8_gb4.combined.final.vcf.gz.tbi","gs://broad-cil-devel-bucket/input_data/hb3_dd2.combined.final.vcf.gz.tbi","gs://broad-cil-devel-bucket/input_data/3d7_hb3.combined.final.vcf.gz.tbi"],
        'gatk_bqsr.gatk_bqsr_task.sample_name':inputs['sample_name'],
    }
    if True:
        bqsr_outputs = check_run_wdl(bqsr_wdl_path, bqsr_inputs)
    else:
        bqsr_outpus = {
            'gatk_bqsr.gatk_bqsr_task.out_bam':'gs://broad-cil-devel-bucket/gatk_bqsr/aa98f8e0-7a8e-4df6-9700-d88b1b870875/call-gatk_bqsr_task/Candida_Auris.bqsr.bam',
            'gatk_bqsr.gatk_bqsr_task.out_bam_index':'gs://broad-cil-devel-bucket/gatk_bqsr/aa98f8e0-7a8e-4df6-9700-d88b1b870875/call-gatk_bqsr_task/Candida_Auris.bqsr.bam.bai',
            'gatk_bqsr.gatk_bqsr_task.out_bqsr_table':'gs://broad-cil-devel-bucket/gatk_bqsr/aa98f8e0-7a8e-4df6-9700-d88b1b870875/call-gatk_bqsr_task/Candida_Auris.bqsr.table',
        }

    haplotypecaller_wdl_path = 'btl_gatk_haplotypecaller/taskdef.btl_gatk_haplotypecaller.wdl'
    haplotypecaller_inputs = {
        'gatk_haplotypecaller.gatk_haplotypecaller_task.reference_tgz':inputs['reference_tgz'],
        'gatk_haplotypecaller.gatk_haplotypecaller_task.in_bam':bqsr_outputs['gatk_bqsr.gatk_bqsr_task.out_bam'],
        'gatk_haplotypecaller.gatk_haplotypecaller_task.in_bam_index':bqsr_outputs['gatk_bqsr.gatk_bqsr_task.out_bam_index'],
        'gatk_haplotypecaller.gatk_haplotypecaller_task.bqsr_table':bqsr_outputs['gatk_bqsr.gatk_bqsr_task.out_bqsr_table']
        'gatk_haplotypecaller.gatk_haplotypecaller_task.output_disk_gb':'10',
        'gatk_haplotypecaller.gatk_haplotypecaller_task.sample_name':inputs['sample_name'],
        'gatk_haplotypecaller.gatk_haplotypecaller_task.debug_dump_flag':'onfail',

        'gatk_haplotypecaller.gatk_haplotypecaller_task.known_sites_vcf_tbis':["gs://broad-cil-devel-bucket/input_data/7g8_gb4.combined.final.vcf.gz.tbi","gs://broad-cil-devel-bucket/input_data/hb3_dd2.combined.final.vcf.gz.tbi","gs://broad-cil-devel-bucket/input_data/3d7_hb3.combined.final.vcf.gz.tbi"],
        'gatk_haplotypecaller.gatk_haplotypecaller_task.known_sites_vcfs':["gs://broad-cil-devel-bucket/input_data/7g8_gb4.combined.final.vcf.gz","gs://broad-cil-devel-bucket/input_data/hb3_dd2.combined.final.vcf.gz","gs://broad-cil-devel-bucket/input_data/3d7_hb3.combined.final.vcf.gz"],

        }
    if True:
        haplotypecaller_outputs = check_run_wdl(haplotypecaller_wdl_path, haplotypecaller_inputs)
    else:
        haplotypecaller_outputs = {
            'gatk_haplotypecaller.gatk_haplotypecaller_task.out_gvcf':'gs://broad-cil-devel-bucket/gatk_haplotypecaller/dbf79433-58e5-4931-b77c-3fcf801f3db2/call-gatk_haplotypecaller_task/Candida_Auris.gvcf'
        }


    outputs = {'out_gvcf':haplotypecaller_outputs['gatk_haplotypecaller.gatk_haplotypecaller_task.out_gvcf']}
    return outputs

def run_process_cohort(inputs):

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
#need to add R to docker to enable plots

if __name__ == '__main__':
    index_reference_inputs = {
        'ref_fasta':'gs://broad-cil-devel-bucket/input_data/minion_illumina_hybrid_clean_MT.fasta',
        'ref_name':'minion_illumina_hybrid_clean_MT'
        }
    if False:
        index_reference_outputs = run_index_reference(index_reference_inputs)
    else:
        index_reference_outputs = {
            'reference_tgz':'gs://broad-cil-devel-bucket/gatk_indexref/5e5363b4-4ff5-4e49-a8eb-7db0218c4125/call-IndexReference/minion_illumina_hybrid_clean_MT.tgz'
        }



    input_bams_by_sample_name={
        'Candida_Auris':'gs://broad-cil-devel-bucket/input_data/Candida_Auris.bam'
    }

    for sample_name in input_bams_by_sample_name:
        bam = input_bams_by_sample_name[sample_name]
        process_sample_inputs = {
            'reference_tgz':index_reference_outputs['reference_tgz'],
            'in_bam':bam,
            'sample_name':sample_name
        }
        process_sample_outputs = run_process_sample(process_sample_inputs)