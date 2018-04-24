#!/usr/bin/python

import os
import sys
import collections





import pipetteClient

def run_pcawg_pipeline(communicationDirBase, pipelineOutDir, pipelinePriority, adaptor_dir, module_dir, refdata_dir, tmp_dir, gnos_outdir, broad_outdir, indiv_id, bam_tumor, bam_normal, sub_workflow, fail_mode):

    #count number of CPU's
    fid = open('/proc/cpuinfo')
    cpu_count = 0
    for line in fid:
        if 'processor' in line:
            cpu_count = cpu_count + 1
    fid.close()

    pipeline = initialize_pipeline(communicationDirBase, pipelineOutDir, pipelinePriority, fail_mode)
    
    annotations = pcawg_pipeline_pp(pipeline, adaptor_dir, module_dir, refdata_dir, tmp_dir, gnos_outdir, broad_outdir, indiv_id, bam_tumor, bam_normal, sub_workflow, cpu_count)
    
    pipeline.go()

    #add to annotation variable
    #dump annotation variable to disk    
    
    
def pcawg_pipeline_pp (pipeline, adaptor_dir, module_dir, refdata_dir, tmp_dir, gnos_outdir, broad_outdir,
                       indiv_id, bam_tumor, bam_normal, sub_workflow, cpu_count):

    # if False:
    #     #omit for the initial pipeline
    #     module_subdir = 'coclean'
    #     (bam_cleaned_tumor, bam_cleaned_normal) = coclean_pp(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, indiv_id, bam_tumor, bam_normal)
    #

    # if True:
    #     module_subdir = 'fix_header_t'
    #     bam_tumor = fix_header_pp(pipeline,module_subdir,adaptor_dir, module_dir, refdata_dir,  indiv_id, bam_tumor_arg)
    #
    #     module_subdir = 'fix_header_n'
    #     bam_normal = fix_header_pp(pipeline,module_subdir,adaptor_dir, module_dir, refdata_dir,  indiv_id, bam_normal_arg)
    # else:
    #     bam_tumor = bam_tumor_arg
    #     bam_normal = bam_normal_arg

    sub_workflow_set = set(sub_workflow.split(','))
    use_pcawg_contigs = False
    # OXoq
    # Mutect
    # ContEst
    # ReCapseg coverage
    # collapseHetSitesIntervals
    # genotypeGVCF

    if sub_workflow_set.intersection({'hello'}):

        module_subdir = 'hello'
        hello_outfile = hello(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, indiv_id)
        make_external_links_pm(pipeline,'hello_extlinks',module_dir, broad_outdir, [hello_outfile] )


    if sub_workflow_set.intersection({'contest','contest_mutect','contest_mutect_variantbam','contest_M2','all'}):
        module_subdir = 'contest'
        contest_percent_value_file = contest_pp(pipeline, module_subdir,adaptor_dir, module_dir, refdata_dir,  indiv_id, bam_tumor, bam_normal, tmp_dir, use_pcawg_contigs)

        make_external_links_pm(pipeline,'contest_extlinks',module_dir, broad_outdir, [contest_percent_value_file] )

        module_subdir = 'contest_div100'
        contest_fraction_value_file = div100_pm(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, indiv_id, contest_percent_value_file)

    else:
        contest_percent_value_file = None
        contest_fraction_value_file = None
        
    if sub_workflow_set.intersection({'nocont_mutect','contest_mutect','contest_mutect_variantbam','all'}):


        #tbd ensure region is set to all for production, not test.
        module_subdir = 'mutect'
        region = 'all'
        (call_stats, coverage, power) = mutect1_sg_pp(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, indiv_id, bam_tumor, bam_normal, contest_fraction_value_file,region, use_pcawg_contigs)
        make_external_links_pm(pipeline, 'mutect_extlinks', module_dir, broad_outdir, [call_stats, coverage, power] )


        module_subdir = 'callstats_to_maflite'
        maflite_file = callstats_to_maflite(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, indiv_id, call_stats)

        module_subdir = 'oncotator'
        mutect_vcf_output = oncotator(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, indiv_id, maflite_file)

        module_subdir = 'tabix_mutect'
        tabix_mutect_outputs = tabix_pm(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, indiv_id, mutect_vcf_output,
                                        'broad-mutect.DATECODE.somatic.snv_mnv')


        make_external_links_pm(pipeline,'mutect_tabix_extlinks',module_dir, gnos_outdir, tabix_mutect_outputs )

    else:
        call_stats = None

    if sub_workflow_set.intersection({'contest_M2','all'}):
     
        #tbd ensure region is set to all for production, not test.
        module_subdir = 'M2_scatter'
        region = 'all'
        (M2_vcf,M2_pass_vcf) = mutect2_sg_pp(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, indiv_id, bam_tumor, bam_normal, contest_fraction_value_file,region, use_pcawg_contigs)

        #make_external_links_pm(pipeline, 'mutect2_extlinks', module_dir, broad_outdir, [M2_vcf,M2_pass_vcf] )

        module_subdir = 'tabix_mutect2'
        tabix_mutect2_outputs = tabix_pm(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, indiv_id, M2_pass_vcf,
                                        'broad-mutect2.DATECODE.somatic')

        make_external_links_pm(pipeline,'mutect2_tabix_extlinks',module_dir, gnos_outdir, tabix_mutect2_outputs )

    else:
        M2_vcf = None
        M2_pass_vcf = None

        
        
        
    if sub_workflow_set.intersection({'contest_mutect_variantbam','all'}):
        # TBD note new version available from Jeremiah
        module_subdir = 'variant_bam_tumor'
        variant_bam_tumor_outputs = variant_bam_pm(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, bam_tumor, call_stats,M2_pass_vcf)
        make_external_links_pm(pipeline,'variant_bam_tumor_extlinks',module_dir, broad_outdir, variant_bam_tumor_outputs )

        module_subdir = 'variant_bam_normal'
        variant_bam_normal_outputs = variant_bam_pm(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, bam_normal, call_stats,M2_pass_vcf)
        make_external_links_pm(pipeline,'variant_bam_normal_extlinks',module_dir, broad_outdir, variant_bam_normal_outputs )


    if sub_workflow_set.intersection({'svaba','dranger_svaba_mergesvcalls'}):

        if sub_workflow_set.intersection({'svaba@allcores'}):
            num_cores = str(cpu_count)
        else:
            num_cores = '3'

        module_subdir = 'svaba'
        (module_svaba_somatic_vcf_outfile,module_svaba_germline_vcf_outfile,module_svaba_somatic_indel_vcf_outfile,module_svaba_germline_indel_vcf_outfile,svaba_broad_files) = 	\
        	svaba_pm(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, bam_tumor, bam_normal, indiv_id, num_cores)

        module_subdir = 'tabix_svaba_somatic_sv'
        tabix_svaba_somatic_sv_outputs = tabix_pm(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, indiv_id,
                                                module_svaba_somatic_vcf_outfile, 'broad-svaba.DATECODE.somatic.sv')
        module_subdir = 'tabix_svaba_germline_sv'
        tabix_svaba_germline_sv_outputs = tabix_pm(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, indiv_id,
                                                module_svaba_germline_vcf_outfile, 'broad-svaba.DATECODE.germline.sv')
        module_subdir = 'tabix_svaba_somatic_indel'
        tabix_svaba_somatic_indel_outputs = tabix_pm(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, indiv_id,
                                                module_svaba_somatic_indel_vcf_outfile, 'broad-svaba.DATECODE.somatic.indel')
        module_subdir = 'tabix_svaba_germline_indel'
        tabix_svaba_germline_indel_outputs = tabix_pm(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, indiv_id,
                                                module_svaba_germline_indel_vcf_outfile, 'broad-svaba.DATECODE.germline.indel')

        svaba_tabixed_vcfs = tabix_svaba_somatic_sv_outputs+tabix_svaba_germline_sv_outputs+ \
                               tabix_svaba_somatic_indel_outputs+tabix_svaba_germline_indel_outputs

        make_external_links_pm(pipeline,'svaba_vcf_extlinks',module_dir, gnos_outdir, svaba_tabixed_vcfs)
        make_external_links_pm(pipeline,'svaba_broad_extlinks',module_dir, broad_outdir, svaba_broad_files)

    # if sub_workflow_set.intersection({'svaba','dranger_svaba_mergesvcalls'}):
    #     # TBD note new version available from Jeremiah
    #     #TBD - the vcf files need to be tabixed still

    #     if sub_workflow_set.intersection({'svaba@allcores'}):
    #         num_cores = str(cpu_count)
    #     else:
    #         num_cores = '3'

    #     module_subdir = 'svaba'
    #     (module_svaba_somatic_vcf_outfile,module_svaba_germline_vcf_outfile,module_svaba_somatic_indel_vcf_outfile,module_svaba_germline_indel_vcf_outfile,svaba_broad_files) = 	\
    #     	svaba_pm(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, bam_tumor, bam_normal, indiv_id, num_cores)

    #     module_subdir = 'tabix_svaba_somatic_sv'
    #     tabix_svaba_somatic_sv_outputs = tabix_pm(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, indiv_id,
    #                                             module_svaba_somatic_vcf_outfile, 'broad-svaba.DATECODE.somatic.sv')
    #     module_subdir = 'tabix_svaba_germline_sv'
    #     tabix_svaba_germline_sv_outputs = tabix_pm(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, indiv_id,
    #                                             module_svaba_germline_vcf_outfile, 'broad-svaba.DATECODE.germline.sv')
    #     module_subdir = 'tabix_svaba_somatic_indel'
    #     tabix_svaba_somatic_indel_outputs = tabix_pm(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, indiv_id,
    #                                             module_svaba_somatic_indel_vcf_outfile, 'broad-svaba.DATECODE.somatic.indel')
    #     module_subdir = 'tabix_svaba_germline_indel'
    #     tabix_svaba_germline_indel_outputs = tabix_pm(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, indiv_id,
    #                                             module_svaba_germline_indel_vcf_outfile, 'broad-svaba.DATECODE.germline.indel')

    #     svaba_tabixed_vcfs = tabix_svaba_somatic_sv_outputs+tabix_svaba_germline_sv_outputs+ \
    #                            tabix_svaba_somatic_indel_outputs+tabix_svaba_germline_indel_outputs

    #     make_external_links_pm(pipeline,'svaba_vcf_extlinks',module_dir, gnos_outdir, svaba_tabixed_vcfs)
    #     make_external_links_pm(pipeline,'svaba_broad_extlinks',module_dir, broad_outdir, svaba_broad_files)

    if sub_workflow_set.intersection({'fragcounter','all'}):
        module_subdir = 'fragcounter_tumor'
        JaBbA_cov_rds_tumor = fragCounter_pm(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, bam_tumor)
        make_external_links_pm(pipeline,'fragcounter_tumor_extlinks',module_dir, broad_outdir, [JaBbA_cov_rds_tumor])

        module_subdir = 'fragcounter_normal'
        JaBbA_cov_rds_normal = fragCounter_pm(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, bam_normal)
        make_external_links_pm(pipeline,'fragcounter_normal_extlinks',module_dir, broad_outdir, [JaBbA_cov_rds_normal])


    if sub_workflow_set.intersection({'oxoq','all'}):
        module_subdir = 'oxoq'
        oxoq_file = OxoQ_pm(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir,  indiv_id, bam_tumor, bam_normal, use_pcawg_contigs)
        make_external_links_pm(pipeline,'oxoq_extlinks',module_dir, broad_outdir, [oxoq_file])

        
    if sub_workflow_set.intersection({'tokens','all'}):
        module_subdir = 'tokens'
        pileup_tokens_outfile = tokens_pm(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, indiv_id, bam_tumor, bam_normal)
        make_external_links_pm(pipeline,'tokens_extlinks',module_dir, broad_outdir, [pileup_tokens_outfile])



    if sub_workflow_set.intersection({'haplotypecaller','haplotypecaller_forcecallhets','all'}):
        module_subdir = 'extract_bam_id'
        bam_id_file = extract_bam_id_pp(pipeline, module_subdir, indiv_id, adaptor_dir, module_dir,refdata_dir, bam_normal)
        make_external_links_pm(pipeline,'extract_bam_id_extlinks',module_dir, broad_outdir, [bam_id_file])

        module_subdir = 'haplotypecaller_sg'
        (haplotype_caller_gvcf, haplotype_caller_tbi) = haplotype_caller_sg_pp(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, indiv_id, bam_normal)
        #The vcf files from haplotype caller are not ready for PCAWG consumption, they need further processing via VQSR
        make_external_links_pm(pipeline, 'haplotypecaller_sg_extlinks' ,module_dir, broad_outdir, (haplotype_caller_gvcf, haplotype_caller_tbi))


    if sub_workflow_set.intersection({'haplotypecaller_forcecallhets'}):  # removed from 'all'
        module_subdir = 'genotype_gvcf'
        single_sample_gvcf_filename = haplotype_caller_gvcf
        single_sample_vcf_filename = genotypeGVCF(pipeline, module_dir, module_subdir, indiv_id, single_sample_gvcf_filename, use_pcawg_contigs)
        make_external_links_pm(pipeline, 'genotype_gvcf_extlinks', module_dir, broad_outdir, [single_sample_vcf_filename])

        module_subdir = 'collapse_het_sites_intervals'
        force_call_het_sites_intervals_filename = collapseHetSitesIntervals(pipeline, module_dir, module_subdir, indiv_id, call_stats, single_sample_vcf_filename, use_pcawg_contigs)
        make_external_links_pm(pipeline, 'collapse_het_sites_intervals_extlinks', module_dir, broad_outdir, [force_call_het_sites_intervals_filename])

        # Module is needed
        module_subdir = 'mutect_het_sites'
        region = 'force_call'
        (callstats_het_sites, coverage_het_sites) = mutect_sg_pp(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, indiv_id, bam_tumor, bam_normal,
        contest_value_file, region, use_pcawg_contigs, force_call_het_sites_intervals_filename)
        make_external_links_pm(pipeline, 'mutect_het_sites_extlinks', module_dir, broad_outdir, [callstats_het_sites, coverage_het_sites])

    if sub_workflow_set.intersection({'recapseg','all'}):
        chroms = map(str, range(1, 23)+['X', 'Y'])
        re_capseg_cov_chrom_normals = []
        re_capseg_cov_chrom_tumors = []
        for chrom in chroms:
            module_subdir = 're_capseg_coverage_normal/chr%s' % chrom
            indiv_id_suffix = 'normal_chr%s' % chrom
            re_capseg_cov_chrom_normal = ReCapSeg_Coverage_pp(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, indiv_id, bam_normal, chrom, indiv_id_suffix, use_pcawg_contigs)
            re_capseg_cov_chrom_normals += [re_capseg_cov_chrom_normal]

            module_subdir = 're_capseg_coverage_tumor/chr%s' % chrom
            indiv_id_suffix = 'tumor_chr%s' % chrom
            re_capseg_cov_chrom_tumor = ReCapSeg_Coverage_pp(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, indiv_id, bam_tumor, chrom, indiv_id_suffix, use_pcawg_contigs)
            re_capseg_cov_chrom_tumors += [re_capseg_cov_chrom_tumor]

        module_subdir = 're_capseg_coverage_normal/merged'
        indiv_id_suffix = 'normal'
        re_capseg_cov_normal = concatenate_coverage_files_pm(pipeline, module_subdir, indiv_id, indiv_id_suffix, re_capseg_cov_chrom_normals)
        make_external_links_pm(pipeline,'re_capseg_coverage_normal_extlinks', module_dir, broad_outdir, [re_capseg_cov_normal])


        module_subdir = 're_capseg_coverage_tumor/merged'
        indiv_id_suffix = 'tumor'
        re_capseg_cov_tumor = concatenate_coverage_files_pm(pipeline, module_subdir, indiv_id, indiv_id_suffix, re_capseg_cov_chrom_tumors)

        make_external_links_pm(pipeline,'re_capseg_coverage_tumor_extlinks', module_dir, broad_outdir, [re_capseg_cov_tumor])


    if sub_workflow_set.intersection({'dranger','dranger_svaba_mergesvcalls','all'}):
        module_subdir = 'dRangerPreProcess_Tumor'
        (tdir,tisz) = dRangerPreprocess_pp(pipeline, module_subdir,adaptor_dir, module_dir, refdata_dir, indiv_id, bam_tumor, bam_type="Tumor", save_intermediate_files="True")
        make_external_links_pm(pipeline,'dRangerPreProcess_Tumor_extlinks', module_dir, broad_outdir, [tisz])

        module_subdir = 'dRangerPreProcess_Normal'
        (ndir,nisz) = dRangerPreprocess_pp(pipeline, module_subdir,adaptor_dir, module_dir, refdata_dir, indiv_id, bam_normal, bam_type="Normal", save_intermediate_files="True")
        make_external_links_pm(pipeline,'dRangerPreProcess_Normal_extlinks', module_dir, broad_outdir, [nisz])

        module_subdir = 'dRangerRun'
        drrbp,drrmat,stdout,stderr = dRangerRun_pp(pipeline, module_subdir,adaptor_dir, module_dir, refdata_dir, indiv_id, bam_tumor, bam_normal,tdir,ndir,save_intermediate_files='true',use_PoN_JH='true')
        make_external_links_pm(pipeline,'dRangerRun_extlinks', module_dir, broad_outdir, [drrbp,drrmat,stdout,stderr])

        module_subdir = 'BreakPointer_Tumor'
        (tbp,tsam) = dRanger_Breakpointer_pp(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, indiv_id,bam_tumor,drrbp,tisz)
        make_external_links_pm(pipeline,'BreakPointer_Tumor_extlinks', module_dir, broad_outdir, [tbp,tsam])

        module_subdir = 'BreakPointer_Normal'
        (nbp,nsam) = dRanger_Breakpointer_pp(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, indiv_id,bam_normal,drrbp,nisz)
        make_external_links_pm(pipeline,'BreakPointer_Normal_extlinks', module_dir, broad_outdir, [nbp,nsam])

        module_subdir = 'dRanger_Finalize'
        (drrfinalized,drrmatfinalized,da,ds) = dRangerFinalize_pp(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, indiv_id,tbp,nbp,drrmat)
        make_external_links_pm(pipeline,'dRanger_Finalize_extlinks', module_dir, broad_outdir, [drrfinalized,drrmatfinalized,da,ds])

        module_subdir = 'getdRangerSupportingReads_Tumor'
        dRanger_reads = getdRangerSupportingReads_pp(pipeline, module_subdir, adaptor_dir, module_dir,indiv_id,bam_tumor,tsam,drrmatfinalized,drrfinalized)
        make_external_links_pm(pipeline,'getdRangerSupportingReads_Tumor_extlinks', module_dir, broad_outdir, [dRanger_reads])

        module_subdir = 'getdRangerSupportingReads_Normal'
        dRanger_reads_normal = getdRangerSupportingReads_pp(pipeline, module_subdir, adaptor_dir, module_dir,indiv_id,bam_normal,nsam,drrmatfinalized,drrfinalized)
        make_external_links_pm(pipeline,'getdRangerSupportingReads_Normal_extlinks', module_dir, broad_outdir, [dRanger_reads_normal])

        module_subdir = 'dRanger2VCF'
        dRangerVCF,dranger_detail = dRanger2VCF_pp(pipeline, module_subdir, adaptor_dir, module_dir, drrmatfinalized, refdata_dir, indiv_id,tsam,dRanger_reads,dRanger_reads_normal)
        make_external_links_pm(pipeline,'dRanger2VCF_extlinks', module_dir, broad_outdir, [dranger_detail])

        module_subdir = 'tabix_dRanger'
        tabix_dranger_outputs = tabix_pm(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, indiv_id, dRangerVCF,
                                        'broad-dRanger.DATECODE.somatic.sv')
        make_external_links_pm(pipeline,'tabix_dRanger_extlinks',module_dir, gnos_outdir, tabix_dranger_outputs)

    if sub_workflow_set.intersection({'dranger_svaba_mergesvcalls'}):
        module_subdir = 'merge_sv_vcf'
        sv_merged_vcf = sv_merge_pm(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, module_svaba_somatic_vcf_outfile, dRangerVCF, indiv_id)

        module_subdir = 'tabix_merge_sv_vcf'        
        tabix_sv_merge_outputs = tabix_pm(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, indiv_id, sv_merged_vcf,
                                          'broad-dRanger_svaba.DATECODE.somatic.sv')
        make_external_links_pm(pipeline,'tabix_merge_sv_vcf_extlinks',module_dir, gnos_outdir, tabix_sv_merge_outputs)

    annotations = {}
    return annotations


def hello(pipeline,module_subdir,adaptor_dir, module_dir, refdata_dir,  indiv_id):

    out_fn = 'outfile.txt'
    out_module_path = os.path.join('$MODULEOUTDIR',out_fn)

    cmdStr = " ".join(['python', os.path.join(module_dir,'hello','hello.py'), indiv_id])



    moduleSubDir = module_subdir
    resources = {'maxmem':1, 'maxtime':3600}
    jobName = 'hello'
    inputFiles = []
    filesToBeOutput = [out_module_path]
    filesToBeDeleted = []

    filesToBeOutput_expanded = pipeline.dispense(
        moduleSubDir = moduleSubDir, #unique subdirectory under PIPELINEOUTDIR
        cmdStr = cmdStr, #command-line to run
        resources = resources, #dict with maxmem (in GB) and maxtime (in secs)
        jobName= jobName, #non-unique job-type name
        inputFiles=inputFiles, #list of full paths
        filesToBeOutput=filesToBeOutput, #list of full paths, using MODULEOUTDIR
        filesToBeDeleted=filesToBeDeleted, #list of full paths
        caching='PipelineDefault', #True, False, PipelineDefault
        cleanUpJobFilesOnFail='PipelineDefault', #True, False, PipelineDefault
        cleanUpPipelineJobsOnFail='PipelineDefault', #True, False, PipelineDefault
        executionEngine='PipelineDefault' #lsf, local, mockPrint, PipelineDefault
    )
    print ('dispense ' + jobName)

    outpath = filesToBeOutput_expanded[0]
    return outpath


def fix_header_pp(pipeline,module_subdir,adaptor_dir, module_dir, refdata_dir,  indiv_id, bam_arg):

    out_fn = os.path.basename(bam_arg[:-4])+'.reheadered.bam'
    out_module_path = os.path.join('$MODULEOUTDIR',out_fn)

    cmdStr = " ".join(['python', os.path.join(module_dir,'fix_header','fix_header.py'), bam_arg, out_module_path])



    moduleSubDir = module_subdir
    resources = {'maxmem':1, 'maxtime':3600};
    jobName = 'fix_header'
    inputFiles = [bam_arg]
    filesToBeOutput = [out_module_path]
    filesToBeDeleted = []

    filesToBeOutput_expanded = pipeline.dispense(
        moduleSubDir = moduleSubDir, #unique subdirectory under PIPELINEOUTDIR
        cmdStr = cmdStr, #command-line to run
        resources = resources, #dict with maxmem (in GB) and maxtime (in secs)
        jobName= jobName, #non-unique job-type name
        inputFiles=inputFiles, #list of full paths
        filesToBeOutput=filesToBeOutput, #list of full paths, using MODULEOUTDIR
        filesToBeDeleted=filesToBeDeleted, #list of full paths
        caching='PipelineDefault', #True, False, PipelineDefault
        cleanUpJobFilesOnFail='PipelineDefault', #True, False, PipelineDefault
        cleanUpPipelineJobsOnFail='PipelineDefault', #True, False, PipelineDefault
        executionEngine='PipelineDefault' #lsf, local, mockPrint, PipelineDefault
    )
    print ('dispense ' + jobName)

    outpath = filesToBeOutput_expanded[0]
    return outpath



def collapseHetSitesIntervals(pipeline, module_dir, module_subdir, indiv_id, sample_call_stats_filename, single_sample_vcf_filename, use_pcawg_contigs):
    if use_pcawg_contigs:
        reference_file = 'human_g1k_v37_decoy.fasta'
    else:
        reference_file = 'Homo_sapiens_assembly19.fasta'

    force_call_het_sites_intervals_filename = '%s.het_sites_2_forcecall.intervals' % indiv_id
    #calling the module using the firehose names as the parameters:
    module_libdir = os.path.join(module_dir, 'CollapseHetSitesIntervals')
    cmdStr = ' '.join(['python', os.path.join(module_libdir, 'collapse_het_sites_intervals.py')])

    referenceFilename = os.path.join(refdata_dir, 'public', reference_file)
    outputFilename = '$MODULEOUTDIR/' + force_call_het_sites_intervals_filename
    inputFiles = [sample_call_stats_filename, single_sample_vcf_filename]

    cmdStr += ' '.join(['',
    			'--gatk-jar', os.path.join(*[module_dir, 'haplotypecaller_sg', 'GenomeAnalysisTK-3.3-0-g37228af', 'GenomeAnalysisTK.jar']),
    			'--reference', referenceFilename,
			    '--sample-call-stats-filename', sample_call_stats_filename,
                '--sample-vcf-filename', single_sample_vcf_filename,
                '--disable_auto_index_creation_and_locking_when_reading_rods',
                '--force-call-het-intervals-filename', force_call_het_sites_intervals_filename,
			    '--module-dir', '$MODULEOUTDIR/'])
    moduleSubDir = module_subdir

    resources = {'maxmem': 10, 'maxtime': 14400}; #memory and time limits
    jobName = 'CollapseHetSitesIntervals' #job name
    moduleSubDir = module_subdir
    filesToBeOutput = [outputFilename]
    filesToBeDeleted = []

    filesToBeOutput_expanded = pipeline.dispense(
        moduleSubDir = moduleSubDir, #unique subdirectory under PIPELINEOUTDIR
        cmdStr = cmdStr, #command-line to run
        resources = resources, #dict with maxmem (in GB) and maxtime (in secs)
        jobName= jobName, #non-unique job-type name
        inputFiles=inputFiles, #list of full paths
        filesToBeOutput=filesToBeOutput, #list of full paths, using MODULEOUTDIR
        filesToBeDeleted=filesToBeDeleted, #list of full paths
        caching='PipelineDefault', #True, False, PipelineDefault
        cleanUpJobFilesOnFail='PipelineDefault', #True, False, PipelineDefault
        cleanUpPipelineJobsOnFail='PipelineDefault', #True, False, PipelineDefault
        executionEngine='PipelineDefault' #lsf, local, mockPrint, PipelineDefault
    )
    print ('dispense ' + jobName)

    return filesToBeOutput_expanded[0]

def genotypeGVCF(pipeline, module_dir, module_subdir, indiv_id, single_sample_gvcf_filename, use_pcawg_contigs):
    single_sample_vcf_filename = '%s.single_sample.vcf' % indiv_id
    #calling the module using the firehose names as the parameters:
    module_libdir = os.path.join(module_dir, 'haplotypecaller_sg')

    if use_pcawg_contigs:
        reference_file = 'human_g1k_v37_decoy.fasta'
    else:
        reference_file = 'Homo_sapiens_assembly19.fasta'

    cmdStr = ' '.join(['java', '-Xmx2g', '-jar', os.path.join(*[module_libdir, 'GenomeAnalysisTK-3.3-0-g37228af', 'GenomeAnalysisTK.jar'])])

    outputFilename = '$MODULEOUTDIR/' + single_sample_vcf_filename
    referenceFilename = os.path.join(refdata_dir, 'public', reference_file)
    inputFiles = [single_sample_gvcf_filename]

    cmdStr += ' '.join(['',
    			        '-T', 'GenotypeGVCFs',
                        '--variant', single_sample_gvcf_filename,
                        '--disable_auto_index_creation_and_locking_when_reading_rods',
                        '--out', outputFilename,
                        '--reference_sequence', referenceFilename])
    moduleSubDir = module_subdir

    resources = {'maxmem': 12, 'maxtime': 20000}; #memory and time limits
    jobName = 'GenotypeGVCFs' #job name
    moduleSubDir = module_subdir
    filesToBeOutput = [outputFilename]
    filesToBeDeleted = []

    filesToBeOutput_expanded = pipeline.dispense(
        moduleSubDir = moduleSubDir, #unique subdirectory under PIPELINEOUTDIR
        cmdStr = cmdStr, #command-line to run
        resources = resources, #dict with maxmem (in GB) and maxtime (in secs)
        jobName= jobName, #non-unique job-type name
        inputFiles=inputFiles, #list of full paths
        filesToBeOutput=filesToBeOutput, #list of full paths, using MODULEOUTDIR
        filesToBeDeleted=filesToBeDeleted, #list of full paths
        caching='PipelineDefault', #True, False, PipelineDefault
        cleanUpJobFilesOnFail='PipelineDefault', #True, False, PipelineDefault
        cleanUpPipelineJobsOnFail='PipelineDefault', #True, False, PipelineDefault
        executionEngine='PipelineDefault' #lsf, local, mockPrint, PipelineDefault
    )
    print ('dispense ' + jobName)

    return filesToBeOutput_expanded[0]


def concatenate_coverage_files_pm(pipeline, module_subdir, indiv_id, indiv_id_suffix, coverage_files):
    concatenated_coverage_filename = '%s.%s.uncorrected_target_order.coverage' % (indiv_id, indiv_id_suffix)
    cmdStr = 'cat ' + ' '.join(coverage_files) + '| sort | uniq > %s' % concatenated_coverage_filename

    moduleSubDir = module_subdir
    resources = {'maxmem':1, 'maxtime':300};
    jobName = 'concatenate_coverage' + "_" + indiv_id_suffix
    inputFiles = coverage_files
    filesToBeOutput = ['$MODULEOUTDIR/' + concatenated_coverage_filename]
    filesToBeDeleted = []

    filesToBeOutput_expanded = pipeline.dispense(
        moduleSubDir = moduleSubDir, #unique subdirectory under PIPELINEOUTDIR
        cmdStr = cmdStr, #command-line to run
        resources = resources, #dict with maxmem (in GB) and maxtime (in secs)
        jobName= jobName, #non-unique job-type name
        inputFiles=inputFiles, #list of full paths
        filesToBeOutput=filesToBeOutput, #list of full paths, using MODULEOUTDIR
        filesToBeDeleted=filesToBeDeleted, #list of full paths
        caching='PipelineDefault', #True, False, PipelineDefault
        cleanUpJobFilesOnFail='PipelineDefault', #True, False, PipelineDefault
        cleanUpPipelineJobsOnFail='PipelineDefault', #True, False, PipelineDefault
        executionEngine='PipelineDefault' #lsf, local, mockPrint, PipelineDefault
    )
    print ('dispense ' + jobName)

    outpath = filesToBeOutput_expanded[0]
    return outpath

def variant_bam_pm(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, bam, call_stats,M2_pass_vcf):

    #TBD add call_stats to the cmdStr and inputFiles
    ################################ CUSTOM FOR VARIANT BAM ##################
    rule = 'data/rules2.json'
    bam_filename = os.path.basename(bam)
    bam_basename = bam_filename[:-4]
    bam_outname = bam_basename + '.var.bam'
    module_libdir = os.path.join(module_dir,'VariantBam') ## replace VariantBam with your task folder name
    rule_file = os.path.join(module_libdir, rule) 

    ld_path = 'export LD_LIBRARY_PATH=' + module_libdir + '/bamtools-2.1.0/lib; echo $LD_LIBRARY_PATH; '
    cmdStr = ld_path
    cmdStr += os.path.join(module_dir,'VariantBam','variant')
    #cmdStr += ' '.join([' -i',bam,'-f',rule_file,'-o',bam_outname,'-k',os.path.join(module_libdir,'')]
    #cmdStr += ' '.join([' -i',bam,'-f',rule_file,'-o',bam_outname])

    if call_stats == None:
        cmdStr += ' '.join([' ', bam,'-r',rule_file, '-o', bam_outname])
    elif M2_pass_vcf == None:
        cmdStr += ' '.join([' ', bam,'-r',rule_file, '-c', call_stats,' -o', bam_outname])
    else:
        cmdStr += ' '.join([' ',bam,'-r',rule_file, '-c', call_stats,'-g',M2_pass_vcf,' -o', bam_outname])
    
    module_variant_bam_output = '$MODULEOUTDIR/' + bam_basename + '.var.bam'
    #module_variant_bam_output_bai = module_variant_bam_output + '.bai'
    module_qc_stats = '$MODULEOUTDIR/qcreport.txt'
    module_merged_callstats = '$MODULEOUTDIR/merged_rules.bed'
    ######################################################
    
    moduleSubDir = module_subdir

    ####################### set the resource requirements here
    resources = {'maxmem':1, 'maxtime':2000 * 3600};
    jobName = 'variant_bam'
    ###################### INPUT FILES
    inputFiles = [bam, call_stats,M2_pass_vcf]
    ##### LIST OF FILES TO BE OUTPUT
    filesToBeOutput = [module_variant_bam_output, module_qc_stats, module_merged_callstats]
    filesToBeDeleted = []
    
    filesToBeOutput_expanded = pipeline.dispense(
        moduleSubDir = moduleSubDir, #unique subdirectory under PIPELINEOUTDIR
        cmdStr = cmdStr, #command-line to run
        resources = resources, #dict with maxmem (in GB) and maxtime (in secs)
        jobName= jobName, #non-unique job-type name
        inputFiles=inputFiles, #list of full paths
        filesToBeOutput=filesToBeOutput, #list of full paths, using MODULEOUTDIR
        filesToBeDeleted=filesToBeDeleted, #list of full paths
        caching='PipelineDefault', #True, False, PipelineDefault
        cleanUpJobFilesOnFail='PipelineDefault', #True, False, PipelineDefault
        cleanUpPipelineJobsOnFail='PipelineDefault', #True, False, PipelineDefault
        executionEngine='PipelineDefault' #lsf, local, mockPrint, PipelineDefault
    )
    
    variant_bam_output = filesToBeOutput_expanded[0]
    qc_stats = filesToBeOutput_expanded[1]
    merged_callstats = filesToBeOutput_expanded[2]
    return (variant_bam_output, qc_stats, merged_callstats)


def svaba_pm(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, bam_tumor, bam_normal, indiv_id, cores): 
    svaba_module_dir = "/opt/svaba_install"
    cmdStr = os.path.join(svaba_module_dir,'svaba','bin','svaba')
    
    svaba_exclusions = os.path.join(refdata_dir,'public','svaba_exclusions.bed')
    dbsnp_indel = os.path.join(refdata_dir,'public','dbsnp_indel.vcf')   
    ref_genome  = os.path.join(refdata_dir,'public','Homo_sapiens_assembly19.fasta')
    
    cmdStr += ' '.join([' run -t',bam_tumor,'-n',bam_normal,'-p',cores, '-a', indiv_id, '-B', svaba_exclusions, '-D', dbsnp_indel, '-G', ref_genome ])



    module_svaba_somatic_vcf = '$MODULEOUTDIR/' + indiv_id + '.svaba.somatic.sv.vcf'
    module_svaba_germline_vcf = '$MODULEOUTDIR/' + indiv_id + '.svaba.germline.sv.vcf'
    module_svaba_somatic_indel_vcf = '$MODULEOUTDIR/' + indiv_id + '.svaba.somatic.indel.vcf'
    module_svaba_germline_indel_vcf = '$MODULEOUTDIR/' + indiv_id + '.svaba.germline.indel.vcf'

    module_svaba_contigs_bam = '$MODULEOUTDIR/' + indiv_id + '.contigs.bam'


    # module_svaba_somatic_vcf = '$MODULEOUTDIR/' + indiv_id + '.broad-svaba.DATECODE.somatic.sv.vcf'
    # module_svaba_germline_vcf = '$MODULEOUTDIR/' + indiv_id + '.broad-svaba.DATECODE.germline.sv.vcf'
    # module_svaba_somatic_indel_vcf = '$MODULEOUTDIR/' + indiv_id + '.broad-svaba.DATECODE.somatic.indel.vcf'
    # module_svaba_germline_indel_vcf = '$MODULEOUTDIR/' + indiv_id + '.broad-svaba.DATECODE.germline.indel.vcf'


    # module_svaba_r2c_bam = '$MODULEOUTDIR/' + indiv_id + '.r2c.bam'
    # module_svaba_contigs_bam = '$MODULEOUTDIR/' + indiv_id + '.contigs.bam'
    # module_svaba_contigs_bai = '$MODULEOUTDIR/' + indiv_id + '.contigs.bam.bai'
    # module_svaba_bps = '$MODULEOUTDIR/' + indiv_id + '.bps.txt.gz'
    # module_svaba_cigar = '$MODULEOUTDIR/' + indiv_id + '.cigarmap.txt.gz'
    # module_svaba_aln = '$MODULEOUTDIR/' + indiv_id + '.alignments.txt.gz'
    # module_svaba_dsc = '$MODULEOUTDIR/' + indiv_id + '.discordant.txt.gz'
    # module_svaba_contigs_all = '$MODULEOUTDIR/' + indiv_id + '.contigs_all.sam.gz'
    # module_svaba_contigs_del = '$MODULEOUTDIR/' + indiv_id + '.contigs.sam'
    ######################################################
    
    moduleSubDir = module_subdir

    ####################### set the resource requirements here
    #observed 6gb on bam4
    resources = {'maxmem':15, 'maxtime':14400};
    jobName = 'svaba'
    ###################### INPUT FILES
    inputFiles = [bam_tumor, bam_normal]
    ##### LIST OF FILES TO BE OUTPUT
    filesToBeOutput = [
                        module_svaba_somatic_vcf, module_svaba_germline_vcf, module_svaba_somatic_indel_vcf, module_svaba_germline_indel_vcf,
                        module_svaba_contigs_bam
                       ]
    filesToBeDeleted = [] #[module_svaba_contigs_del]
    
    filesToBeOutput_expanded = pipeline.dispense(
        moduleSubDir = moduleSubDir, #unique subdirectory under PIPELINEOUTDIR
        cmdStr = cmdStr, #command-line to run
        resources = resources, #dict with maxmem (in GB) and maxtime (in secs)
        jobName= jobName, #non-unique job-type name
        inputFiles=inputFiles, #list of full paths
        filesToBeOutput=filesToBeOutput, #list of full paths, using MODULEOUTDIR
        filesToBeDeleted=filesToBeDeleted, #list of full paths
        caching='PipelineDefault', #True, False, PipelineDefault
        cleanUpJobFilesOnFail='PipelineDefault', #True, False, PipelineDefault
        cleanUpPipelineJobsOnFail='PipelineDefault', #True, False, PipelineDefault
        executionEngine='PipelineDefault' #lsf, local, mockPrint, PipelineDefault
    )


    module_svaba_somatic_vcf_outfile = filesToBeOutput_expanded[0]
    module_svaba_germline_vcf_outfile = filesToBeOutput_expanded[1]
    module_svaba_somatic_indel_vcf_outfile = filesToBeOutput_expanded[2]
    module_svaba_germline_indel_vcf_outfile = filesToBeOutput_expanded[3]
    broad_files = filesToBeOutput_expanded[4:]
    return (module_svaba_somatic_vcf_outfile,module_svaba_germline_vcf_outfile,module_svaba_somatic_indel_vcf_outfile,module_svaba_germline_indel_vcf_outfile,
        broad_files)




# def svaba_pm(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, bam_tumor, bam_normal, indiv_id, cores): 

#     ld_path_list = [os.path.join(module_dir, 'Snowman', 'hts'),
#                     os.path.join(module_dir, 'Snowman', 'bamtools-2.1.0','lib'),
#                     ]
#     ld_path = 'export LD_LIBRARY_PATH=' + ':'.join(ld_path_list) + '; '
#     cmdStr = ld_path
#     cmdStr += os.path.join(module_dir,'Snowman','svaba run')
    
#     ponStr = os.path.join(module_dir,'Snowman','lung_snow24_pon.txt.gz')
#     mskStr = os.path.join(module_dir,'Snowman','um75-hs37d5.bed.gz')

#     refr  = refdata_dir + '/public/Homo_sapiens_assembly19.fasta'
#     cmdStr += ' '.join([' -t',bam_tumor,'-n',bam_normal,'-p',cores, '-g', refr, '-e', '0.05', '-a', indiv_id, '-c 500000 --no-zip', '-q', ponStr, '-m', mskStr])

#     module_svaba_somatic_vcf = '$MODULEOUTDIR/' + indiv_id + '.broad-svaba.DATECODE.somatic.sv.vcf'
#     module_svaba_germline_vcf = '$MODULEOUTDIR/' + indiv_id + '.broad-svaba.DATECODE.germline.sv.vcf'
#     module_svaba_somatic_indel_vcf = '$MODULEOUTDIR/' + indiv_id + '.broad-svaba.DATECODE.somatic.indel.vcf'
#     module_svaba_germline_indel_vcf = '$MODULEOUTDIR/' + indiv_id + '.broad-svaba.DATECODE.germline.indel.vcf'


#     module_svaba_r2c_bam = '$MODULEOUTDIR/' + indiv_id + '.r2c.bam'
#     module_svaba_contigs_bam = '$MODULEOUTDIR/' + indiv_id + '.contigs.bam'
#     module_svaba_contigs_bai = '$MODULEOUTDIR/' + indiv_id + '.contigs.bam.bai'
#     module_svaba_bps = '$MODULEOUTDIR/' + indiv_id + '.bps.txt.gz'
#     module_svaba_cigar = '$MODULEOUTDIR/' + indiv_id + '.cigarmap.txt.gz'
#     module_svaba_aln = '$MODULEOUTDIR/' + indiv_id + '.alignments.txt.gz'
#     module_svaba_dsc = '$MODULEOUTDIR/' + indiv_id + '.discordant.txt.gz'
#     module_svaba_contigs_all = '$MODULEOUTDIR/' + indiv_id + '.contigs_all.sam.gz'
#     module_svaba_contigs_del = '$MODULEOUTDIR/' + indiv_id + '.contigs.sam'


#     ######################################################
    
#     moduleSubDir = module_subdir

#     ####################### set the resource requirements here
#     #observed 6gb on bam4
#     resources = {'maxmem':15, 'maxtime':14400};
#     jobName = 'svaba'
#     ###################### INPUT FILES
#     inputFiles = [bam_tumor, bam_normal]
#     ##### LIST OF FILES TO BE OUTPUT
#     filesToBeOutput = [
#                         module_svaba_somatic_vcf, module_svaba_germline_vcf, module_svaba_somatic_indel_vcf, module_svaba_germline_indel_vcf,
#                         module_svaba_r2c_bam,
#                         module_svaba_contigs_bam, module_svaba_contigs_bai,
#                         module_svaba_bps, module_svaba_aln, module_svaba_dsc, module_svaba_contigs_all, module_svaba_cigar ## JEREMIAH 150410
#                        ]
#     filesToBeDeleted = [] #[module_svaba_contigs_del]
    
#     filesToBeOutput_expanded = pipeline.dispense(
#         moduleSubDir = moduleSubDir, #unique subdirectory under PIPELINEOUTDIR
#         cmdStr = cmdStr, #command-line to run
#         resources = resources, #dict with maxmem (in GB) and maxtime (in secs)
#         jobName= jobName, #non-unique job-type name
#         inputFiles=inputFiles, #list of full paths
#         filesToBeOutput=filesToBeOutput, #list of full paths, using MODULEOUTDIR
#         filesToBeDeleted=filesToBeDeleted, #list of full paths
#         caching='PipelineDefault', #True, False, PipelineDefault
#         cleanUpJobFilesOnFail='PipelineDefault', #True, False, PipelineDefault
#         cleanUpPipelineJobsOnFail='PipelineDefault', #True, False, PipelineDefault
#         executionEngine='PipelineDefault' #lsf, local, mockPrint, PipelineDefault
#     )


#     module_svaba_somatic_vcf_outfile = filesToBeOutput_expanded[0]
#     module_svaba_germline_vcf_outfile = filesToBeOutput_expanded[1]
#     module_svaba_somatic_indel_vcf_outfile = filesToBeOutput_expanded[2]
#     module_svaba_germline_indel_vcf_outfile = filesToBeOutput_expanded[3]
#     broad_files = filesToBeOutput_expanded[4:]
#     return (module_svaba_somatic_vcf_outfile,module_svaba_germline_vcf_outfile,module_svaba_somatic_indel_vcf_outfile,module_svaba_germline_indel_vcf_outfile,broad_files)

def sv_merge_pm(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, svaba_vcf, dranger_vcf, indiv_id):

    ld_path_list = module_dir + '/Snowman/hts:' + module_dir + '/Snowman/bamtools-2.1.0/lib'
    ld_path = 'LD_LIBRARY_PATH=' + ld_path_list + '; export LD_LIBRARY_PATH; '
    cmdStr = ld_path
    #note - sv_merge is part of the Snowman binary.
    cmdStr += os.path.join(module_dir,'Snowman','svaba vcf')

    cmdStr += ' '.join([' -r',dranger_vcf,'-s',svaba_vcf,'-a', indiv_id, '--no-zip'])

    module_merged_vcf = '$MODULEOUTDIR/' + indiv_id + '.broad-merged.DATECODE.somatic.sv.vcf'
    module_merged_tbi = '$MODULEOUTDIR/' + indiv_id + '.broad-merged.DATECODE.somatic.sv.vcf'
    
    ######################################################
    
    moduleSubDir = module_subdir
    
    ####################### set the resource requirements here
    resources = {'maxmem':4, 'maxtime':14400};
    jobName = 'sv_merge'
    ###################### INPUT FILES
    inputFiles = [svaba_vcf, dranger_vcf]
    ##### LIST OF FILES TO BE OUTPUT
    filesToBeOutput = [module_merged_vcf]
    filesToBeDeleted = []
    
    filesToBeOutput_expanded = pipeline.dispense(
        moduleSubDir = moduleSubDir, #unique subdirectory under PIPELINEOUTDIR
        cmdStr = cmdStr, #command-line to run
        resources = resources, #dict with maxmem (in GB) and maxtime (in secs)
        jobName= jobName, #non-unique job-type name
        inputFiles=inputFiles, #list of full paths
        filesToBeOutput=filesToBeOutput, #list of full paths, using MODULEOUTDIR
        filesToBeDeleted=filesToBeDeleted, #list of full paths
        caching='PipelineDefault', #True, False, PipelineDefault
        cleanUpJobFilesOnFail='PipelineDefault', #True, False, PipelineDefault
        cleanUpPipelineJobsOnFail='PipelineDefault', #True, False, PipelineDefault
        executionEngine='PipelineDefault' #lsf, local, mockPrint, PipelineDefault
        )
    
    output_merged_vcf =filesToBeOutput_expanded[0]

    return output_merged_vcf


def oncotator(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, indiv_id, maflite_file):
    
    outfile = indiv_id + '.vcf'

    #building command line string:
    #cmdStr = os.path.join(module_dir,'oncotator','oncotator')
    cmdStr = 'oncotator'
    cmdStr += ' '.join([' --db-dir', os.path.join(module_dir, 'empty/'), '-o VCF', maflite_file, outfile, 'hg19'])

    mutect_vcf = '$MODULEOUTDIR/' + indiv_id + '.vcf'  
    moduleSubDir = module_subdir

    ####################### set the resource requirements here
    resources = {'maxmem':1, 'maxtime':14400};
    jobName = 'oncotator'
    ###################### INPUT FILES
    inputFiles = [maflite_file]
    ##### LIST OF FILES TO BE OUTPUT
    filesToBeOutput = [mutect_vcf]
    filesToBeDeleted = []
    
    filesToBeOutput_expanded = pipeline.dispense(
        moduleSubDir = moduleSubDir, #unique subdirectory under PIPELINEOUTDIR
        cmdStr = cmdStr, #command-line to run
        resources = resources, #dict with maxmem (in GB) and maxtime (in secs)
        jobName= jobName, #non-unique job-type name
        inputFiles=inputFiles, #list of full paths
        filesToBeOutput=filesToBeOutput, #list of full paths, using MODULEOUTDIR
        filesToBeDeleted=filesToBeDeleted, #list of full paths
        caching='PipelineDefault', #True, False, PipelineDefault
        cleanUpJobFilesOnFail='PipelineDefault', #True, False, PipelineDefault
        cleanUpPipelineJobsOnFail='PipelineDefault', #True, False, PipelineDefault
        executionEngine='PipelineDefault' #lsf, local, mockPrint, PipelineDefault
    )
    
    mutect_vcf_output = filesToBeOutput_expanded[0]
    return (mutect_vcf_output)

def callstats_to_maflite(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, indiv_id, call_stats): 

 #calling the module using the firehose names as the parameters:
    module_libdir = os.path.join(module_dir,'callstats_to_maflite')
    cmdStr = os.path.join(adaptor_dir,'run_module.py') #don't change
    cmdStr += ' --module_libdir ' + module_libdir + ' ' #don't change

    #Add in the parameters as the module is called from firehose:
    cmdStr += ' '.join(['--input.call_stats.file', call_stats, 
                        '--genome.build 37',
                        '--mode FSTAR', 
                        '--output.prefix',indiv_id, 
                        '--triallelic_mode_KEEP_or_REJECT REJECT',
                        '--extra.columns tumor_f,init_t_lod,t_lod_fstar,t_alt_count,t_ref_count,judgement,n_alt_count,n_ref_count,observed_in_normals_count',
                        '--f_threshold','0'])

    #change these variables to the firehose parameters. 
    maflite = '$MODULEOUTDIR/' + indiv_id + '.maf' # maflite output file.
        
    moduleSubDir = module_subdir
        
    resources = {'maxmem':1, 'maxtime':14400}; #memory and time limits
    jobName = 'callstats_to_maflite' #job name
    inputFiles = [call_stats] #parameters from the function to be passed in. 
    #change to include the variables you listed above to save from your function:
    filesToBeOutput = [maflite]
    
    #STOP CHANGING HERE
    #########################
    filesToBeDeleted = []
    
    filesToBeOutput_expanded = pipeline.dispense(
        moduleSubDir = moduleSubDir, #unique subdirectory under PIPELINEOUTDIR
        cmdStr = cmdStr, #command-line to run
        resources = resources, #dict with maxmem (in GB) and maxtime (in secs)
        jobName= jobName, #non-unique job-type name
        inputFiles=inputFiles, #list of full paths
        filesToBeOutput=filesToBeOutput, #list of full paths, using MODULEOUTDIR
        filesToBeDeleted=filesToBeDeleted, #list of full paths
        caching='PipelineDefault', #True, False, PipelineDefault
        cleanUpJobFilesOnFail='PipelineDefault', #True, False, PipelineDefault
        cleanUpPipelineJobsOnFail='PipelineDefault', #True, False, PipelineDefault
        executionEngine='PipelineDefault' #lsf, local, mockPrint, PipelineDefault
    )
    
    maflite_file = filesToBeOutput_expanded[0]
    return maflite_file


def ReCapSeg_Coverage_pp(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir,  indiv_id, bam, chrom, indiv_id_suffix, use_pcawg_contigs):

    if use_pcawg_contigs:
        reference_file = 'human_g1k_v37_decoy.fasta'
    else:
        reference_file = 'Homo_sapiens_assembly19.fasta'

    #calling the module using the firehose names as the parameters:
    module_libdir = os.path.join(module_dir, 'ReCapSegCoverageWGS')

    # commandLine=sh <libdir>wrapper.sh <libdir> "--jar " <libdir>/coverage.jar <sample.ID> <sample.BAM> <bed.track> <reference.sequence>

    #############################################################
    #START CHANGING HERE:
    #Add in the parameters as the module is called from firehose:
    #note memory is set both in cmd line and in resources
    cmdStr = ' '.join(['java', '-Xmx10g', '-jar', os.path.join(module_libdir, 'coverage.jar')]) 
    cmdStr += ' '.join(['', '--analysis_type', 'BaitDepthWalker'])

    sample_id = indiv_id + '_' + indiv_id_suffix
    outputFilename = '$MODULEOUTDIR/' + sample_id + '.coverage'
    intervalsFilename = os.path.join(refdata_dir, 'public', 'wgs_1kbp_intervals.bed')
    referenceFilename = os.path.join(refdata_dir, 'public', reference_file)

    cmdStr += ' '.join(['', 
                        '--input_file', bam, 
                        '--intervals', chrom,  # process one chromosome at a time
                        '--out', outputFilename, 
                        '--bed', intervalsFilename, 
                        '--reference_sequence', referenceFilename])
    moduleSubDir = module_subdir

    resources = {'maxmem': 12, 'maxtime': 14400}; #memory and time limits
    jobName = 'ReCapSeg_Coverage' + '_' + indiv_id_suffix #job name

    inputFiles = [bam] #parameters from the function to be passed in.
    #changed based on fh previous executions
    #change to include the variables you listed above to save from your function:
    filesToBeOutput = [outputFilename]

    #STOP CHANGING HERE
    #########################
    filesToBeDeleted = []
    filesToBeOutput_expanded = pipeline.dispense(
        moduleSubDir=moduleSubDir, #unique subdirectory under PIPELINEOUTDIR
        cmdStr=cmdStr, #command-line to run
        resources=resources, #dict with maxmem (in GB) and maxtime (in secs)
        jobName=jobName, #non-unique job-type name
        inputFiles=inputFiles, #list of full paths
        filesToBeOutput=filesToBeOutput, #list of full paths, using MODULEOUTDIR
        filesToBeDeleted=filesToBeDeleted, #list of full paths
        caching='PipelineDefault', #True, False, PipelineDefault
        cleanUpJobFilesOnFail='PipelineDefault', #True, False, PipelineDefault
        cleanUpPipelineJobsOnFail='PipelineDefault', #True, False, PipelineDefault
        executionEngine='PipelineDefault' #lsf, local, mockPrint, PipelineDefault
    )
    print ('dispense ' + jobName)

    re_capseg_cov = filesToBeOutput_expanded[0]
    return re_capseg_cov





def contest_pp(pipeline, module_subdir,adaptor_dir, module_dir, refdata_dir,  indiv_id, bam_tumor, bam_normal, tmp_dir, use_pcawg_contigs):
    #TBD remove hardcoded tmpdir from manifest
    
    #calling the module using the firehose names as the parameters:
    module_libdir = os.path.join(module_dir,'contest')
    cmdStr = os.path.join(adaptor_dir,'run_module.py') #don't change
    cmdStr += ' --module_libdir ' + module_libdir + ' ' #don't change

    if use_pcawg_contigs:
        reference_file = 'human_g1k_v37_decoy.fasta'
    else:
        reference_file = 'Homo_sapiens_assembly19.fasta'

    #Add in the parameters as the module is called from firehose:
    cmdStr += ' '.join(['--reference.file', os.path.join(refdata_dir,'public',reference_file),
                        '--intervals', os.path.join(refdata_dir,'public','gaf_20111020+broad_wex_1.1_hg19.bed'), 
                        '--sample.id', os.path.basename(bam_tumor), 
                        '--clean.bam.file', bam_tumor, 
                        '--normal.bam',bam_normal, 
                        '--genotypes.vcf none', 
                        '--pop.vcf', os.path.join(refdata_dir,'public','hg19_population_stratified_af_hapmap_3.3.fixed.vcf'), 
                        '--force.array.free true',
                        '--snp.six.bed', os.path.join(refdata_dir,'public','SNP6.hg19.interval_list'),
                        '--job.spec.memory','2', '--tmp.dir', tmp_dir #TBD job spec memory was 2; manifest also changed.
                        ])

    tumor_base_name = os.path.basename(bam_tumor)
    
    #change these variables to the firehose parameters. 
    contest_value = '$MODULEOUTDIR/' + tumor_base_name + '.contamination.txt.firehose' #this is a file, not a single value. 
        
    moduleSubDir = module_subdir
        
    resources = {'maxmem':4, 'maxtime':14400}; #memory and time limits
    jobName = 'contest' #job name
    inputFiles = [bam_tumor, bam_normal] #parameters from the function to be passed in. 
    #change to include the variables you listed above to save from your function:
    filesToBeOutput = [contest_value]
    
    #STOP CHANGING HERE
    #########################
    filesToBeDeleted = []
    
    filesToBeOutput_expanded = pipeline.dispense(
        moduleSubDir = moduleSubDir, #unique subdirectory under PIPELINEOUTDIR
        cmdStr = cmdStr, #command-line to run
        resources = resources, #dict with maxmem (in GB) and maxtime (in secs)
        jobName= jobName, #non-unique job-type name
        inputFiles=inputFiles, #list of full paths
        filesToBeOutput=filesToBeOutput, #list of full paths, using MODULEOUTDIR
        filesToBeDeleted=filesToBeDeleted, #list of full paths
        caching='PipelineDefault', #True, False, PipelineDefault
        cleanUpJobFilesOnFail='PipelineDefault', #True, False, PipelineDefault
        cleanUpPipelineJobsOnFail='PipelineDefault', #True, False, PipelineDefault
        executionEngine='PipelineDefault' #lsf, local, mockPrint, PipelineDefault
    )
    
    contest_value_file = filesToBeOutput_expanded[0]
    return contest_value_file


def extract_bam_id_pp(pipeline, module_subdir, indiv_id, adaptor_dir, module_dir,refdata_dir, bam_normal):

    module_libdir = os.path.join(module_dir, 'extract_bam_id')
    cmdStr = os.path.join(adaptor_dir,'run_module.py') #don't change
    cmdStr += ' --module_libdir ' + module_libdir + ' ' #don't change

    normal_base_name = os.path.basename(bam_normal)
    normal_base_name = normal_base_name[:-4]

    #Add in the parameters as the module is called from firehose:
    cmdStr += ' '.join([#'--sample_id', normal_base_name, 
                        #'--annotation_nam bamfile_id ',
                        '--bam.file', bam_normal,
                        ])

    #change these variables to the firehose parameters. 
    bam_id_file = '$MODULEOUTDIR/upload.txt' #this is a file, not a single value. How pass first line to Haplotypecaller
    moduleSubDir = module_subdir
        
    resources = {'maxmem':1, 'maxtime':14400}; #memory and time limits
    jobName = 'contest' #job name
    inputFiles = [bam_normal] #parameters from the function to be passed in. 
    #change to include the variables you listed above to save from your function:
    filesToBeOutput = [bam_id_file]
    
    #STOP CHANGING HERE
    #########################
    filesToBeDeleted = []
    
    filesToBeOutput_expanded = pipeline.dispense(
        moduleSubDir = moduleSubDir, #unique subdirectory under PIPELINEOUTDIR
        cmdStr = cmdStr, #command-line to run
        resources = resources, #dict with maxmem (in GB) and maxtime (in secs)
        jobName= jobName, #non-unique job-type name
        inputFiles=inputFiles, #list of full paths
        filesToBeOutput=filesToBeOutput, #list of full paths, using MODULEOUTDIR
        filesToBeDeleted=filesToBeDeleted, #list of full paths
        caching='PipelineDefault', #True, False, PipelineDefault
        cleanUpJobFilesOnFail='PipelineDefault', #True, False, PipelineDefault
        cleanUpPipelineJobsOnFail='PipelineDefault', #True, False, PipelineDefault
        executionEngine='PipelineDefault' #lsf, local, mockPrint, PipelineDefault
    )
    
    bam_id_file = filesToBeOutput_expanded[0]
    return bam_id_file


def haplotype_caller_sg_pp(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, indiv_id, bam_normal): #get the sample id from extract_bam_id_pp
    
    
    module_libdir = os.path.join(module_dir,'haplotypecaller_sg')

    normal_base_name = os.path.basename(bam_normal)
    normal_base_name = normal_base_name[:-4]
    dbsnp_vcf = os.path.join(refdata_dir,'public','dbsnp_134_b37.leftAligned.vcf')
    #TBD - bam_id_file will not exist at the time this function is called, it must be read at job execution time from the shell.
    #bamfile_id = open(bam_id_file,'r').read().split('\n')[0]
    reference_sequence = os.path.join(refdata_dir, 'public', 'Homo_sapiens_assembly19.fasta')
    #reference_sequence = os.path.join(refdata_dir,'public','human_g1k_v37_decoy.fasta')

    #NOTE - prepare.sh script sets JAVA_HOME explicitly, to use Java1.8 for itself.
    #Add in the parameters as the module is called from firehose:
    prepare_args = {'reference_sequence':reference_sequence,
                    'number_scatter_jobs':'50',
                    'jobspecmemory':'5',
                    'bamfile':bam_normal,
                    #'sampleID':'`cat ' + bam_id_file + ' | awk 1 ORS=" "`', #file contents with no linefeed
                    'sampleID':indiv_id,
                    'dbsnp_file':dbsnp_vcf,
                    }

    #change these variables to the firehose parameters. 
    gvcf = '$MODULEOUTDIR/' + indiv_id + '.gvcf.gz'
    gvcf_tbi = '$MODULEOUTDIR/' + indiv_id + '.gvcf.gz.tbi'

    moduleSubDir = module_subdir
    scatter_maxmem = 7
    gather_maxmem = 2
    sg_module_subdir = os.path.join(module_subdir,'sg')
    scatter_width = 50

    #resources = {'maxmem':20, 'maxtime':14400}; #memory and time limits
    #jobName = 'haplotypecaller_sg' #job name
    input_files = [bam_normal] #parameters from the function to be passed in. 
    #change to include the variables you listed above to save from your function:
    output_files = [gvcf, gvcf_tbi]

    haplotype_caller_output_files = scatter_gather_pp(pipeline, sg_module_subdir, module_libdir, adaptor_dir, input_files, output_files, prepare_args, scatter_width, scatter_maxmem, gather_maxmem)

    haplotype_caller_gvcf = haplotype_caller_output_files[0]
    haplotype_caller_tbi = haplotype_caller_output_files[1]
    return (haplotype_caller_gvcf, haplotype_caller_tbi)




def haplotype_caller_pp(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, bam_id_file, bam_normal): #get the sample id from extract_bam_id_pp
    #DEPRECATED - no longer supported
    module_libdir = os.path.join(module_dir,'haplotypecaller')
    #calling the module using the firehose names as the parameters:
    cmdStr = os.path.join(adaptor_dir,'run_module.py') #don't change
    cmdStr += ' --module_libdir ' + module_libdir + ' ' #don't change

    normal_base_name = os.path.basename(bam_normal)
    normal_base_name = normal_base_name[:-4]
    dbsnp_vcf = os.path.join(refdata_dir,'public','dbsnp_134_b37.leftAligned.vcf')

    #Add in the parameters as the module is called from firehose:
    cmdStr += ' '.join(['--reference', os.path.join(refdata_dir,'public','human_g1k_v37_decoy.fasta') , 
                        '--bam.file', bam_normal, 
                        '--dbsnp', dbsnp_vcf, 
                        '--sample', '`cat ' + bam_id_file + '`'],
                       )

    #change these variables to the firehose parameters. 
    gvcf = '$MODULEOUTDIR/' + normal_base_name + '.gvcf.gz'
    gvcf_tbi = '$MODULEOUTDIR/' + normal_base_name + '.gvcf.gz.tbi'

    moduleSubDir = module_subdir
        
    resources = {'maxmem':1, 'maxtime':14400}; #memory and time limits
    jobName = 'contest' #job name
    inputFiles = [bam_id_file,bam_normal] #parameters from the function to be passed in. 
    #change to include the variables you listed above to save from your function:
    filesToBeOutput = [gvcf, gvcf_tbi]
    
    #STOP CHANGING HERE
    #########################
    filesToBeDeleted = []
    
    filesToBeOutput_expanded = pipeline.dispense(
        moduleSubDir = moduleSubDir, #unique subdirectory under PIPELINEOUTDIR
        cmdStr = cmdStr, #command-line to run
        resources = resources, #dict with maxmem (in GB) and maxtime (in secs)
        jobName= jobName, #non-unique job-type name
        inputFiles=inputFiles, #list of full paths
        filesToBeOutput=filesToBeOutput, #list of full paths, using MODULEOUTDIR
        filesToBeDeleted=filesToBeDeleted, #list of full paths
        caching='PipelineDefault', #True, False, PipelineDefault
        cleanUpJobFilesOnFail='PipelineDefault', #True, False, PipelineDefault
        cleanUpPipelineJobsOnFail='PipelineDefault', #True, False, PipelineDefault
        executionEngine='PipelineDefault' #lsf, local, mockPrint, PipelineDefault
    )
    
    return filesToBeOutput_expanded



##########
#####JaBbA preprocessing
##########
def fragCounter_pm(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, input_bam):

    #calling the module using the firehose names as the parameters:
    module_libdir = os.path.join(module_dir,'fragCounter')
    cmdStr = os.path.join(adaptor_dir,'run_module.py') #don't change
    cmdStr += ' --module_libdir ' + module_libdir + ' ' #don't change


    #############################################################
    #START CHANGING HERE:
    #Add in the parameters as the module is called from firehose:
    memory = 10 # value needed both in cmdline and for Pipette resources; set it here once to avoid c/p errors
    cmdStr += ' '.join(['--Bam', input_bam ,
                        '--Cov', 'NULL',
                        '--WindowSize', '200',
                        '--job.spec.memory', str(memory)])        
    moduleSubDir = module_subdir
    
    resources = {'maxmem':memory+2, 'maxtime':20000}; #memory and time limits
    jobName = 'fragCounter' #job name
    inputFiles = [input_bam] #parameters from the function to be passed in.
    #changed based on fh previous executions
    #change to include the variables you listed above to save from your function:
    filesToBeOutput = ['$MODULEOUTDIR/cov.rds']
    
    #STOP CHANGING HERE
    #########################
    filesToBeDeleted = []
    
    filesToBeOutput_expanded = pipeline.dispense(
        moduleSubDir = moduleSubDir, #unique subdirectory under PIPELINEOUTDIR
        cmdStr = cmdStr, #command-line to run
        resources = resources, #dict with maxmem (in GB) and maxtime (in secs)
        jobName= jobName, #non-unique job-type name
        inputFiles=inputFiles, #list of full paths
        filesToBeOutput=filesToBeOutput, #list of full paths, using MODULEOUTDIR
        filesToBeDeleted=filesToBeDeleted, #list of full paths
        caching='PipelineDefault', #True, False, PipelineDefault
        cleanUpJobFilesOnFail='PipelineDefault', #True, False, PipelineDefault
        cleanUpPipelineJobsOnFail='PipelineDefault', #True, False, PipelineDefault
        executionEngine='PipelineDefault' #lsf, local, mockPrint, PipelineDefault
    )
    print ('dispense ' + jobName)

    JaBbA_cov_rds = filesToBeOutput_expanded[0]
    return JaBbA_cov_rds



##########
#####OxoQ
##########
def OxoQ_pm(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir,  indiv_id, bam_tumor, bam_normal, use_pcawg_contigs):
    if use_pcawg_contigs:
        reference_file = 'human_g1k_v37_decoy.fasta'
    else:
        reference_file = 'Homo_sapiens_assembly19.fasta'


    #calling the module using the firehose names as the parameters:
    module_libdir = os.path.join(module_dir,'OxoQ')
    cmdStr = 'python ' + os.path.join(adaptor_dir,'run_module.py') #don't change
    cmdStr += ' --module_libdir ' + module_libdir + ' ' #don't change

    #############################################################
    #START CHANGING HERE:
    #Add in the parameters as the module is called from firehose:
    cmdStr += ' '.join(['--id', indiv_id ,
                        '--bam_file', bam_tumor,
                        '--reference_fasta', os.path.join(refdata_dir,'public',reference_file),
                        '--dbSNP_vcf', os.path.join(refdata_dir,'public','dbsnp_134_b37.leftAligned.vcf'),
                        '--context', 'CCG',
                        '--output','.'])

    #STOP_AFTER = 100000000

    moduleSubDir = module_subdir


    resources = {'maxmem':3, 'maxtime':40000}; #memory and time limits
    jobName = 'OxoQ_A' #job name
    inputFiles = [bam_tumor] #parameters from the function to be passed in.
    #change to include the variables you listed above to save from your function:
    filesToBeOutput = ['$MODULEOUTDIR/' + indiv_id + '.oxoQ.txt']

    #STOP CHANGING HERE
    #########################
    filesToBeDeleted = []

    filesToBeOutput_expanded = pipeline.dispense(
        moduleSubDir = moduleSubDir, #unique subdirectory under PIPELINEOUTDIR
        cmdStr = cmdStr, #command-line to run
        resources = resources, #dict with maxmem (in GB) and maxtime (in secs)
        jobName= jobName, #non-unique job-type name
        inputFiles=inputFiles, #list of full paths
        filesToBeOutput=filesToBeOutput, #list of full paths, using MODULEOUTDIR
        filesToBeDeleted=filesToBeDeleted, #list of full paths
        caching='PipelineDefault', #True, False, PipelineDefault
        cleanUpJobFilesOnFail='PipelineDefault', #True, False, PipelineDefault
        cleanUpPipelineJobsOnFail='PipelineDefault', #True, False, PipelineDefault
        executionEngine='PipelineDefault' #lsf, local, mockPrint, PipelineDefault
    )
    print ('dispense ' + jobName)

    oxoq_file = filesToBeOutput_expanded[0]

    return oxoq_file






###########
#####tokens
###########
def tokens_pm(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, indiv_id, bam_tumor, bam_normal):

    #calling the module using the firehose names as the parameters:
    module_libdir = os.path.join(module_dir,'tokens')
    cmdStr = os.path.join(adaptor_dir,'run_module.py') #don't change
    cmdStr += ' --module_libdir ' + module_libdir + ' ' #don't change

    #Add in the parameters as the module is called from firehose:
    pileup_tokens = indiv_id
    pileup_tokens_moddir = '$MODULEOUTDIR/' + pileup_tokens + '.tok'
    cmdStr += ' '.join(['--bam', bam_normal , 
                        '--output_name', pileup_tokens, 
                        '--refdir',refdata_dir + '/public/hg19_fasta'])


    moduleSubDir = module_subdir

    resources = {'maxmem':8, 'maxtime':40000}; #memory and time limits
    jobName = 'tokens_A' #job name
    inputFiles = [bam_normal] #parameters from the function to be passed in.
    #change to include the variables you listed above to save from your function:
    filesToBeOutput = [pileup_tokens_moddir]

    #STOP CHANGING HERE
    #########################
    filesToBeDeleted = []

    filesToBeOutput_expanded = pipeline.dispense(
        moduleSubDir = moduleSubDir, #unique subdirectory under PIPELINEOUTDIR
        cmdStr = cmdStr, #command-line to run
        resources = resources, #dict with maxmem (in GB) and maxtime (in secs)
        jobName= jobName, #non-unique job-type name
        inputFiles=inputFiles, #list of full paths
        filesToBeOutput=filesToBeOutput, #list of full paths, using MODULEOUTDIR
        filesToBeDeleted=filesToBeDeleted, #list of full paths
        caching='PipelineDefault', #True, False, PipelineDefault
        cleanUpJobFilesOnFail='PipelineDefault', #True, False, PipelineDefault
        cleanUpPipelineJobsOnFail='PipelineDefault', #True, False, PipelineDefault
        executionEngine='PipelineDefault' #lsf, local, mockPrint, PipelineDefault
    )
    print ('dispense ' + jobName)
    
    #start changing again

    pileup_tokens_outfile = filesToBeOutput_expanded[0]
    return pileup_tokens_outfile


def div100_pm(pipeline,module_subdir,adaptor_dir, module_dir, refdata_dir,  indiv_id, contest_percent_value_file):

    out_module_path = os.path.join('$MODULEOUTDIR','contest_fraction_value.txt')

    cmdStr = " ".join(['python', os.path.join(module_dir,'div100','div100.py'), contest_percent_value_file, out_module_path])



    moduleSubDir = module_subdir
    resources = {'maxmem':1, 'maxtime':60};
    jobName = 'div100'
    inputFiles = [contest_percent_value_file]
    filesToBeOutput = [out_module_path]
    filesToBeDeleted = []

    filesToBeOutput_expanded = pipeline.dispense(
        moduleSubDir = moduleSubDir, #unique subdirectory under PIPELINEOUTDIR
        cmdStr = cmdStr, #command-line to run
        resources = resources, #dict with maxmem (in GB) and maxtime (in secs)
        jobName= jobName, #non-unique job-type name
        inputFiles=inputFiles, #list of full paths
        filesToBeOutput=filesToBeOutput, #list of full paths, using MODULEOUTDIR
        filesToBeDeleted=filesToBeDeleted, #list of full paths
        caching='PipelineDefault', #True, False, PipelineDefault
        cleanUpJobFilesOnFail='PipelineDefault', #True, False, PipelineDefault
        cleanUpPipelineJobsOnFail='PipelineDefault', #True, False, PipelineDefault
        executionEngine='PipelineDefault' #lsf, local, mockPrint, PipelineDefault
    )
    print ('dispense ' + jobName)

    contest_fraction_value_file = filesToBeOutput_expanded[0]
    return contest_fraction_value_file


def mutect_sg_pp(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, indiv_id, bam_tumor, bam_normal, contest_value_file, region, use_pcawg_contigs, intervals_filename=None):
    scatter_width = 200

    if use_pcawg_contigs:
        reference_file = 'human_g1k_v37_decoy.fasta'
    else:
        reference_file = 'Homo_sapiens_assembly19.fasta'
    
    normal_base_name = os.path.basename(bam_normal)
    normal_base_name = normal_base_name[:-4]
    tumor_base_name = os.path.basename(bam_tumor)
    tumor_base_name = tumor_base_name[:-4]
    ref_seq_fasta = os.path.join(refdata_dir,'public',reference_file)
    dbsnp_vcf = os.path.join(refdata_dir,'public','dbsnp_134_b37.leftAligned.vcf')
    #
    normal_panel = os.path.join(refdata_dir,'public','refseq_exome_10bp_hg19_300_1kg_normal_panel.vcf')
    cosmic_vcf = os.path.join(refdata_dir,'public','hg19_cosmic_v54_120711.vcf')
    

    prepare_args = {'scatter.jobs':str(scatter_width),
                    'normal.name':normal_base_name,
                    'normal.bam':bam_normal,
                    'tumor.name':tumor_base_name,
                    'tumor.bam':bam_tumor,
                    'reference.sequence.fasta':ref_seq_fasta,
                    'dbsnp.vcf':dbsnp_vcf,
                    'normal.panel':normal_panel,
                    'cosmic.vcf':cosmic_vcf,
                    'output.prefix':indiv_id,
                    'downsample.to.coverage':'10000',
                    }
    if contest_value_file is not None:
        # use contents of file as an arg, clipping off \n
        prepare_args['fraction.contamination'] = '`cat ' + contest_value_file + ' | awk 1 ORS=" "`'
    else:
        prepare_args['fraction.contamination'] = '0.02'
        
    
    if region == 'test':
        prepare_args['targets.interval.list'] = os.path.join(refdata_dir,'public','test.bed')  #tbd target.interval.list arg needs to be removed for production
    elif region == 'all':
        pass
    elif region == 'force_call':
        prepare_args['optional.parameter.1'] = '--force_output'
    else:
        raise Exception('unrecognized value for region arg: ' + region)


        
    
    module_callstats = os.path.join('$MODULEOUTDIR',indiv_id + '.call_stats.txt')
    wig_file = os.path.join('$MODULEOUTDIR',indiv_id + '.coverage.wig.txt.gz')
    
    module_libdir = os.path.join(module_dir,'mutect')
    input_files = [bam_normal,bam_tumor] # input to prepare
    if contest_value_file is not None:
        input_files.append(contest_value_file)
    if intervals_filename is not None:
        input_files.append(intervals_filename)
        prepare_args['targets.interval.list'] = intervals_filename
    else:
        pass
        #commented out for now, default interval list appears broken
        #wgs_interval_list = os.path.join(refdata_dir,'public','mutect_wgs_intervals.interval_list')
        #prepare_args['targets.interval.list'] = wgs_interval_list

    output_files = [module_callstats,wig_file] #output from gather
    sg_module_subdir = os.path.join(module_subdir,'sg')
    
    scatter_maxmem = 4
    gather_maxmem = 2

    expanded_output_files = scatter_gather_pp(pipeline, sg_module_subdir, module_libdir, adaptor_dir, input_files, output_files, prepare_args, scatter_width, scatter_maxmem, gather_maxmem)
    
    call_stats = expanded_output_files[0]
    coverage = expanded_output_files[1]
    return (call_stats,coverage)
    
    


def mutect1_sg_pp(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, indiv_id, bam_tumor, bam_normal, contest_fraction_value_file, region, use_pcawg_contigs):

    normal_base_name = os.path.basename(bam_normal)
    normal_base_name = normal_base_name[:-4]
    tumor_base_name = os.path.basename(bam_tumor)
    tumor_base_name = tumor_base_name[:-4]
    dbsnp_vcf = os.path.join(refdata_dir,'public','dbsnp_134_b37.leftAligned.vcf')
    #
    normal_panel = os.path.join(refdata_dir,'public','refseq_exome_10bp_hg19_300_1kg_normal_panel.vcf')
    cosmic_vcf = os.path.join(refdata_dir,'public','hg19_cosmic_v54_120711.vcf')

    if use_pcawg_contigs:
        reference_file = 'human_g1k_v37_decoy.fasta'
        interval_file = 'mutect_wgs_intervals.intervals'
        scatter_width = 195
    else:
        reference_file = 'Homo_sapiens_assembly19.fasta'
        interval_file = 'mutect_wgs_intervals.standard_contigs.intervals'
        scatter_width = 123
    ref_seq_fasta = os.path.join(refdata_dir,'public',reference_file)
    wgs_interval_list = os.path.join(refdata_dir,'public',interval_file)

    prepare_args = {'normal.name':normal_base_name,
                    'normal.bam':bam_normal,
                    'tumor.name':tumor_base_name,
                    'tumor.bam':bam_tumor,
                    'reference.sequence.fasta':ref_seq_fasta,
                    'dbsnp.vcf':dbsnp_vcf,
                    'normal.panel.vcf':normal_panel,
                    'cosmic.vcf':cosmic_vcf,
                    'output.prefix':indiv_id,
                    'downsample_to_coverage':'10000',
                    'targets.interval.list':wgs_interval_list
                    }
    if contest_fraction_value_file is not None:
        # use contents of file as an arg,
        prepare_args['fraction.contamination'] = '`cat ' + contest_fraction_value_file + '`'
        # # use contents of file as an arg, clipping off \n
        # prepare_args['fraction.contamination'] = '`cat ' + contest_value_file + ' | awk 1 ORS=" "`'
    else:
        prepare_args['fraction.contamination'] = '0.02'


    if region == 'test':
        prepare_args['targets.interval.list'] = os.path.join(refdata_dir,'public','test.bed')  #tbd target.interval.list arg needs to be removed for production
    elif region == 'all':
        pass
    elif region == 'dbsnp':
        #prepare_args['targets.interval.list'] = os.path.join(refdata_dir,'public','dbsnp_134_b37.leftAligned.intervals')
        prepare_args['targets.interval.list'] = os.path.join(refdata_dir,'public','1000genomes.common_variants.minor_allele_freq_5_percent.intervals')
        prepare_args['optional.parameter.1'] = '--force_output'
    else:
        raise Exception('unrecognized value for region arg: ' + region)


    module_callstats = os.path.join('$MODULEOUTDIR',indiv_id + '.call_stats.txt')
    wig_file = os.path.join('$MODULEOUTDIR',indiv_id + '.coverage.wig.txt')
    power_file = os.path.join('$MODULEOUTDIR',indiv_id + '.power.wig.txt')

    module_libdir = os.path.join(module_dir,'mutect1')
    input_files = [bam_normal,bam_tumor] # input to prepare
    if contest_fraction_value_file is not None:
        input_files.append(contest_fraction_value_file)
    output_files = [module_callstats,wig_file,power_file] #output from gather
    sg_module_subdir = os.path.join(module_subdir,'sg')

    scatter_maxmem = 4
    gather_maxmem = 2

    expanded_output_files = scatter_gather_pp(pipeline, sg_module_subdir, module_libdir, adaptor_dir, input_files, output_files, prepare_args, scatter_width, scatter_maxmem, gather_maxmem)

    call_stats = expanded_output_files[0]
    coverage = expanded_output_files[1]
    power = expanded_output_files[2]
    return (call_stats,coverage,power)


def mutect2_sg_pp(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, indiv_id, bam_tumor, bam_normal, contest_value_file, region, use_pcawg_contigs, intervals_filename=None):
    #scatter_width = 200
    #scatter_width = 626
    scatter_width = 482 ### match with targets_interval_list
    targets_interval_list = os.path.join(refdata_dir,'public','wgs_hg19_482_quasiuniform_intervals.bed')

    if use_pcawg_contigs:
        reference_file = 'human_g1k_v37_decoy.fasta'
    else:
        reference_file = 'Homo_sapiens_assembly19.fasta'
    
    normal_base_name = os.path.basename(bam_normal)
    normal_base_name = normal_base_name[:-4]
    tumor_base_name = os.path.basename(bam_tumor)
    tumor_base_name = tumor_base_name[:-4]
    ref_seq_fasta = os.path.join(refdata_dir,'public',reference_file)
    dbsnp_vcf = os.path.join(refdata_dir,'public','dbsnp_134_b37.leftAligned.vcf')
    #targets_interval_list = os.path.join(refdata_dir,'public','wgs_calling_regions.v1.time.M2order.bed')
    #
    normal_panel = os.path.join(refdata_dir,'public','refseq_exome_10bp_hg19_300_1kg_normal_panel.vcf')
    cosmic_vcf = os.path.join(refdata_dir,'public','hg19_cosmic_v54_120711.vcf')
    

    prepare_args = {'targets.interval.list':targets_interval_list,
                    'normal.name':normal_base_name,
                    'normal.bam':bam_normal,
                    'tumor.name':tumor_base_name,
                    'tumor.bam':bam_tumor,
                    'reference.sequence.fasta':ref_seq_fasta,
                    'dbsnp.vcf':dbsnp_vcf,
                    'cosmic.vcf':cosmic_vcf,
                    'normal.panel.vcf':normal_panel,
                    'output.prefix':indiv_id,
                    'maxReadsInMemoryPerSample':'1000'
                    }
    if contest_value_file is not None:
        # use contents of file as an arg, clipping off \n
        prepare_args['fraction.contamination'] = '`cat ' + contest_value_file + ' | awk 1 ORS=" "`'
    else:
        prepare_args['fraction.contamination'] = '0.02'
        
    
    if region == 'test':
        prepare_args['targets.interval.list'] = os.path.join(refdata_dir,'public','test.bed')  #tbd target.interval.list arg needs to be removed for production
    elif region == 'all':
        pass
    else:
        raise Exception('unrecognized value for region arg: ' + region)


        
    
    module_vcf = os.path.join('$MODULEOUTDIR',indiv_id + '.vcf')
    pass_vcf = os.path.join('$MODULEOUTDIR',indiv_id + '.PASS.vcf')
    
    module_libdir = os.path.join(module_dir,'M2_scatter')
    input_files = [bam_normal,bam_tumor] # input to prepare
    if contest_value_file is not None:
        input_files.append(contest_value_file)
    if intervals_filename is not None:
        input_files.append(intervals_filename)
        prepare_args['targets.interval.list'] = intervals_filename
    else:
        pass
        #commented out for now, default interval list appears broken
        #wgs_interval_list = os.path.join(refdata_dir,'public','mutect_wgs_intervals.interval_list')
        #prepare_args['targets.interval.list'] = wgs_interval_list

    output_files = [module_vcf,pass_vcf] #output from gather
    sg_module_subdir = os.path.join(module_subdir,'sg')
    
    scatter_maxmem = 4
    gather_maxmem = 4

    expanded_output_files = scatter_gather_pp(pipeline, sg_module_subdir, module_libdir, adaptor_dir, input_files, output_files, prepare_args, scatter_width, scatter_maxmem, gather_maxmem)
    
    M2_vcf = expanded_output_files[0]
    M2_pass_vcf = expanded_output_files[1]
    return (M2_vcf,M2_pass_vcf)
    
    


    
    
def coclean_pp(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, indiv_id, bam_tumor, bam_normal):
    module_subdir = 'input_bams_indiv'
    input_bams_fn = list_to_file_pm(pipeline,module_subdir,[bam_tumor, bam_normal],'input_bams.list')
    
    normal_base_name = os.path.basename(bam_normal)
    tumor_base_name = os.path.basename(bam_tumor)
    
    module_subdir = 'bam_map'
    str_list_list = [[tumor_base_name,tumor_base_name],[normal_base_name,normal_base_name]]
    map_file = list_of_lists_to_file_pm(pipeline,module_subdir,str_list_list,'input_bams.map')
    
    
    #module_subdir = 'tsvtolist'
    #(mapfile, listfile) = tsvtolist_pm(pipeline, module_subdir,adaptor_dir, module_dir, refdata_dir, indiv_id, input_bams_fn)
    
    module_subdir = 'RealignerTargetCreator_indiv'
    target_intervals = realignerTargetCreator_pp(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, indiv_id, input_bams_fn)

    module_subdir = 'IndelRealigner_indiv'
    (bam_cleaned_tumor, bam_cleaned_normal) = IndelRealigner_pp(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, indiv_id, input_bams_fn, target_intervals, bam_tumor, bam_normal, map_file)
    
    return (bam_cleaned_tumor, bam_cleaned_normal)
    
    
    
    
def realignerTargetCreator_pp(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, indiv_id, input_bams_fn):

#--module_libdir /cga/fh/pcawg_pipeline/modules/RealignerTargetCreator 
#--job.count 24 
#--input.bam /xchip/cga_home/gsaksena/prj/2015/docker_bringup_2015-01-26/target_gen_scatter_gather/input_bams.list 
#--reference.genome /seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta 
#--base.name pairname 
    scatter_width = 24
    prepare_args = {'job.count':str(scatter_width),
                    'input.bam':input_bams_fn,
                    'reference.genome':os.path.join(refdata_dir,'public/Homo_sapiens_assembly19.fasta'),
                    'base.name':indiv_id
                    }
    module_libdir = os.path.join(module_dir,'RealignerTargetCreator')
    input_files = [input_bams_fn] # input to prepare
    output_files = ['$MODULEOUTDIR/'+indiv_id+'.merged.intervals'] #output from gather
    sg_module_subdir = os.path.join(module_subdir,'sg')
    
    scatter_maxmem = 2
    gather_maxmem = 10

    expanded_output_files = scatter_gather_pp(pipeline, sg_module_subdir, module_libdir, adaptor_dir, input_files, output_files, prepare_args, scatter_width, scatter_maxmem, gather_maxmem)
    
    target_intervals = expanded_output_files[0]
    return target_intervals



def IndelRealigner_pp(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, indiv_id, input_bams_fn, target_intervals,bam_tumor, bam_normal, map_file):


    scatter_width = 50
    prepare_args = {'job.count':str(scatter_width),
                    'input.bam':input_bams_fn,
                    'reference.genome':os.path.join(refdata_dir,'public/Homo_sapiens_assembly19.fasta'),
                    'base.name':indiv_id,
                    'max.reads':'1000',
                    'merged.intervals':target_intervals,
                    'output.map':map_file 
                    }
    module_libdir = os.path.join(module_dir,'IndelRealigner')
    input_files = [target_intervals, input_bams_fn, map_file] # input to prepare
    output_files = ['$MODULEOUTDIR/'+os.path.basename(bam_tumor), '$MODULEOUTDIR/'+os.path.basename(bam_normal)] #output from gather
    sg_module_subdir = os.path.join(module_subdir,'sg')
    
    scatter_maxmem = 3
    gather_maxmem = 11

    expanded_output_files = scatter_gather_pp(pipeline, sg_module_subdir, module_libdir, adaptor_dir, input_files, output_files, prepare_args, scatter_width, scatter_maxmem, gather_maxmem)
    
    bam_cleaned_tumor = expanded_output_files[0]
    bam_cleaned_normal = expanded_output_files[1]
    return (bam_cleaned_tumor, bam_cleaned_normal)

#Needs to be called twice for normals and tumors
def dRangerPreprocess_pp(pipeline, module_subdir,adaptor_dir, module_dir, refdata_dir, indiv_id, bam,bam_type,save_intermediate_files):
    scatter_width = 24 #needs to be filled in with right number
    blacklist = os.path.join(refdata_dir,'public','lane_blacklist.txt')
    if bam_type is 'Normal':
        prepare_args = {    'bam':bam, #Putting arguments for both Normal PreProcess and Tumor Preprocess here
                            'blacklist':blacklist,
                            'refdir' :refdata_dir,
                            'cutoffSW':'150',
                            'cutoffSN':'450',
                            'MAX_ISZ':'2000',
                            'WEIRDPAIR_THRESHOLD':'10000',
                            'WEIRDPAIR_DISTINCTNESS_THRESHOLD':'2000',
                            'sample.id':indiv_id,
                            'save_intermediate_files':save_intermediate_files #TBD - set to False for productin
                        }
    elif bam_type is 'Tumor':
        prepare_args = {    'bam':bam, #Putting arguments for both Normal PreProcess and Tumor Preprocess here
                            'blacklist':blacklist,
                            'refdir' :refdata_dir,
                            'cutoffSW':'200',
                            'cutoffSN':'600',
                            'MAX_ISZ':'2000',
                            'WEIRDPAIR_THRESHOLD':'10000',
                            'WEIRDPAIR_DISTINCTNESS_THRESHOLD':'2000',
                            'sample.id':indiv_id,
                            'save_intermediate_files':save_intermediate_files #TBD - set to False for productin
                        }

    else:
        sys.exit('Invalid bam_type.  Bam_type for dRangerPreprocess must be Normal or Tumor.')

    module_libdir = os.path.join(module_dir,'dRangerPreprocess')
    input_files = [bam,blacklist] # input to prepare
    output_files = ['$MODULEOUTDIR/'+indiv_id+'.dRanger_input','$MODULEOUTDIR/'+indiv_id+'.all.isz'] #output from gather
    sg_module_subdir = os.path.join(module_subdir,'sg')

    scatter_maxmem = 1
    gather_maxmem = 24

    expanded_output_files = scatter_gather_pp(pipeline, sg_module_subdir, module_libdir, adaptor_dir, input_files, output_files, prepare_args, scatter_width, scatter_maxmem, gather_maxmem)

    dranger_input = expanded_output_files[0]
    breakpoints = expanded_output_files[1]
    return dranger_input,breakpoints


#This needs to be run on normals and tumors
def dRanger_Breakpointer_pp(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, indiv_id,bam,dRanger_results,isz):

    #calling the module using the firehose names as the parameters:
    module_libdir = os.path.join(module_dir,'dRanger_BreakPointer')
    cmdStr = os.path.join(adaptor_dir,'run_module.py') #don't change
    cmdStr += ' --module_libdir ' + module_libdir + ' ' #don't change

    blacklist = os.path.join(refdata_dir,'public','lane_blacklist.txt')
    #  number of scatter jobs
    jobsize=500
    scatter_jobs = 50
        
    #############################################################
    input_files = [blacklist,dRanger_results] # input to prepare
    output_files = ['$MODULEOUTDIR/'+indiv_id+'.breakpoints.txt','$MODULEOUTDIR/'+indiv_id+'.matched.sam'] #output from gather
    sg_module_subdir = os.path.join(module_subdir,'sg')

    scatter_maxmem = 6
    gather_maxmem = 2
    prepare_args = {
                        'sample.id':indiv_id,
                        'dRanger_results.txt':dRanger_results,
                        'bam':bam,
                        'lane.blacklist':blacklist,
                        'refdir' :os.path.join(refdata_dir,'public','hg19'),
                        'scatter_jobs':str(scatter_jobs),
                        'insertionsize':isz,
                        'fish_no_supp_reads_thres':'6',
                        'fish_low_confidence_fixsidewithread':'80',
                        'fish_low_confidence_fixsidewithoutread':'200',
                        'fish_high_confidence_fixsidewithread':'40',
                        'fish_high_confidence_fixsidewithoutread':'50',
                        'fish_allowedmm':'80',
                        'fish_tipsize':'15',
                        'fish_minmmintip':'5',
                        'fish_maxN':'10',
                        'fish_expand_pairs_extraction':'500',
                        'fish_maxreads':'50000',
                        'readlen':'101',
                        'align_enough_reads':'20',
                        'align_wsplit':'8',
                        'align_min_qual':'0.75'
                    }
    adaptor_dir = 'python3 ' + adaptor_dir
    expanded_output_files = scatter_gather_pp(pipeline, sg_module_subdir, module_libdir, adaptor_dir, input_files, output_files, prepare_args, scatter_jobs, scatter_maxmem, gather_maxmem)

    breakpoints = expanded_output_files[0]
    sam = expanded_output_files[1]
    return breakpoints,sam


#Needs tdir, output directory of dRangerPreprocess run on tumors as input paramter
#Also needs ndir
def dRangerRun_pp(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, indiv_id, bam_tumor, bam_normal,tdir,ndir,save_intermediate_files,use_PoN_JH):

    #calling the module using the firehose names as the parameters:
    module_libdir = os.path.join(module_dir,'dRangerRun')
    cmdStr = os.path.join(adaptor_dir,'run_module.py') #don't change
    cmdStr += ' --module_libdir ' + module_libdir + ' ' #don't change


    #############################################################
    #START CHANGING HERE:
    #Add in the parameters as the module is called from firehose:
    envstr = 'export LD_LIBRARY_PATH=$MCR_LD_LIBRARY_PATH; python3 '
    cmdStr = envstr + cmdStr
    cmdStr += ' '.join(['--individual',indiv_id,
                        '--tdir',tdir, #Needs tdir, output directory of dRangerPreprocess run on tumors as input paramter
                        '--ndir',ndir,#Needs ndir, output directory of dRangerPreprocess run on Normals as input paramter
                        '--tminmapq','5',
                        '--minpairs','2',
                        '--windowsize','2000',
                        '--nminwindow','2000',
                        '--normpaneldb',os.path.join(refdata_dir,'public','dRanger_PoN_JH'),
                        '--minsomratio','50',
                        '--nminspanfrac','0.5',
                        '--minscoreforbp','0.01',
                        '--min_ignore_matchingnormal','Inf', #Not sure if this is right way to designate infinity in python
                        '--build',os.path.join(refdata_dir,'public','R.mat'),
                        '--save_intermediate_files',save_intermediate_files,
			            '--use_PoN_JH',use_PoN_JH])




    moduleSubDir = module_subdir

    resources = {'maxmem':30, 'maxtime':14400}; #memory and time limits
    jobName = 'dRangerRun_A' #job name
    inputFiles = [tdir,ndir] #parameters from the function to be passed in.
    #change to include the variables you listed above to save from your function:
    filesToBeOutput = ['$MODULEOUTDIR/' + indiv_id + '.dRanger_results.forBP.txt','$MODULEOUTDIR/' + indiv_id + '.dRanger_results.mat','$MODULEOUTDIR/stdout.txt','$MODULEOUTDIR/stderr.txt']

    #STOP CHANGING HERE
    #########################
    filesToBeDeleted = []

    filesToBeOutput_expanded = pipeline.dispense(
        moduleSubDir = moduleSubDir, #unique subdirectory under PIPELINEOUTDIR
        cmdStr = cmdStr, #command-line to run
        resources = resources, #dict with maxmem (in GB) and maxtime (in secs)
        jobName= jobName, #non-unique job-type name
        inputFiles=inputFiles, #list of full paths
        filesToBeOutput=filesToBeOutput, #list of full paths, using MODULEOUTDIR
        filesToBeDeleted=filesToBeDeleted, #list of full paths
        caching='PipelineDefault', #True, False, PipelineDefault
        cleanUpJobFilesOnFail='PipelineDefault', #True, False, PipelineDefault
        cleanUpPipelineJobsOnFail='PipelineDefault', #True, False, PipelineDefault
        executionEngine='PipelineDefault' #lsf, local, mockPrint, PipelineDefault
    )
    print ('dispense ' + jobName)


    forbp = filesToBeOutput_expanded[0] #text file for breakpointer
    mat= filesToBeOutput_expanded[1] #mat file for dRangerFinalize
    stdout = filesToBeOutput_expanded[2]
    stderr = filesToBeOutput_expanded[3]
    return (forbp,mat,stdout,stderr)


def dRangerFinalize_pp(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, indiv_id,tbp,nbp,drrmat):

    #calling the module using the firehose names as the parameters:
    module_libdir = os.path.join(module_dir,'dRangerFinalize')
    cmdStr = os.path.join(adaptor_dir,'run_module.py') #don't change
    cmdStr += ' --module_libdir ' + module_libdir + ' ' #don't change


    #############################################################
    #START CHANGING HERE:
    #Add in the parameters as the module is called from firehose:
    envstr = 'export LD_LIBRARY_PATH=$MCR_LD_LIBRARY_PATH; python3 '
    cmdStr = envstr + cmdStr
    cmdStr += ' '.join(['--individual',indiv_id,
                        '--dRmatfile',drrmat,
                        '--circospng','None',
                        '--BPt',tbp,
                        '--BPn',nbp,
                        '--diffthresh','200'])




    #change these variables to the firehose parameters.
    #forbp = '$MODULEOUTDIR/' + indiv_id + '.dRanger_results.forBP.txt'
    #mat = '$MODULEOUTDIR/' + indiv_id + '.dRanger_results.mat'


    moduleSubDir = module_subdir

    resources = {'maxmem':8, 'maxtime':14400}; #memory and time limits
    jobName = 'dRangerFinalize_A' #job name
    inputFiles = [tbp,nbp,drrmat] #parameters from the function to be passed in.
    #change to include the variables you listed above to save from your function:
    filesToBeOutput = ['$MODULEOUTDIR/' + indiv_id + '.dRanger_results.detail.somatic.txt',
                       '$MODULEOUTDIR/' + indiv_id + '.dRanger_results.detail.all.mat',
                       '$MODULEOUTDIR/' + indiv_id + '.dRanger_results.detail.all.txt',
                       '$MODULEOUTDIR/' + indiv_id + '.dRanger_results.somatic.txt',
                       ]

    #STOP CHANGING HERE
    #########################
    filesToBeDeleted = []

    filesToBeOutput_expanded = pipeline.dispense(
        moduleSubDir = moduleSubDir, #unique subdirectory under PIPELINEOUTDIR
        cmdStr = cmdStr, #command-line to run
        resources = resources, #dict with maxmem (in GB) and maxtime (in secs)
        jobName= jobName, #non-unique job-type name
        inputFiles=inputFiles, #list of full paths
        filesToBeOutput=filesToBeOutput, #list of full paths, using MODULEOUTDIR
        filesToBeDeleted=filesToBeDeleted, #list of full paths
        caching='PipelineDefault', #True, False, PipelineDefault
        cleanUpJobFilesOnFail='PipelineDefault', #True, False, PipelineDefault
        cleanUpPipelineJobsOnFail='PipelineDefault', #True, False, PipelineDefault
        executionEngine='PipelineDefault' #lsf, local, mockPrint, PipelineDefault
    )
    print ('dispense ' + jobName)

    drr =filesToBeOutput_expanded[0]
    mat = filesToBeOutput_expanded[1]
    da = filesToBeOutput_expanded[2]
    ds = filesToBeOutput_expanded[3]

    return (drr,mat,da,ds)

def dRanger2VCF_pp(pipeline, module_subdir, adaptor_dir, module_dir, drrmatfinalized, refdata_dir, indiv_id,tsam,dRanger_reads,dRanger_reads_Normal):

    #calling the module using the firehose names as the parameters:
    module_libdir = os.path.join(module_dir,'dRanger2VCF')
    cmdStr = os.path.join(adaptor_dir,'run_module.py') #don't change
    cmdStr += ' --module_libdir ' + module_libdir + ' ' #don't change

    tumor_id = indiv_id + 'T'
    normal_id = indiv_id + 'N'

    #############################################################
    #START CHANGING HERE:
    #Add in the parameters as the module is called from firehose:
    envstr = 'export LD_LIBRARY_PATH=$MCR_LD_LIBRARY_PATH; python3 '
    cmdStr = envstr + cmdStr
    cmdStr += ' '.join(['--dRanger_results.detail.all.mat',drrmatfinalized,
                        '--ref_area',os.path.join(refdata_dir,'public','hg19'),
                        '--build','hg19',
                        '--minSomaticScore' ,'4',
                        '--minBPresult','0',
                        '--out.area','.',
                        '--pair_id',indiv_id,
                        '--Breakpointer.samfile',tsam,
                        '--Dranger_ReadFile',dRanger_reads,
                        '--Breakpointer.sufficient_splitreads','4',
                        '--tumor_id',tumor_id,
                        '--normal_id',normal_id,
                        '--Dranger_Normal_ReadFile',dRanger_reads_Normal])
                        


    #change these variables to the firehose parameters.
    #forbp = '$MODULEOUTDIR/' + indiv_id + '.dRanger_results.forBP.txt'
    #mat = '$MODULEOUTDIR/' + indiv_id + '.dRanger_results.mat'


    moduleSubDir = module_subdir
    # inflate maxtime for priority
    resources = {'maxmem':4, 'maxtime':1 * 3600}; #memory and time limits
    jobName = 'dRanger2VCF_A' #job name
    inputFiles = [drrmatfinalized,tsam,dRanger_reads,dRanger_reads_Normal] #parameters from the function to be passed in.
    #change to include the variables you listed above to save from your function:
    filesToBeOutput = ['$MODULEOUTDIR/' + indiv_id + '.dRanger_results.adjacency.vcf',
                       '$MODULEOUTDIR/' + indiv_id + '.dRanger_results.detail.txt',
                       ]
    
    #STOP CHANGING HERE
    #########################
    filesToBeDeleted = []

    filesToBeOutput_expanded = pipeline.dispense(
        moduleSubDir = moduleSubDir, #unique subdirectory under PIPELINEOUTDIR
        cmdStr = cmdStr, #command-line to run
        resources = resources, #dict with maxmem (in GB) and maxtime (in secs)
        jobName= jobName, #non-unique job-type name
        inputFiles=inputFiles, #list of full paths
        filesToBeOutput=filesToBeOutput, #list of full paths, using MODULEOUTDIR
        filesToBeDeleted=filesToBeDeleted, #list of full paths
        caching='PipelineDefault', #True, False, PipelineDefault
        cleanUpJobFilesOnFail='PipelineDefault', #True, False, PipelineDefault
        cleanUpPipelineJobsOnFail='PipelineDefault', #True, False, PipelineDefault
        executionEngine='PipelineDefault' #lsf, local, mockPrint, PipelineDefault
    )
    print ('dispense ' + jobName)

    vcf =filesToBeOutput_expanded[0]
    detail = filesToBeOutput_expanded[1]
 
    return (vcf, detail)


def getdRangerSupportingReads_pp(pipeline, module_subdir, adaptor_dir, module_dir, indiv_id,bam,sam,matfinalized,drrfinalized):


    read_scope = 'all' #'all' used for official pcawg pipeline and should give more accurate results, while 'somatic' should be faster.

    #calling the module using the firehose names as the parameters:
    module_libdir = os.path.join(module_dir,'getdRangerSupportingReads')
    cmdStr = os.path.join(adaptor_dir,'run_module.py') #don't change
    cmdStr += ' --module_libdir ' + module_libdir + ' ' #don't change


    #############################################################
    #START CHANGING HERE:
    #Add in the parameters as the module is called from firehose:
    cmdStr += ' '.join(['--pair_id',indiv_id,
                        '--tumor.bam',bam,
                        '--dRanger_results.all.matfile',matfinalized,
                        '--output_dir','.',
                        '--output_sam','dRanger.sam',
                        '--minmapq','5',
                        '--dRanger_results.all.or.somatic.',read_scope,
                        '--BreakPointer.SAM.file',sam,
                        '--FragLen_window','300',
                        '--SplitRead_window','10',
                        '--MaxReads_per_break','5000',
                        '--MinScore','2'])


    #change these variables to the firehose parameters.
    forbp = '$MODULEOUTDIR/' + indiv_id + '.dRanger_results.forBP.txt'
    mat = '$MODULEOUTDIR/' + indiv_id + '.dRanger_results.mat'


    moduleSubDir = module_subdir

    resources = {'maxmem':8, 'maxtime':14400}; #memory and time limits
    jobName = 'getdRangerSupportingReads_A' #job name
    inputFiles = [bam,matfinalized,sam,drrfinalized] #parameters from the function to be passed in.
    #change to include the variables you listed above to save from your function:
    #filesToBeOutput = ['$MODULEOUTDIR/' + indiv_id + '.dRanger.supporting_reads.txt']
    filesToBeOutput = ['$MODULEOUTDIR/' + indiv_id + '.dRanger.all_reads.txt']

    #STOP CHANGING HERE
    #########################
    filesToBeDeleted = []

    filesToBeOutput_expanded = pipeline.dispense(
        moduleSubDir = moduleSubDir, #unique subdirectory under PIPELINEOUTDIR
        cmdStr = cmdStr, #command-line to run
        resources = resources, #dict with maxmem (in GB) and maxtime (in secs)
        jobName= jobName, #non-unique job-type name
        inputFiles=inputFiles, #list of full paths
        filesToBeOutput=filesToBeOutput, #list of full paths, using MODULEOUTDIR
        filesToBeDeleted=filesToBeDeleted, #list of full paths
        caching='PipelineDefault', #True, False, PipelineDefault
        cleanUpJobFilesOnFail='PipelineDefault', #True, False, PipelineDefault
        cleanUpPipelineJobsOnFail='PipelineDefault', #True, False, PipelineDefault
        executionEngine='PipelineDefault' #lsf, local, mockPrint, PipelineDefault
    )
    print ('dispense ' + jobName)

    reads = filesToBeOutput_expanded[0]

    return (reads)


def tabix_pm(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, indiv_id, vcf_file, filetag):
    #calling the module using the firehose names as the parameters:
    module_libdir = os.path.join(module_dir,'tabix')
    cmdStr = os.path.join(adaptor_dir,'run_module.py') #don't change
    cmdStr += ' --module_libdir ' + module_libdir + ' ' #don't change

    #Add in the parameters as the module is called from firehose:
    cmdStr += ' '.join([
                        '--input_file', vcf_file,
                        '--output_base_name', indiv_id,
                        '--output_extension', '.' + filetag+'.vcf'])

    #change these variables to the firehose parameters.
    vcf_gz = '$MODULEOUTDIR/' + indiv_id + '.' + filetag + '.vcf.gz'
    vcf_gz_tbi = '$MODULEOUTDIR/' + indiv_id + '.' + filetag + '.vcf.gz.tbi'

    moduleSubDir = module_subdir

    resources = {'maxmem':1, 'maxtime':14400}; #memory and time limits
    jobName = 'tabix' #job name
    inputFiles = [vcf_file] #parameters from the function to be passed in.
    #change to include the variables you listed above to save from your function:
    filesToBeOutput = [vcf_gz, vcf_gz_tbi]

    #STOP CHANGING HERE
    #########################
    filesToBeDeleted = []

    filesToBeOutput_expanded = pipeline.dispense(
        moduleSubDir = moduleSubDir, #unique subdirectory under PIPELINEOUTDIR
        cmdStr = cmdStr, #command-line to run
        resources = resources, #dict with maxmem (in GB) and maxtime (in secs)
        jobName= jobName, #non-unique job-type name
        inputFiles=inputFiles, #list of full paths
        filesToBeOutput=filesToBeOutput, #list of full paths, using MODULEOUTDIR
        filesToBeDeleted=filesToBeDeleted, #list of full paths
        caching='PipelineDefault', #True, False, PipelineDefault
        cleanUpJobFilesOnFail='PipelineDefault', #True, False, PipelineDefault
        cleanUpPipelineJobsOnFail='PipelineDefault', #True, False, PipelineDefault
        executionEngine='PipelineDefault' #lsf, local, mockPrint, PipelineDefault
    )

    vcf_gz_outfile = filesToBeOutput_expanded[0]
    vcf_gz_tbi_outfile = filesToBeOutput_expanded[1]
    return (vcf_gz_outfile, vcf_gz_tbi_outfile)




def tsvtolist_pm(pipeline, module_subdir,adaptor_dir, module_dir, refdata_dir, indiv_id, input_bams_fn):

    #calling the module using the firehose names as the parameters:
    module_libdir = os.path.join(module_dir,'tsvtolist')
    cmdStr = os.path.join(adaptor_dir,'run_module.py') #don't change
    cmdStr += ' --module_libdir ' + module_libdir + ' ' #don't change


    #Add in the parameters as the module is called from firehose:
    cmdStr += ' '.join(['--tsv.file', input_bams_fn, 
                        '--base.name', indiv_id, 
                        '--map.modifier', '.cleaned.bam', 
                        '--list.extension', 'list', 
                        '--map.extension', 'map', 
                        ])
        
    moduleSubDir = module_subdir
    mapfile_moddir = '$MODULEOUTDIR/'+indiv_id+'.map'
    listfile_moddir = '$MODULEOUTDIR/'+indiv_id+'.list'
    
    resources = {'maxmem':1, 'maxtime':14400}; #memory and time limits
    jobName = 'tsvtolist' #job name
    inputFiles = [input_bams_fn] #parameters from the function to be passed in. 
    #change to include the variables you listed above to save from your function:
    filesToBeOutput = [mapfile_moddir, listfile_moddir]
    
    #STOP CHANGING HERE
    #########################
    filesToBeDeleted = []
    
    filesToBeOutput_expanded = pipeline.dispense(
        moduleSubDir = moduleSubDir, #unique subdirectory under PIPELINEOUTDIR
        cmdStr = cmdStr, #command-line to run
        resources = resources, #dict with maxmem (in GB) and maxtime (in secs)
        jobName= jobName, #non-unique job-type name
        inputFiles=inputFiles, #list of full paths
        filesToBeOutput=filesToBeOutput, #list of full paths, using MODULEOUTDIR
        filesToBeDeleted=filesToBeDeleted, #list of full paths
        caching='PipelineDefault', #True, False, PipelineDefault
        cleanUpJobFilesOnFail='PipelineDefault', #True, False, PipelineDefault
        cleanUpPipelineJobsOnFail='PipelineDefault', #True, False, PipelineDefault
        executionEngine='PipelineDefault' #lsf, local, mockPrint, PipelineDefault
    )
    
    mapfile = filesToBeOutput_expanded[0]
    listfile = filesToBeOutput_expanded[1]
    return (mapfile, listfile)

def list_of_lists_to_file_pm(pipeline,list_to_file_module_subdir,str_list_list,out_fn):
    cmdStr = ''
    
    num_cols = len(str_list_list[0])
    fns = []
    for colnum in range(num_cols):
        fn = 'col' + str(colnum)
        fns.append(fn)
        for s_list in str_list_list:
            cmdStr += 'echo %s >> %s; '%(s_list[colnum],fn)
    cmdStr += 'paste ' + ' '.join(fns) + ' > ' + out_fn
    
    moduleSubDir = list_to_file_module_subdir
    resources = {'maxmem':1, 'maxtime':300};
    jobName = 'list_to_file'
    inputFiles = []
    filesToBeOutput = ['$MODULEOUTDIR/'+out_fn]
    filesToBeDeleted = []
    
    filesToBeOutput_expanded = pipeline.dispense(
        moduleSubDir = moduleSubDir, #unique subdirectory under PIPELINEOUTDIR
        cmdStr = cmdStr, #command-line to run
        resources = resources, #dict with maxmem (in GB) and maxtime (in secs)
        jobName= jobName, #non-unique job-type name
        inputFiles=inputFiles, #list of full paths
        filesToBeOutput=filesToBeOutput, #list of full paths, using MODULEOUTDIR
        filesToBeDeleted=filesToBeDeleted, #list of full paths
        caching='PipelineDefault', #True, False, PipelineDefault
        cleanUpJobFilesOnFail='PipelineDefault', #True, False, PipelineDefault
        cleanUpPipelineJobsOnFail='PipelineDefault', #True, False, PipelineDefault
        executionEngine='PipelineDefault' #lsf, local, mockPrint, PipelineDefault
    )
    print ('dispense ' + jobName)

    outpath = filesToBeOutput_expanded[0]
    return outpath

def list_to_file_pm(pipeline,list_to_file_module_subdir,str_list,fn):
    cmdStr = ''
    for s in str_list:
        cmdStr += 'echo %s >> %s; '%(s,fn)
    moduleSubDir = list_to_file_module_subdir
    resources = {'maxmem':1, 'maxtime':300};
    jobName = 'list_to_file'
    inputFiles = []
    filesToBeOutput = ['$MODULEOUTDIR/'+fn]
    filesToBeDeleted = []
    
    filesToBeOutput_expanded = pipeline.dispense(
        moduleSubDir = moduleSubDir, #unique subdirectory under PIPELINEOUTDIR
        cmdStr = cmdStr, #command-line to run
        resources = resources, #dict with maxmem (in GB) and maxtime (in secs)
        jobName= jobName, #non-unique job-type name
        inputFiles=inputFiles, #list of full paths
        filesToBeOutput=filesToBeOutput, #list of full paths, using MODULEOUTDIR
        filesToBeDeleted=filesToBeDeleted, #list of full paths
        caching='PipelineDefault', #True, False, PipelineDefault
        cleanUpJobFilesOnFail='PipelineDefault', #True, False, PipelineDefault
        cleanUpPipelineJobsOnFail='PipelineDefault', #True, False, PipelineDefault
        executionEngine='PipelineDefault' #lsf, local, mockPrint, PipelineDefault
    )
    print ('dispense ' + jobName)

    outfn = filesToBeOutput_expanded[0]
    return outfn
    
def initialize_pipeline(communicationDirBase, pipelineOutDir, pipelinePriority, fail_mode):
    if fail_mode == 'abort_on_fail':
        cleanUpPipelineJobsOnFail = "True"
    elif fail_mode == 'continue_on_fail':
        cleanUpPipelineJobsOnFail = "False"
    else:
        raise Exception('unexpected value for failMode: %s'%failMode)
    pipeline = pipetteClient.PipetteClient(
        pipelineOutDir, #root dir for pipeline output, must contain the word Pipette somewhere
        pipelineName='pcawg_pipeline', #non-unique pipeline name
        defaultCaching='True',  #True, False
        defaultCleanUpJobFilesOnFail='False', #True, False
        defaultCleanUpPipelineJobsOnFail=cleanUpPipelineJobsOnFail, #True, False
        defaultExecutionEngine='local', #lsf, local, mockPrint
        pipelinePriority=pipelinePriority,  #integer 0-100, higher value = more important
        injectionMap={}, #substitution dict for dispense parameters, plus MODULEOUTDIR and PIPELINEOUTDIR
        communicationDirBase=communicationDirBase #must match corresponding Pipette Server
        )
    return pipeline


def scatter_gather_pp(pipeline, module_subdir, module_libdir, adaptor_dir, input_files, output_files, prepare_args, scatter_width, scatter_maxmem, gather_maxmem):
    
    prepare_module_subdir = os.path.join(module_subdir,'prepare')
    prepare_out_file = sg_prepare_pm(pipeline, prepare_module_subdir, module_libdir, adaptor_dir, input_files, prepare_args)
    
    scatter_stdouts=[]
    for i in range(scatter_width):
        jobNumber = i+1
        scatter_module_subdir = os.path.join(module_subdir,"scatter.%010d" % (jobNumber))
        scatter_stdout = sg_scatter_pm(pipeline, scatter_module_subdir, jobNumber, module_libdir, adaptor_dir, prepare_out_file, scatter_maxmem)
        scatter_stdouts.append(scatter_stdout)
    
    gather_module_subdir = os.path.join(module_subdir,'gather')

    filesToBeOutput_expanded = sg_gather_pm(pipeline, gather_module_subdir, module_libdir, adaptor_dir, prepare_out_file, scatter_stdouts, gather_maxmem, output_files)
    return filesToBeOutput_expanded

def sg_prepare_pm(pipeline, prepare_module_subdir, module_libdir, adaptor_dir, input_files, prepare_args):

    cmdStr = os.path.join(adaptor_dir,'run_sg_prepare.py')
    cmdStr += ' --module_libdir ' + module_libdir 
    for param in prepare_args:
        cmdStr += ' --' + param + ' ' + prepare_args[param]
        
    moduleSubDir = prepare_module_subdir
        
    resources = {'maxmem':1, 'maxtime':300};
    jobName = 'prepare'
    inputFiles = input_files
    filesToBeOutput = ['$MODULEOUTDIR/prepareResults.out']
    filesToBeDeleted = []
    
    filesToBeOutput_expanded = pipeline.dispense(
        moduleSubDir = moduleSubDir, #unique subdirectory under PIPELINEOUTDIR
        cmdStr = cmdStr, #command-line to run
        resources = resources, #dict with maxmem (in GB) and maxtime (in secs)
        jobName= jobName, #non-unique job-type name
        inputFiles=inputFiles, #list of full paths
        filesToBeOutput=filesToBeOutput, #list of full paths, using MODULEOUTDIR
        filesToBeDeleted=filesToBeDeleted, #list of full paths
        caching='PipelineDefault', #True, False, PipelineDefault
        cleanUpJobFilesOnFail='PipelineDefault', #True, False, PipelineDefault
        cleanUpPipelineJobsOnFail='PipelineDefault', #True, False, PipelineDefault
        executionEngine='PipelineDefault' #lsf, local, mockPrint, PipelineDefault
    )
    print ('dispense ' + jobName)

    outfn = filesToBeOutput_expanded[0]
    return outfn


def sg_scatter_pm(pipeline, scatter_module_subdir, jobNumber, module_libdir, adaptor_dir, prepare_out_file, scatter_maxmem):

    cmdStr = os.path.join(adaptor_dir,'run_sg_scatter.py')
    cmdStr += ' --module_libdir ' + module_libdir 
    cmdStr += ' --prepare_outdir ' + os.path.dirname(prepare_out_file) 
    cmdStr += ' --scatter_index ' + str(jobNumber) 
        
    moduleSubDir = scatter_module_subdir
        
    resources = {'maxmem':scatter_maxmem, 'maxtime':14400};
    jobName = 'scatter'
    inputFiles = [prepare_out_file]
    filesToBeOutput = ['$MODULEOUTDIR/stdout.txt']
    filesToBeDeleted = []
    
    filesToBeOutput_expanded = pipeline.dispense(
        moduleSubDir = moduleSubDir, #unique subdirectory under PIPELINEOUTDIR
        cmdStr = cmdStr, #command-line to run
        resources = resources, #dict with maxmem (in GB) and maxtime (in secs)
        jobName= jobName, #non-unique job-type name
        inputFiles=inputFiles, #list of full paths
        filesToBeOutput=filesToBeOutput, #list of full paths, using MODULEOUTDIR
        filesToBeDeleted=filesToBeDeleted, #list of full paths
        caching='PipelineDefault', #True, False, PipelineDefault
        cleanUpJobFilesOnFail='PipelineDefault', #True, False, PipelineDefault
        cleanUpPipelineJobsOnFail='PipelineDefault', #True, False, PipelineDefault
        executionEngine='PipelineDefault' #lsf, local, mockPrint, PipelineDefault
    )
    print ('dispense ' + jobName)

    
    outfn = filesToBeOutput_expanded[0]
    return outfn


def sg_gather_pm(pipeline, gather_module_subdir, module_libdir, adaptor_dir, prepare_out_file, scatter_stdouts, gather_maxmem, output_files):
    
    cmdStr = os.path.join(adaptor_dir,'run_sg_gather.py')
    cmdStr += ' --module_libdir ' + module_libdir 
    cmdStr += ' --prepare_outdir ' + os.path.dirname(prepare_out_file) 
    for scatter_stdout in scatter_stdouts:
        cmdStr += ' --scatter_outdir ' + os.path.dirname(scatter_stdout)
        
    moduleSubDir = gather_module_subdir
        
    resources = {'maxmem':gather_maxmem, 'maxtime':14400};
    jobName = 'gather'
    inputFiles = scatter_stdouts
    filesToBeOutput = output_files
    filesToBeDeleted = []
    
    filesToBeOutput_expanded = pipeline.dispense(
        moduleSubDir = moduleSubDir, #unique subdirectory under PIPELINEOUTDIR
        cmdStr = cmdStr, #command-line to run
        resources = resources, #dict with maxmem (in GB) and maxtime (in secs)
        jobName= jobName, #non-unique job-type name
        inputFiles=inputFiles, #list of full paths
        filesToBeOutput=filesToBeOutput, #list of full paths, using MODULEOUTDIR
        filesToBeDeleted=filesToBeDeleted, #list of full paths
        caching='PipelineDefault', #True, False, PipelineDefault
        cleanUpJobFilesOnFail='PipelineDefault', #True, False, PipelineDefault
        cleanUpPipelineJobsOnFail='PipelineDefault', #True, False, PipelineDefault
        executionEngine='PipelineDefault' #lsf, local, mockPrint, PipelineDefault
    )
    print ('dispense ' + jobName)

    #returns a list of files in the same order as were passed in as output_files.
    return filesToBeOutput_expanded


def make_external_links_pm(pipeline, module_subdir, module_dir, link_dir, data_files):

    module_libdir = os.path.join(module_dir,'make_external_links')

    cmdStr = 'python %s/make_external_links.py $PIPELINEOUTDIR %s '%(module_libdir, link_dir)
    cmdStr += ' '.join(data_files)


    moduleSubDir = module_subdir

    resources = {'maxmem':1, 'maxtime':14400};
    jobName = 'make_external_links'
    inputFiles = data_files
    filesToBeOutput = [] #Nothing listed, to permit the links to land outside the tree
    filesToBeDeleted = []

    filesToBeOutput_expanded = pipeline.dispense(
        moduleSubDir = moduleSubDir, #unique subdirectory under PIPELINEOUTDIR
        cmdStr = cmdStr, #command-line to run
        resources = resources, #dict with maxmem (in GB) and maxtime (in secs)
        jobName= jobName, #non-unique job-type name
        inputFiles=inputFiles, #list of full paths
        filesToBeOutput=filesToBeOutput, #list of full paths, using MODULEOUTDIR
        filesToBeDeleted=filesToBeDeleted, #list of full paths
        caching='PipelineDefault', #True, False, PipelineDefault
        cleanUpJobFilesOnFail='PipelineDefault', #True, False, PipelineDefault
        cleanUpPipelineJobsOnFail='PipelineDefault', #True, False, PipelineDefault
        executionEngine='PipelineDefault' #lsf, local, mockPrint, PipelineDefault
    )
    print ('dispense ' + jobName)

    #returns a list of files in the same order as were passed in as output_files.
    return None



if __name__ == '__main__':

    communicationDirBase = sys.argv[1]
    pipelineOutDir = sys.argv[2]
    indiv_id = sys.argv[3]
    bam_tumor = sys.argv[4]
    bam_normal = sys.argv[5]
    sub_workflow = sys.argv[6]
    refdata_dir = sys.argv[7]

    fail_mode = 'continue_on_fail' #abort_on_fail or continue_on_fail

    adaptor_dir = '/opt/src/algutil/firehose_module_adaptor'
    module_dir = '/opt/src/modules'
    tmp_dir = pipelineOutDir + '/tmp'
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)
    pipelinePriority = 50

    #make historical tmpdir also
    if not os.path.exists('/cga/fh/pcawg_pipeline'):
        os.makedirs('/cga/fh/pcawg_pipeline')
        os.symlink(tmp_dir, '/cga/fh/pcawg_pipeline/tmp')

    gnos_outdir = '$PIPELINEOUTDIR/links_for_gnos'
    broad_outdir = '$PIPELINEOUTDIR/links_for_broad'

    run_pcawg_pipeline(communicationDirBase, pipelineOutDir, pipelinePriority, adaptor_dir, module_dir, refdata_dir, 
        tmp_dir, gnos_outdir, broad_outdir, indiv_id, bam_tumor, bam_normal, sub_workflow, fail_mode)


