#!/usr/bin/python

import os
import sys
import collections





import pipetteClient

def run_hc_pipeline(communicationDirBase, pipelineOutDir, pipelinePriority, adaptor_dir, module_dir, refdata_dir, tmp_dir, broad_outdir, indiv_id, bam_normal, sub_workflow, fail_mode):

    #count number of CPU's
    fid = open('/proc/cpuinfo')
    cpu_count = 0
    for line in fid:
        if 'processor' in line:
            cpu_count = cpu_count + 1
    fid.close()

    pipeline = initialize_pipeline(communicationDirBase, pipelineOutDir, pipelinePriority, fail_mode)
    
    annotations = hc_pipeline_pp(pipeline, adaptor_dir, module_dir, refdata_dir, tmp_dir, broad_outdir, indiv_id, bam_normal, sub_workflow, cpu_count)
    
    pipeline.go()


def hc_pipeline_pp (pipeline, adaptor_dir, module_dir, refdata_dir, tmp_dir, broad_outdir,
                    indiv_id, bam_normal, sub_workflow, cpu_count):

    module_subdir = 'extract_bam_id'
    bam_id_file = extract_bam_id_pp(pipeline, module_subdir, indiv_id, adaptor_dir, module_dir, refdata_dir, bam_normal)
    make_external_links_pm(pipeline,'extract_bam_id_extlinks',module_dir, broad_outdir, [bam_id_file])

    module_subdir = 'haplotypecaller_sg'
    (haplotype_caller_gvcf, haplotype_caller_tbi) = haplotype_caller_sg_pp(pipeline, module_subdir, adaptor_dir, module_dir, refdata_dir, indiv_id, bam_normal)
    #The vcf files from haplotype caller are not ready for PCAWG consumption, they need further processing via VQSR
    make_external_links_pm(pipeline, 'haplotypecaller_sg_extlinks' ,module_dir, broad_outdir, (haplotype_caller_gvcf, haplotype_caller_tbi))

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
    bam_normal = sys.argv[4]
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

    run_hc_pipeline(communicationDirBase, pipelineOutDir, pipelinePriority, adaptor_dir, module_dir, refdata_dir, 
        tmp_dir, broad_outdir, indiv_id, bam_normal, sub_workflow, fail_mode)


