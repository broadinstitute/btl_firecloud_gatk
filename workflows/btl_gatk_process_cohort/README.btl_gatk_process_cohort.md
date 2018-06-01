# README for Process Cohort

Important! All examples below shown use widdler_dev. Please only use this version until the gatk2cloud changes have been pushed to production.

## TLDR

0. Log in to Stout or an interactive UGES host (widdler_dev will only work here).
1. Prepare your json input file (see examples below). You can use the process_samples gvcf output with GS urls 
for the value for the gatk_process_cohort.gvcf_fofn json key.
2. Run widdler upload like so:

```sh /cil/shed/apps/internal/widdler_dev/widdler.sh upload /cil/shed/apps/internal/wdl_cloud/taskdef.btl_gatk_process_cohort.wdl </path/to/json/input>```

3. Run widdler run like so:

```sh /cil/shed/apps/internal/widdler_dev/widdler.sh run /cil/shed/apps/internal/wdl_cloud/taskdef.btl_gatk_process_cohort.wdl </path/to/json/input> -S gscid-cloud -d /cil/shed/apps/internal/wdl_cloud/wdl_bundle_process_cohort.zip ```

4. Wait for widdler email.

5. If successful, find your downloaded files in the directory specified in the json input.

## Introduction

The Process Cohort portion of the GATK2CLOUD pipeline allows users to perform the following
steps, in sequence:

1. Call joint genotyping on a cohort of gvcf files.
2. Optionally run VQSR (variant quality score recalibration) on the joint calling output.
3. Call variant filtration on joint genotyping (or VQSR output if VQSR is used)
4. Optionally run filter genotypes.
5. Call snpeff.

## JSON Input File Format

The following is a template you can use to create the input json file for the workflow:
```
{
  "gatk_process_cohort.onprem_download_path": "(optional) String?",
  "gatk_process_cohort.snpeff_db_name": "String",
  "gatk_process_cohort.snp_max_gaussians": "Int",
  "gatk_process_cohort.snpeff_db_tgz": "File",
  "gatk_process_cohort.cohort_name": "String",
  "gatk_process_cohort.indel_max_gaussians": "Int",
  "gatk_process_cohort.known_sites_vcfs": "Array[File]",
  "gatk_process_cohort.reference_tgz": "File",
  "gatk_process_cohort.handoff_files": "(optional) Map[String, String]?",
  "gatk_process_cohort.mq_cap_snp": "Int",
  "gatk_process_cohort.run_gatk_vqsr": "(optional) Boolean",
  "gatk_process_cohort.indel_annotation": "Array[String]",
  "gatk_process_cohort.gvcf_fofn": "File",
  "gatk_process_cohort.gatk_joint_genotype_task.extra_gg_params": "(optional) String?",
  "gatk_process_cohort.indel_resource_params": "Array[String]",
  "gatk_process_cohort.gatk_joint_genotype_task.all_sites": "(optional) Boolean?",
  "gatk_process_cohort.gatk_snpeff_task.snpeff_extra_params": "(optional) String?",
  "gatk_process_cohort.ts_filter_indel": "Float",
  "gatk_process_cohort.snp_resource_params": "Array[String]",
  "gatk_process_cohort.mq_cap_indel": "Int",
  "gatk_process_cohort.extra_vr_params": "(optional) String?",
  "gatk_process_cohort.run_gatk_filter_genotypes": "Boolean",
  "gatk_process_cohort.known_sites_vcf_tbis": "Array[File]",
  "gatk_process_cohort.snp_annotation": "Array[String]",
  "gatk_process_cohort.ts_filter_snp": "Float",
  "gatk_process_cohort.output_disk_gb": "String"
}
```
Notes: 
* Optional items can be omitted. 
* If not using VQSR, write false as its value.
* Due to a quirk with the WDL language, the two know_sites parameters can't be set to optional even though the tasks
that use them are being made optional. So if not using known_sites, just pass an empty array([]). 

An example of a workflow input that is not using VQSR/known_sites and related parameters. Note that omitting the 
gatk_process_cohort.run_gatk_vqsr line from the file is equivalent to setting it to false:
```
{
  "gatk_process_cohort.onprem_download_path": "/cil/shed/sandboxes/amr/gatk_cloud/output",
  "gatk_process_cohort.run_gatk_vqsr": false,
  "gatk_process_cohort.ts_filter_snp": 0.0,
  "gatk_process_cohort.ts_filter_indel": 0.0,
  "gatk_process_cohort.snpeff_db_name": "CNA2",
  "gatk_process_cohort.snpeff_db_tgz": "/cil/shed/sandboxes/amr/gatk_cloud/CNA2.tgz",
  "gatk_process_cohort.snp_max_gaussians": 0,
  "gatk_process_cohort.snp_resource_params": [],
  "gatk_process_cohort.indel_resource_params": [],
  "gatk_process_cohort.cohort_name": "CandidaAurisCohort",
  "gatk_process_cohort.indel_max_gaussians": 0,
  "gatk_process_cohort.known_sites_vcfs":[],
  "gatk_process_cohort.reference_tgz": "gs://4b66fc8a-tmp/cromwell-executions/gatk_indexref/14d24780-51c9-4932-bdd9-9562c4f7d493/call-gatk_indexref_task/minion_illumina_hybrid_clean_MT.tgz",
  "gatk_process_cohort.mq_cap_snp": 0,
  "gatk_process_cohort.mq_cap_indel": 0,
  "gatk_process_cohort.indel_annotation": [],
  "gatk_process_cohort.gvcf_fofn": "/cil/shed/sandboxes/amr/gatk_cloud/process_cohort.fofn",
  "gatk_process_cohort.run_gatk_filter_genotypes": true,
  "gatk_process_cohort.known_sites_vcf_tbis": [],
  "gatk_process_cohort.snp_annotation": []
}
```

An example of a workflow that does use VQSR/known_sites is as follows. Note that the other optional tasks (snpeff, filter_genotypes)
are not specified, and such, will not run:
```
{
  "gatk_process_cohort.onprem_download_path": "/cil/shed/sandboxes/amr/gatk_cloud/output",
  "gatk_process_cohort.run_gatk_vqsr": true,
  "gatk_process_cohort.ts_filter_snp": 99.5,
  "gatk_process_cohort.ts_filter_indel": 99.0,
  "gatk_process_cohort.snpeff_db_name": "",
  "gatk_process_cohort.snpeff_db_tgz": "",
  "gatk_process_cohort.snp_max_gaussians": "2",
  "gatk_process_cohort.cohort_name": "PlasmodiumFalciparum",
  "gatk_process_cohort.indel_max_gaussians": "2",
  "gatk_process_cohort.reference_tgz": "gs://4b66fc8a-tmp/cromwell-executions/gatk_indexref/14d24780-51c9-4932-bdd9-9562c4f7d493/call-gatk_indexref_task/PlasmoDB-28_Pfalciparum3D7_Genome.tgz",
  "gatk_process_cohort.mq_cap_snp": "70",
  "gatk_process_cohort.mq_cap_indel": "70",
  "gatk_process_cohort.indel_annotation": [
      "QD",
      "FS",
      "MQ"
   ],
   gatk_process_cohort.snp_resource_params":[  
      "7g8_gb4,known=false,training=true,truth=true,prior=15.0 /gsap/garage-protistvector/U19_Aim4/Pf3K/7g8_gb4.combined.final.vcf.gz",
      "hb3_dd2,known=false,training=true,truth=true,prior=15.0 /gsap/garage-protistvector/U19_Aim4/Pf3K/hb3_dd2.combined.final.vcf.gz",
      "3d7_hb3,known=false,training=true,truth=true,prior=15.0 /gsap/garage-protistvector/U19_Aim4/Pf3K/3d7_hb3.combined.final.vcf.gz"
   ],
   "gatk_process_cohort.indel_resource_params":[  
      "7g8_gb4,known=false,training=true,truth=true,prior=12.0 /gsap/garage-protistvector/U19_Aim4/Pf3K/7g8_gb4.combined.final.vcf.gz",
      "hb3_dd2,known=false,training=true,truth=true,prior=12.0 /gsap/garage-protistvector/U19_Aim4/Pf3K/hb3_dd2.combined.final.vcf.gz",
      "3d7_hb3,known=false,training=true,truth=true,prior=12.0 /gsap/garage-protistvector/U19_Aim4/Pf3K/3d7_hb3.combined.final.vcf.gz"
   ],
  "gatk_process_cohort.known_sites_vcf_tbis": [
                "gs://broad-cil-devel-bucket/input_data/7g8_gb4.combined.final.vcf.gz.tbi",
                "gs://broad-cil-devel-bucket/input_data/hb3_dd2.combined.final.vcf.gz.tbi",
                "gs://broad-cil-devel-bucket/input_data/3d7_hb3.combined.final.vcf.gz.tbi"
  ],
  "gatk_process_cohort.known_sites_vcfs": [
                "gs://broad-cil-devel-bucket/input_data/7g8_gb4.combined.final.vcf.gz",
                "gs://broad-cil-devel-bucket/input_data/hb3_dd2.combined.final.vcf.gz",
                "gs://broad-cil-devel-bucket/input_data/3d7_hb3.combined.final.vcf.gz"
  ],
  "gatk_process_cohort.gvcf_fofn": "/cil/shed/sandboxes/amr/gatk_cloud/output/vcfs.out",
  "gatk_process_cohort.run_gatk_filter_genotypes": false,
  "gatk_process_cohort.snp_annotation": [
      "QD",
      "FS",
      "SOR",
      "DP",
      "MQ"
   ]
}
```
Note that the arrays of files for known sites are Google Storage URLs, but you may provide on-prem file paths as well
and these will be uploaded by Widdler upload. 


Now let's explain some of the other inputs inputs.
### gatk_process_cohort.run_gatk_vqsr

Including this in the json file with a value of true will run VQSR on the results of the joint genotype task. Ommit
or set to 'false' to disable VQSR.

If setting to false, the parameters that VQSR uses should be provided empty arrays ([]) or empty strings("") as 
appropriate. See the first example above.

### gatk_process_cohort.run_gatk_filter_genotypes

Including this in the json file with a value of true will run the filter genotypes script on the results of the
variant filtration task. Omit or set to 'false' to disable filter genotypes.

### gatk_process_cohort.run_snpeff

Including this in the json file with a value of true will run the snpeff task on the results filter genotypes if
that task was used, otherwise on the results of variant filtration. Omit or set to 'false' to snpeff.

Note that if using snpeff, then gatk_process_cohort.snpeff_db_name and gatk_process_cohort.snpeff_db_tgz must also be 
specified, as in the first json example above.

### gatk_process_cohort.reference_tgz

The process samples workflow, which would've been run prior to process_cohort, produces a reference TGZ file. Supply
the google storage path to the TGZ file produced by process samples here. If it was downloaded, you can also provide
the on-prem path for the reference tgz file, though this will lead to the file being uploaded to the Cloud again.

### gatk_process_cohort.output_disk_gb

This specifies the size in gigabytes for the output disk that will be assigned to store the outputs. If not 
specified, a default of 10GB is used. You may find you need to increase this when working with larger data sets.

### gatk_process_cohort.onprem_download_path

This optional parameter, when specified, will indicate where the outputs of the workflow should be downloaded. When
omitted, the outputs will not be downloaded but are still accessible in the Google bucket. The GCID bucket can be found
here:

https://console.cloud.google.com/storage/browser/4b66fc8a-tmp/

### gatk_process_cohort.handoff_files

If you want to only download some of the output files, you may specify that here by passing a json map structure in 
which the key is the task output and the value is an empty string, as seen in the example above. 

### gatk_process_cohort.gvcf_fofn

For the "gatk_process_cohort.gvcf_fofn" input, list all gvcf files like so:

```
/cil/shed/sandboxes/amr/gatk_cloud/Candida_Auris.gvcf
/cil/shed/sandboxes/amr/gatk_cloud/Candida_Auris2.gvcf
```

Note that a fofn using Google Storage URLS is possible and possibly preferable here since no upload would then be 
required for the files contained in the fofn. 

The process samples workflow produces this fofn and it can be passed directly to the process cohort workflow.

###


## USAGE WITH WIDDLER

Using Widdler, the Process Cohort Pipeline is executed in two phases:

### Widdler Upload

Data must first be uploaded to the cloud. The command format is as follows:

```sh widdler.sh upload <wdl file> <input_json>```

example:

```sh /cil/shed/apps/internal/widdler_dev/widdler.sh upload /cil/shed/apps/internal/wdl_cloud/taskdef.btl_gatk_process_cohort.wdl /cil/shed/sandboxes/amr/gatk_cloud/process_cohort.json```

Since we've set up defaults to be the settings for the GCID team, additional parameters do not need to be specified for
GCID people. Please note, however, that because this workflow uses subworkflows, all the subworkflows are provided for
you in the same directory as the master wdl. If you choose to use your own version of the WDL in a different location
you will also need to copy and unzip the "wdl_bundle.zip" file located in the same directory as the master WDL. The
presence of both the zip file and its unpacked contents in the same directory as the master WDL are prerequisites for
running the workflow with Widdler.

### Widdler Run

Once files are uploaded to the Google Cloud bucket, the workflow may be executed as follows:
```sh /cil/shed/apps/internal/widdler_dev/widdler.sh run /cil/shed/apps/internal/wdl_cloud/taskdef.btl_gatk_process_cohort.wdl </path/to/json/input> -S gscid-cloud -d /cil/shed/apps/internal/wdl_cloud/wdl_bundle_process_cohort.zip ```

You should see the usual widdler output providing you with the workflow ID of your workflow. For example:

```
-------------Cromwell Links-------------
   http://35.184.36.201:8000/api/workflows/v1/b460abc0-4377-4fbe-8ee8-3767b1084876/metadata
   http://35.184.36.201:8000/api/workflows/v1/b460abc0-4377-4fbe-8ee8-3767b1084876/timing
   {
       "status": "Submitted", 
       "id": "b460abc0-4377-4fbe-8ee8-3767b1084876"
   }
```

### Widdler Results

You can expect the same widdler e-mail that you would get with an on-prem run, except that your workflow root path
will be a clickable link that will take you directly to the execution directory for your gatk2cloud workflow.

```
Dear amr,

Your workflow with ID 4c873261-3c44-4e38-b47d-a259912b9ef7 has finished with a status of Succeeded. Attached are some files you might find useful.
Here's a brief summary:

Workflow Name: gatk_process_cohort
Started: 2018-05-18T18:41:28.589Z
Ended: 2018-05-18T18:42:30.291Z
Duration: 0 hours, 1 minutes, 1 seconds
workflowRoot: gs://4b66fc8a-tmp/cromwell-executions/gatk_process_cohort/4c873261-3c44-4e38-b47d-a259912b9ef7/ 
Timing graph: http://35.184.36.201:8000/api/workflows/v2/4c873261-3c44-4e38-b47d-a259912b9ef7/timing

Sincerely,
The Widdler
```

You will receive the same attachments as with on-prem workflows, but please note that the stderr and stdout files will not contain any information.
These will eventually be removed from the email for cloud workflows since the bucket URL is provided instead. 

### Getting Help

Feel free to use any of the following options:

1. Join the #gatk2cloud Slack channel.
2. E-mail amr@broadinstitute.org
3. Create a Taiga ticket.