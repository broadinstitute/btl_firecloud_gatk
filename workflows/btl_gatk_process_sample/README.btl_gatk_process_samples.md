# README for Process Samples

Important! All examples below shown use widdler_dev. Please only use this version until the gatk2cloud changes have been pushed to production.

##TLDR

1. Prepare your json input file (see examples below).
2. Run widdler upload like so:
```sh /cil/shed/apps/internal/widdler_dev/widdler.sh upload /cil/shed/apps/internal/wdl_cloud/taskdef.btl_gatk_process_samples.wdl </path/to/json/input>```
3. Run widdler run like so:
```sh /cil/shed/apps/internal/widdler_dev/widdler.sh run /cil/shed/apps/internal/wdl_cloud/taskdef.btl_gatk_process_samples.wdl </path/to/json/input> -S gscid-cloud -d /cil/shed/apps/internal/wdl_cloud/wdl_bundle.zip```
4. Wait for widdler email.
5. If successful, find your downloaded files in the directory specified in the json input.

## Introduction

The Process Samples portion of the GATK2CLOUD pipeline allows users to perform the following
steps, in sequence:

1. Index the reference fasta file.
2. Optionally align the bam if a prealigned BAM is not provided.
3. Optionally perform TCIR (realigner target creator and indel realigner).
4. Optionally perform BQSR (base quality score recalibration).
5. Call haplotypes. 

## JSON Input File Format

The following is a template you can use to create the input json file for the workflow:
```
{
  "gatk_process_samples.known_sites_vcf_tbis": "Array[File]",
  "gatk_process_samples.use_tcir": "Boolean",
  "gatk_process_samples.gatk_haplotypecaller_task.extra_hc_params": "(optional) String?",
  "gatk_process_samples.gatk_haplotypecaller_task.ploidy": "(optional) String?",
  "gatk_process_samples.samples_tsv_fofn": "File",
  "gatk_process_samples.reference_fasta": "File",
  "gatk_process_samples.gatk_haplotypecaller_task.erc": "(optional) String?",
  "gatk_process_samples.use_bqsr": "Boolean",
  "gatk_process_samples.known_sites_vcfs": "Array[File]",
  "gatk_process_samples.prealigned": "(optional) Boolean",
  "gatk_process_samples.onprem_download_path": "(optional) String?"
}
```
Notes: 
* Optional items can be ommitted. 
* If not using TCIR/BQSR, write false in both of their respective fields.
* Due to a quirk with the WDL language, the two know_sites parameters can't be set to optional even though the tasks
that use them are being made optional. So if not using known_sites, just pass an empty array([]). 
* gatk_process_samples.gatk_haplotypecaller_task.ploidy defaults to '2' when not specified.

An example of a workflow input that is not using TCIR/BQSR/known_sites:
```
{
    "gatk_process_samples.reference_fasta":"/cil/shed/sandboxes/amr/gatk_cloud/minion_illumina_hybrid_clean_MT.fasta",
    "gatk_process_samples.known_sites_vcf_tbis":[],
    "gatk_process_samples.known_sites_vcfs":[],
    "gatk_process_samples.samples_tsv_fofn": "/cil/shed/sandboxes/amr/gatk_cloud/samples.tsv",
    "gatk_process_samples.use_tcir": false,
    "gatk_process_samples.use_bqsr": false,
    "gatk_process_samples.onprem_download_path": "/cil/shed/sandboxes/amr/gatk_cloud/output"
}
```

An example of a workflow that does use TCIR/BQSR/known_sites is as follows:
```
{
  "gatk_process_samples.samples_tsv_fofn": "/cil/shed/sandboxes/amr/gatk_cloud/plasmodium.fofn",
  "gatk_process_samples.use_bqsr": true,
  "gatk_process_samples.use_tcir": true,
  "gatk_process_samples.reference_fasta": "/cil/shed/sandboxes/amr/gatk_cloud/PlasmoDB-28_Pfalciparum3D7_Genome.fasta",
  "gatk_process_samples.onprem_download_path": "/cil/shed/sandboxes/amr/gatk_cloud/output",
  "gatk_process_samples.handoff_files": {
        "gatk_process_samples.out_gvcf":""
        },
        "gatk_process_samples.known_sites_vcf_tbis": [
                        "gs://broad-cil-devel-bucket/input_data/7g8_gb4.combined.final.vcf.gz.tbi",
                        "gs://broad-cil-devel-bucket/input_data/hb3_dd2.combined.final.vcf.gz.tbi",
                        "gs://broad-cil-devel-bucket/input_data/3d7_hb3.combined.final.vcf.gz.tbi"
        ],
        "gatk_process_samples.known_sites_vcfs": [
                        "gs://broad-cil-devel-bucket/input_data/7g8_gb4.combined.final.vcf.gz",
                        "gs://broad-cil-devel-bucket/input_data/hb3_dd2.combined.final.vcf.gz",
                        "gs://broad-cil-devel-bucket/input_data/3d7_hb3.combined.final.vcf.gz"
        ]
}
```
Note that the arrays of files for known sites are Google Storage URLs, but you may provide on-prem file paths as well
and these will be uploaded by Widdler upload. 


Now let's explain some of the other inputs inputs.

### gatk_process_samples.onprem_download_path

This optional parameter, when specified, will indicate where the outputs of the workflow should be downloaded. When
omitted, the outputs will not be downloaded but are still accessible in the Google bucket. The GCID bucket can be found
here:

https://console.cloud.google.com/storage/browser/4b66fc8a-tmp/

### gatk_process_samples.handoff_files

If you want to only download some of the output files, you may specify that here by passing a json map structure in 
which the key is the task output and the value is an empty string, as seen in the example above. 

### gatk_process_samples.samples_tsv_fofn

For the "gatk_process_samples.samples_tsv_fofn" input, list all sample names and the paths to their bam files
using tab-delimitted format like so:

```
Candida_Auris   /broad/hptmp/jlchang/wdl_cloud_test/samples/Candida_Auris.bam
Candida_Auris2  /broad/hptmp/jlchang/wdl_cloud_test/samples/Candida_Auris2.bam
```
###


## USAGE WITH WIDDLER

Using Widdler, the Process Samples Pipeline is executed in two phases:

### Widdler Upload

Data must first be uploaded to the cloud. The command format is as follows:

```sh widdler.sh upload <wdl file> <input_json>```

example:
```sh /cil/shed/apps/internal/widdler_dev/widdler.sh upload /cil/shed/apps/internal/wdl_cloud/taskdef.btl_gatk_process_samples.wdl /cil/shed/sandboxes/amr/gatk_cloud/process_samples.json```

Since we've set up defaults to be the settings for the GCID team, additional parameters do not need to be specified for
GCID people. Please note, however, that because this workflow uses subworkflows, all the subworkflows are provided for
you in the same directory as the master wdl. If you choose to use your own version of the WDL in a different location
you will also need to copy and unzip the "wdl_bundle.zip" file located in the same directory as the master WDL. The
presence of both the zip file and its unpacked contents in the same directory as the master WDL are prerequisites for
running the workflow with Widdler.

An example of the widdler output is as follows:

```gcloud auth activate-service-account --key-file /cil/shed/resources/widdler/gcid_service_account.json
   WARNING: Unable to create SSLContext.
   Activated service account credentials for: [gcid-cromwell-bucket-account@gcid-cromwell.iam.gserviceaccount.com]
   WARNING: Unable to create SSLContext.
   
   
   Updates are available for some Cloud SDK components.  To install them,
   please run:
     $ gcloud components update
   
   gcloud config set account gcid-cromwell-bucket-account@gcid-cromwell.iam.gserviceaccount.com
   Updated property [core/account].
   gsutil cp /cil/shed/sandboxes/amr/gatk_cloud/PlasmoDB-28_Pfalciparum3D7_Genome.fasta gs://4b66fc8a-tmp/broad-file-inputs/cil/shed/sandboxes/amr/gatk_cloud/PlasmoDB-28_Pfalciparum3D7_Genome.fasta
   Copying file:///cil/shed/sandboxes/amr/gatk_cloud/PlasmoDB-28_Pfalciparum3D7_Genome.fasta [Content-Type=application/octet-stream]...
   \ [1 files][ 22.6 MiB/ 22.6 MiB]                                                
   Operation completed over 1 objects/22.6 MiB.                                     
   gsutil cp /seq/plasmodium/data/pf3k/5.1/bam/PT0002-CW.bam gs://4b66fc8a-tmp/broad-file-inputs/seq/plasmodium/data/pf3k/5.1/bam/PT0002-CW.bam
   Copying file:///seq/plasmodium/data/pf3k/5.1/bam/PT0002-CW.bam [Content-Type=application/octet-stream]...
   ==> NOTE: You are uploading one or more large file(s), which would run          
   significantly faster if you enable parallel composite uploads. This
   feature can be enabled by editing the
   "parallel_composite_upload_threshold" value in your .boto
   configuration file. However, note that if you do this large files will
   be uploaded as `composite objects
   <https://cloud.google.com/storage/docs/composite-objects>`_,which
   means that any user who downloads such objects will need to have a
   compiled crcmod installed (see "gsutil help crcmod"). This is because
   without a compiled crcmod, computing checksums on composite objects is
   so slow that gsutil disables downloads of composite objects.
   
   | [1 files][  2.6 GiB/  2.6 GiB]   23.7 MiB/s                                   
   Operation completed over 1 objects/2.6 GiB.                                      
   gsutil cp /seq/plasmodium/data/pf3k/5.1/bam/PT0218-C.bam gs://4b66fc8a-tmp/broad-file-inputs/seq/plasmodium/data/pf3k/5.1/bam/PT0218-C.bam
   Copying file:///seq/plasmodium/data/pf3k/5.1/bam/PT0218-C.bam [Content-Type=application/octet-stream]...
   ==> NOTE: You are uploading one or more large file(s), which would run          
   significantly faster if you enable parallel composite uploads. This
   feature can be enabled by editing the
   "parallel_composite_upload_threshold" value in your .boto
   configuration file. However, note that if you do this large files will
   be uploaded as `composite objects
   <https://cloud.google.com/storage/docs/composite-objects>`_,which
   means that any user who downloads such objects will need to have a
   compiled crcmod installed (see "gsutil help crcmod"). This is because
   without a compiled crcmod, computing checksums on composite objects is
   so slow that gsutil disables downloads of composite objects.
   
   | [1 files][  4.1 GiB/  4.1 GiB]   23.1 MiB/s                                   
   Operation completed over 1 objects/4.1 GiB.                                      
   gsutil cp /cil/shed/sandboxes/amr/gatk_cloud/plasmodium.fofn.cloud gs://4b66fc8a-tmp/broad-file-inputs/cil/shed/sandboxes/amr/gatk_cloud/plasmodium.fofn.cloud
   Copying file:///cil/shed/sandboxes/amr/gatk_cloud/plasmodium.fofn.cloud [Content-Type=application/octet-stream]...
   - [1 files][  184.0 B/  184.0 B]                                                
   Operation completed over 1 objects/184.0 B.                                      
   The following files have been uploaded to 4b66fc8a-tmp:
   /cil/shed/sandboxes/amr/gatk_cloud/PlasmoDB-28_Pfalciparum3D7_Genome.fasta
   /seq/plasmodium/data/pf3k/5.1/bam/PT0002-CW.bam
   /seq/plasmodium/data/pf3k/5.1/bam/PT0218-C.bam
   /cil/shed/sandboxes/amr/gatk_cloud/plasmodium.fofn.cloud
```

### Widdler Run

Once files are uploaded to the Google Cloud bucket, the workflow may be executed as follows:
```sh /cil/shed/apps/internal/widdler_dev/widdler.sh run /cil/shed/apps/internal/wdl_cloud/taskdef.btl_gatk_process_samples.wdl </path/to/json/input> -S gscid-cloud -d /cil/shed/apps/internal/wdl_cloud/wdl_bundle.zip```

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

Your workflow with ID 3acbeacc-695a-429d-a914-1f85695fb2a9 has finished with a status of Succeeded. Attached are some files you might find useful.
Here's a brief summary:

Workflow Name: gatk_process_samples
Started: 2018-05-10T13:47:00.591Z
Ended: 2018-05-10T14:45:22.750Z
Duration: 0 hours, 58 minutes, 22 seconds
workflowRoot: gs://broad-cil-devel-bucket/gatk_process_samples/3acbeacc-695a-429d-a914-1f85695fb2a9/ 
Timing graph: http://35.193.85.62:8000/api/workflows/v2/3acbeacc-695a-429d-a914-1f85695fb2a9/timing

Sincerely,
The Widdler
```

You will receive the same attachments, but please note that the stderr and stdout files will not contain any information.
These will eventually be removed from the email for cloud workflows since the bucket URL is provided instead. 

### Getting Help

Feel free to use any of the following options:

1. Join the #gatk2cloud Slack channel.
2. E-mail amr@broadinstitute.org
3. Create a Taiga ticket.