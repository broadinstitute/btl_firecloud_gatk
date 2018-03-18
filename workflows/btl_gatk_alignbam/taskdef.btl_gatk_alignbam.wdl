
workflow gatk_alignbam {
    String? onprem_download_path
    Map[String, String]? handoff_files

    call gatk_alignbam_task 
}

# TODO break off the portion after BAM-to-fastq



task gatk_alignbam_task {
    String picard_path = "/cil/shed/apps/external/picard/current/bin/picard.jar"
    File in_bam
    String sample_name
    String fq1_fn = "${sample_name}.1.fq"
    String fq2_fn = "${sample_name}.2.fq"
    String aligned_bam_fn = "${sample_name}.aligned.bam"
    String sorted_bam_fn = "${sample_name}.sorted.bam"
    String marked_bam_fn = "${sample_name}.marked_duplicates.bam"
    String marked_duplicates_metrics_fn = "${sample_name}.marked_duplicates.metrics"
    String out_bam_fn = "${sample_name}.bam"
    String out_bam_index_fn = "${out_bam_fn}.bai"
    File reference_tgz

    String read_group = "\\'@RG\\\\\\\\tID:FLOWCELL_${sample_name}\\\\\\\\tSM:${sample_name}\\\\\\\\tPL:ILLUMINA\\\\\\\\tLB:LIB_${sample_name}\\'"

    String output_disk_gb 
    String boot_disk_gb = "10"
    String ram_gb = "60"
    String cpu_cores = "1"
    String preemptible = "0"
    String debug_dump_flag

    command {
        set -euo pipefail
        ln -sT `pwd` /opt/execution
        ln -sT `pwd`/../inputs /opt/inputs

        /opt/src/algutil/monitor_start.py
python_cmd="
import subprocess
def run(cmd):
    print (cmd)
    subprocess.check_call(cmd,shell=True)

run('echo STARTING tar xvf to unpack reference')
run('date')
run('tar xvf ${reference_tgz}')

run('echo STARTING SamToFastq')
run('date')
run('java -Xmx12G -jar ${picard_path} SamToFastq INPUT=${in_bam} FASTQ=${fq1_fn} SECOND_END_FASTQ=${fq2_fn} VALIDATION_STRINGENCY=LENIENT')

run('echo STARTING bwa mem')
run('date')
run('bwa mem -t 8 -R ${read_group} ref.fasta ${fq1_fn} ${fq2_fn} | samtools view -bS - > ${aligned_bam_fn}')

run('echo STARTING SortSam')
run('date')
run('java -Xmx8G -jar ${picard_path} SortSam I=${aligned_bam_fn} O=${sorted_bam_fn} SO=coordinate')

run('echo STARTING MarkDuplicates')
run('date')
run('java -Xmx8G -jar ${picard_path} MarkDuplicates I=${sorted_bam_fn} O=${marked_bam_fn} M=${marked_duplicates_metrics_fn}')

run('echo STARTING ReorderSam')
run('date')
run('java -Xmx8G -jar ${picard_path} ReorderSam I=${marked_bam_fn} O=${out_bam_fn} R=ref.fasta')

run('echo STARTING index')
run('date')
run('samtools index ${out_bam_fn}')

"

        echo "$python_cmd"
        set +e
        python -c "$python_cmd"
        export exit_code=$?
        set -e
        echo exit code is $exit_code
        ls

        # create bundle conditional on failure
        if [[ "${debug_dump_flag}" == "always" || ( "${debug_dump_flag}" == "onfail" && $exit_code -ne 0 ) ]]
        then
            echo "Creating debug bundle"
            # tar up the output directory
            touch debug_bundle.tar.gz
            tar cfz debug_bundle.tar.gz --exclude=debug_bundle.tar.gz .
        else
            touch debug_bundle.tar.gz
        fi     
        /opt/src/algutil/monitor_stop.py

        # exit statement must be the last line in the command block 
        exit $exit_code

    }
    output {
        File out_bam = "${out_bam_fn}"
        File out_bam_index = "${out_bam_index_fn}"
        File monitor_start="monitor_start.log"
        File monitor_stop="monitor_stop.log"
        File dstat="dstat.log"
        File debug_bundle="debug_bundle.tar.gz"
    } runtime {
        docker : "gcr.io/btl-dockers/btl_gatk:1"
        memory: "${ram_gb}GB"
        cpu: "${cpu_cores}"
        disks: "local-disk ${output_disk_gb} HDD"
        bootDiskSizeGb: "${boot_disk_gb}"
        preemptible: "${preemptible}"
    }
    parameter_meta {
        picard: "The absolute path to the picard jar to execute."
        in_bam: "The bam file to convert to fastq."
        sample_dir: "The sample-specific directory inside output_dir for each sample."
        sample_name: "The name of the sample as indicated by the 1st column of the gatk.samples_file json input."
        out_fq1: "The fastq file containing the first read of each pair."
        out_fq2: "The fastq file containing the second read of each pair"
    }
}

