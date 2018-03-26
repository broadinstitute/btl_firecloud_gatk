
workflow gatk_alignbam {
    String? onprem_download_path
    Map[String, String]? handoff_files

    File in_bam
    String sample_name
    File reference_tgz

    Boolean done = false
    String debug_dump_flag = "onfail"
    String preemptible = "0"

    if (size(in_bam,"GB")<2) {
        call gatk_alignbam_task {
            inputs: 
            unaligned_bam = in_bam,
            sample_name = sample_name,
            reference_tgz = reference_tgz,
            debug_dump_flag = debug_dump_flag,
            preemptible = preemptible,
        }
        done = true  
    }

    if ( !done) {
        call gatk_bamtofastq_task {
            inputs: 
            unaligned_bam = in_bam,
            sample_name = sample_name,
            debug_dump_flag = debug_dump_flag,
            preemptible = preemptible,
        }
        call gatk_bwamem_task {
            inputs: 
            fq1 = gatk_bamtofastq_task.fq1,
            fq2 = gatk_bamtofastq_task.fq2,
            sample_name = sample_name,
            reference_tgz = reference_tgz,
            debug_dump_flag = debug_dump_flag,
            preemptible = preemptible,
        }
        call gatk_sortsam_task {
            inputs: 
            aligned_bam = gatk_bwamem_task.aligned_bam,
            sample_name = sample_name,
            debug_dump_flag = debug_dump_flag,
            preemptible = preemptible,
        }
        call gatk_markduplicates_task {
            inputs: 
            sorted_bam = gatk_sortsam_task.sorted_bam,
            sample_name = sample_name,
            debug_dump_flag = debug_dump_flag,
            preemptible = preemptible,
        }
        call gatk_reordersam_task {
            inputs: 
            marked_bam = gatk_markduplicates_task.marked_bam,
            sample_name = sample_name,
            reference_tgz = reference_tgz,
            debug_dump_flag = debug_dump_flag,
            preemptible = preemptible,
        }
        done = true
    }

    output {
        File out_bam = select_first([gatk_alignbam_task.out_bam, gatk_reordersam_task.out_bam])
        File out_bam_index = select_first([gatk_alignbam_task.out_bam_index, gatk_reordersam_task.out_bam_index])
    }

}



task gatk_alignbam_task {
    String picard_path = "/cil/shed/apps/external/picard/current/bin/picard.jar"
    File unaligned_bam
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

    String output_disk_gb = "10"
    String boot_disk_gb = "10"
    String ram_gb = "13"
    String cpu_cores = "2"
    String preemptible 
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
#TMP_DIR=.
run('echo STARTING tar xvf to unpack reference')
run('date')
run('tar xvf ${reference_tgz}')

run('echo STARTING SamToFastq')
run('date')
run('java -Xms8G -jar ${picard_path} SamToFastq  INPUT=${unaligned_bam} FASTQ=${fq1_fn} SECOND_END_FASTQ=${fq2_fn} VALIDATION_STRINGENCY=LENIENT')

run('echo STARTING bwa mem')
run('date')
run('bwa mem -t 2 -R ${read_group} ref.fasta ${fq1_fn} ${fq2_fn} | samtools view -bS - > ${aligned_bam_fn}')

run('echo STARTING SortSam')
run('date')
run('java -Xms8G -jar ${picard_path} SortSam  I=${aligned_bam_fn} O=${sorted_bam_fn} SO=coordinate')

run('echo STARTING MarkDuplicates')
run('date')
run('java -Xms8G -jar ${picard_path} MarkDuplicates I=${sorted_bam_fn} O=${marked_bam_fn} M=${marked_duplicates_metrics_fn}')

run('echo STARTING ReorderSam')
run('date')
run('java -Xms8G -jar ${picard_path} ReorderSam   I=${marked_bam_fn} O=${out_bam_fn} R=ref.fasta')

run('echo STARTING index')
run('date')
run('samtools index ${out_bam_fn}')

run('echo DONE')
run('date')


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
####################### bam2fastq ######################


task gatk_bamtofastq_task {
    String picard_path = "/cil/shed/apps/external/picard/current/bin/picard.jar"
    File unaligned_bam
    String sample_name
    String fq1_fn = "${sample_name}.1.fq"
    String fq2_fn = "${sample_name}.2.fq"



    String output_disk_gb = "375"
    String boot_disk_gb = "10"
    String ram_gb = "208"
    String cpu_cores = "32"
    String preemptible 
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


run('echo STARTING SamToFastq')
run('date')
run('java -Xms198G -jar ${picard_path} SamToFastq  INPUT=${unaligned_bam} FASTQ=${fq1_fn} SECOND_END_FASTQ=${fq2_fn} VALIDATION_STRINGENCY=LENIENT')

run('echo DONE')
run('date')


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
        File fq1="${fq1_fn}"
        File fq2="${fq2_fn}"
        File monitor_start="monitor_start.log"
        File monitor_stop="monitor_stop.log"
        File dstat="dstat.log"
        File debug_bundle="debug_bundle.tar.gz"
    } runtime {
        docker : "gcr.io/btl-dockers/btl_gatk:1"
        memory: "${ram_gb}GB"
        cpu: "${cpu_cores}"
        disks: "local-disk ${output_disk_gb} LOCAL"
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

############### bwa mem ####################

task gatk_bwamem_task {
    String sample_name
    File fq1
    File fq2
    String aligned_bam_fn = "${sample_name}.aligned.bam"
    File reference_tgz

    String read_group = "\\'@RG\\\\\\\\tID:FLOWCELL_${sample_name}\\\\\\\\tSM:${sample_name}\\\\\\\\tPL:ILLUMINA\\\\\\\\tLB:LIB_${sample_name}\\'"

    String output_disk_gb = "375"
    String boot_disk_gb = "10"
    String ram_gb = "360"
    String cpu_cores = "96"
    String preemptible
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

run('echo STARTING bwa mem')
run('date')
run('bwa mem -t 96 -R ${read_group} ref.fasta ${fq1} ${fq2} | samtools view -bS - > ${aligned_bam_fn}')

run('echo DONE')
run('date')


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
        File aligned_bam = "${aligned_bam_fn}"
        File monitor_start="monitor_start.log"
        File monitor_stop="monitor_stop.log"
        File dstat="dstat.log"
        File debug_bundle="debug_bundle.tar.gz"
    } runtime {
        docker : "gcr.io/btl-dockers/btl_gatk:1"
        memory: "${ram_gb}GB"
        cpu: "${cpu_cores}"
        disks: "local-disk ${output_disk_gb} LOCAL"
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

################## sort sam #######################

task gatk_sortsam_task {
    String picard_path = "/cil/shed/apps/external/picard/current/bin/picard.jar"
    String sample_name
    File aligned_bam
    String sorted_bam_fn = "${sample_name}.sorted.bam"


    String output_disk_gb = "375"
    String boot_disk_gb = "10"
    String ram_gb = "60"
    String cpu_cores = "16"
    String preemptible 
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

run('echo STARTING SortSam')
run('date')
run('java -Xms50G -jar ${picard_path} SortSam  I=${aligned_bam} O=${sorted_bam_fn} SO=coordinate')


run('echo DONE')
run('date')


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
        File sorted_bam = "${sorted_bam_fn}"
        File monitor_start="monitor_start.log"
        File monitor_stop="monitor_stop.log"
        File dstat="dstat.log"
        File debug_bundle="debug_bundle.tar.gz"
    } runtime {
        docker : "gcr.io/btl-dockers/btl_gatk:1"
        memory: "${ram_gb}GB"
        cpu: "${cpu_cores}"
        disks: "local-disk ${output_disk_gb} LOCAL"
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


########### mark duplicates  ##############


task gatk_markduplicates_task {
    String picard_path = "/cil/shed/apps/external/picard/current/bin/picard.jar"
    String sample_name

    File sorted_bam
    String marked_bam_fn = "${sample_name}.marked_duplicates.bam"
    String marked_duplicates_metrics_fn = "${sample_name}.marked_duplicates.metrics"

    String output_disk_gb = "375"
    String boot_disk_gb = "10"
    String ram_gb = "208"
    String cpu_cores = "32"
    String preemptible 
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

run('echo STARTING MarkDuplicates')
run('date')
run('java -Xms198G -jar ${picard_path} MarkDuplicates I=${sorted_bam} O=${marked_bam_fn} M=${marked_duplicates_metrics_fn}')


run('echo DONE')
run('date')


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
        File marked_bam = "${marked_bam_fn}"
        File marked_duplicates_metrics = "${marked_duplicates_metrics_fn}"
        File monitor_start="monitor_start.log"
        File monitor_stop="monitor_stop.log"
        File dstat="dstat.log"
        File debug_bundle="debug_bundle.tar.gz"
    } runtime {
        docker : "gcr.io/btl-dockers/btl_gatk:1"
        memory: "${ram_gb}GB"
        cpu: "${cpu_cores}"
        disks: "local-disk ${output_disk_gb} LOCAL"
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


############## reorder sam  #####################

task gatk_reordersam_task {
    String picard_path = "/cil/shed/apps/external/picard/current/bin/picard.jar"
    String sample_name
    File marked_bam
    String out_bam_fn = "${sample_name}.bam"
    String out_bam_index_fn = "${out_bam_fn}.bai"
    File reference_tgz

    String output_disk_gb = "376"
    String boot_disk_gb = "10"
    String ram_gb = "60"
    String cpu_cores = "16"
    String preemptible 
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

run('echo STARTING ReorderSam')
run('date')
run('java -Xms50G -jar ${picard_path} ReorderSam   I=${marked_bam} O=${out_bam_fn} R=ref.fasta')

run('echo STARTING index')
run('date')
run('samtools index ${out_bam_fn}')

run('echo DONE')
run('date')


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
        disks: "local-disk ${output_disk_gb} LOCAL"
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