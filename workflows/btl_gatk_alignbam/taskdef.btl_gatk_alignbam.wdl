
workflow gatk_alignbam {

    String picard = "/cil/shed/apps/external/picard/current/bin/picard.jar"
    File reference_tgz

    String sample_name

# TODO verify whether index needs to be passed in
    call SamToFastq {
        input:
        picard = picard,
        in_bam = in_bam,
        sample_name = sample_name,
    }

    call AlignBAM {
        input:
        reference_tgz = reference_tgz,
        fq1 = SamToFastq.fq1,
        fq2 = SamToFastq.fq2
    }
    call SortBAM {
        input:
        picard = picard,
        sample_name = sample_name,
        aligned_bam = AlignBAM.bam
    }

    call MarkDuplicates {
        input:
        picard = picard,
        sample_name = sample_name,
        sorted_bam = SortBAM.bam
    }

    # TODO why reorderbam rather than sortbam??
    call ReorderBAM {
        input:
        picard = picard,
        sample_name = sample_name,
        marked_bam = MarkDuplicates.bam,
        ref = index_reference_out,
        dict = index_reference_dict
    }




}


task SamToFastq {
    String picard
    File in_bam
    String sample_name
    String out_fq1 = "${sample_name}.1.fq"
    String out_fq2 = "${sample_name}.2.fq"

    String output_disk_gb 
    String boot_disk_gb = "10"
    String ram_gb = "3"
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

run('java -Xmx12G -jar ${picard} SamToFastq INPUT=${in_bam} FASTQ=${out_fq1} SECOND_END_FASTQ=${out_fq2} VALIDATION_STRINGENCY=LENIENT')
"

        echo "$python_cmd"
        python -c "$python_cmd"
        export exit_code=$?
        echo exit code is $exit_code

        # create bundle conditional on failure
        if [[ "$debug_dump_flag" == "always" || ( "$debug_dump_flag" == "onfail" && $exit_code -ne 0 ) ]]
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
        String done = "Done"
        File fq1 = out_fq1
        File fq2 = out_fq2
        File monitor_start="monitor_start.log"
        File monitor_stop="monitor_stop.log"
        File dstat="dstat.log"
    } runtime {
        task_name: "SamToFastq"
        docker : "gcr.io/btl-dockers/btl_gatk:1"
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

task AlignBAM {
    String sample_name
    File reference_tgz
    File fq1
    File fq2
    String read_group = "'@RG\\tID:FLOWCELL_${sample_name}\\tSM:${sample_name}\\tPL:ILLUMINA\\tLB:LIB_${sample_name}'"

    String output_disk_gb 
    String boot_disk_gb = "10"
    String ram_gb = "3"
    String cpu_cores = "1"
    String preemptible = "0"
    String debug_dump_flag

    command {
        set -euo pipefail
        ln -sT `pwd` /opt/execution
        ln -sT `pwd`/../inputs /opt/inputs

        /opt/src/algutil/monitor_start.py
        tar xvf ${reference_tgz}

        bwa mem -t 8 -R ${read_group} ref.fasta ${fq1} ${fq2} | samtools view -bS -> ${sample_name}.aligned.bam
        /opt/src/algutil/monitor_stop.py
    }
    output {
        File bam = "${sample_name}.aligned.bam"
        String done = "done"
        File monitor_start="monitor_start.log"
        File monitor_stop="monitor_stop.log"
        File dstat="dstat.log"
    } runtime {
        task_name: "AlignBAM"
        docker : "gcr.io/btl-dockers/btl_gatk:1"
    }
    parameter_meta {
        ref: "fasta file of reference genome"
        sample_dir: "The sample-specific directory inside output_dir for each sample."
        sample_name: "The name of the sample as indicated by the 1st column of the gatk.samples_file json input."
        fq_array: "An array containing the paths to the first and second fastq files."
        read_group: "The read group string that will be included in the bam header."
    }
}

task SortBAM {
    String picard
    String sample_name
    File aligned_bam

    String output_disk_gb 
    String boot_disk_gb = "10"
    String ram_gb = "3"
    String cpu_cores = "1"
    String preemptible = "0"
    String debug_dump_flag

    command {
        set -euo pipefail
        ln -sT `pwd` /opt/execution
        ln -sT `pwd`/../inputs /opt/inputs

        /opt/src/algutil/monitor_start.py
        java -Xmx8G -jar ${picard} SortSam I=${aligned_bam} O=${sample_name}.sorted.bam SO=coordinate
        /opt/src/algutil/monitor_stop.py

    }
    output {
        File bam = "${sample_name}.sorted.bam"
        String done = "Done"
        File monitor_start="monitor_start.log"
        File monitor_stop="monitor_stop.log"
        File dstat="dstat.log"
    } runtime {
        task_name: "SortBAM"
        docker : "gcr.io/btl-dockers/btl_gatk:1"
    }
    parameter_meta {
        picard: "The absolute path to the picard jar to execute."
        aligned_sam: "The sam file to be sorted."
        sample_dir: "The sample-specific directory inside output_dir for each sample."
        bam: "The sorted bam file."
    }
}

task MarkDuplicates {
    String picard
    String sample_name
    File sorted_bam

    String output_disk_gb 
    String boot_disk_gb = "10"
    String ram_gb = "3"
    String cpu_cores = "1"
    String preemptible = "0"
    String debug_dump_flag


    command {
        set -euo pipefail
        ln -sT `pwd` /opt/execution
        ln -sT `pwd`/../inputs /opt/inputs

        /opt/src/algutil/monitor_start.py
        java -Xmx8G -jar ${picard} MarkDuplicates I=${sorted_bam} O=${sample_name}.marked_duplicates.bam M=${sample_name}.marked_duplicates.metrics
        /opt/src/algutil/monitor_stop.py
    }
    output {
        File bam = "${sample_name}.marked_duplicates.bam"
        String done = "done"
        File monitor_start="monitor_start.log"
        File monitor_stop="monitor_stop.log"
        File dstat="dstat.log"
    } runtime {
        task_name: "MarkDuplicates"
        docker : "gcr.io/btl-dockers/btl_gatk:1"
    }
    parameter_meta {
        picard: "The absolute path to the picard jar to execute."
        sorted_bam: "The sorted bam file to mark duplicates on."
        out_bam: "The bam file where duplicates are marked."
        metrics: "The marked duplicates metrics file."
    }
}

task ReorderBAM {
    String picard
    String sample_name
    File marked_bam
    String out_bam = "${sample_name}.reordered.bam"
    File ref
    File dict

    String output_disk_gb 
    String boot_disk_gb = "10"
    String ram_gb = "3"
    String cpu_cores = "1"
    String preemptible = "0"
    String debug_dump_flag


    command {
        set -euo pipefail
        ln -sT `pwd` /opt/execution
        ln -sT `pwd`/../inputs /opt/inputs

        /opt/src/algutil/monitor_start.py
        java -Xmx8G -jar ${picard} ReorderSam I=${marked_bam} O=${out_bam} R=${ref}
        samtools index ${out_bam}
        /opt/src/algutil/monitor_stop.py

    }
    output {
        File bam = out_bam
        File index = "${sample_name}.reordered.bam.bai"
        String done = "done"
        File monitor_start="monitor_start.log"
        File monitor_stop="monitor_stop.log"
        File dstat="dstat.log"
    } runtime {
        task_name: "ReorderBAM"
        docker : "gcr.io/btl-dockers/btl_gatk:1"
    }
    parameter_meta {
        picard: "The absolute path to the picard jar to execute."
        marked_bam: "The bam file where duplicates are marked."
        out_bam: "The reordered bam file."
        ref: "fasta file of reference genome"
    }
}
