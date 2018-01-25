
workflow gatk_alignbam {

    String picard = "/cil/shed/apps/external/picard/current/bin/picard.jar"
    File index_reference_out
    File index_reference_dict
    File index_reference_amb
    File index_reference_ann
    File index_reference_bwt
    File index_reference_fai
    File index_reference_pac
    File index_reference_sa
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
        ref = index_reference_out,
        dict = index_reference_dict,
        amb = index_reference_amb,
        ann = index_reference_ann,
        bwt = index_reference_bwt,
        fai = index_reference_fai,
        pac = index_reference_pac,
        sa = index_reference_sa,
        sample_name = sample_name,
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
    command {
        set -euo pipefail
        /opt/src/algutil/monitor_start.py
        java -Xmx12G -jar ${picard} SamToFastq INPUT=${in_bam} FASTQ=${out_fq1} SECOND_END_FASTQ=${out_fq2} VALIDATION_STRINGENCY=LENIENT
        /opt/src/algutil/monitor_stop.py

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
    File ref
    File dict
    File amb
    File ann
    File bwt
    File fai
    File pac
    File sa
    File fq1
    File fq2
    String read_group = "'@RG\\tID:FLOWCELL_${sample_name}\\tSM:${sample_name}\\tPL:ILLUMINA\\tLB:LIB_${sample_name}'"
    command {
        set -euo pipefail
        /opt/src/algutil/monitor_start.py
        bwa mem -t 8 -R ${read_group} ${ref} ${fq1} ${fq2} | samtools view -bS -> ${sample_name}.aligned.bam
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
    command {
        set -euo pipefail
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
    command {
        set -euo pipefail
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
    command {
        set -euo pipefail
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
