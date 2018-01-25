workflow gatk_tcir {

    File index_reference_out
    File index_reference_dict
    File index_reference_amb
    File index_reference_ann
    File index_reference_bwt
    File index_reference_fai
    File index_reference_pac
    File index_reference_sa
    String sample_name
    String gatk = "/humgen/gsa-hpprojects/GATK/bin/GenomeAnalysisTK-3.7-93-ge9d8068/GenomeAnalysisTK.jar"
    File bam_input
    File bam_input_index


    call RealignerTargetCreator {
        input:
        gatk = gatk,
        ref = index_reference_out,
        dict = index_reference_dict,
        amb = index_reference_amb,
        ann = index_reference_ann,
        bwt = index_reference_bwt,
        fai = index_reference_fai,
        pac = index_reference_pac,
        sa = index_reference_sa,
        sample_name = sample_name,
        bam_input = bam_input,
        bam_input_index = bam_input_index
    }
    call IndelRealigner {
        input:
        gatk = gatk,
        ref = index_reference_out,
        dict = index_reference_dict,
        amb = index_reference_amb,
        ann = index_reference_ann,
        bwt = index_reference_bwt,
        fai = index_reference_fai,
        pac = index_reference_pac,
        sa = index_reference_sa,
        sample_name = sample_name,
        bam_input = bam_input,
        bam_input_index = bam_input_index,
        intervals = RealignerTargetCreator.intervals
    }
}


task RealignerTargetCreator {
    String gatk
    File ref
    File dict
    File amb
    File ann
    File bwt
    File fai
    File pac
    File sa
    File bam_input
    File bam_input_index
    String sample_name
    String out = "${sample_name}.interval_list"
    command {
        set -euo pipefail
        /opt/src/algutil/monitor_start.py
        java -Xmx8G -jar ${gatk} -T RealignerTargetCreator -nct 1 -nt 24 -R ${ref} -I ${bam_input} -o ${out}
        /opt/src/algutil/monitor_stop.py

    }
    output {
        File intervals = out
        File monitor_start="monitor_start.log"
        File monitor_stop="monitor_stop.log"
        File dstat="dstat.log"
    } runtime {
        task_name: "RealignerTargetCreator"
        docker : "gcr.io/btl-dockers/btl_gatk:1"
    }
    parameter_meta {
        gatk: "The absolute path to the gatk executable jar."
        ref: "fasta file of reference genome"
        in_bam: "The input bam for the gatk task"
        out: "The intervals list to be used by IndelRealigner"
    }
}

task IndelRealigner {
    String gatk
    File ref
    File dict
    File amb
    File ann
    File bwt
    File fai
    File pac
    File sa
    File bam_input
    File bam_input_index
    File intervals
    String sample_name
    String out = "${sample_name}.indels_realigned.bam"

    command {
        set -euo pipefail
        /opt/src/algutil/monitor_start.py
        java -Xmx4G -jar ${gatk} -T IndelRealigner -nct 1 -nt 1 -R ${ref} -I ${bam_input} -targetIntervals ${intervals} -o ${out}
        samtools index ${out}
        /opt/src/algutil/monitor_stop.py

    }
    output {
        File bam = "${sample_name}.indels_realigned.bam"
        File index = "${sample_name}.indels_realigned.bam.bai"
        String done = "done"
        File monitor_start="monitor_start.log"
        File monitor_stop="monitor_stop.log"
        File dstat="dstat.log"
    } runtime {
        task_name: "IndelRealigner"
        docker : "gcr.io/btl-dockers/btl_gatk:1"
    }
    parameter_meta {
        gatk: "The absolute path to the gatk executable jar."
        ref: "fasta file of reference genome"
        in_bam: "The input bam for the gatk task"
        intervals: "The intervals list to be used by IndelRealigner"
        out: "the bam including realigned indels."
    }
}