 workflow gatk_bqsr {
 
    String gatk = "/humgen/gsa-hpprojects/GATK/bin/GenomeAnalysisTK-3.7-93-ge9d8068/GenomeAnalysisTK.jar"

    File index_reference_out
    File index_reference_dict
    File index_reference_amb
    File index_reference_ann
    File index_reference_bwt
    File index_reference_fai
    File index_reference_pac
    File index_reference_sa
    String sample_name
    File ? indelrealigner_bam
    File ? uncleaned_bam
    File ? indelrealigner_bam_index
    File ? uncleaned_bam_index

    
    File bqsr_bam = select_first([indelrealigner_bam, uncleaned_bam])
    File bqsr_index = select_first([indelrealigner_bam_index, uncleaned_bam_index])
    #File known_sites

    call BaseRecalibrator as BaseRecal_1 {
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
        bam = bqsr_bam,
        index = bqsr_index,
        #known_sites = known_sites,
        out_file = "recal_data.table"
    }

    call BaseRecalibrator as BaseRecal_2 {
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
        bam = bqsr_bam,
        index = bqsr_index,
        #known_sites = known_sites,
        bqsr = BaseRecal_1.out,
        out_file = "post_recal_data.table"
    }

    call AnalyzeCovariates {
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
        before = BaseRecal_1.out,
        after = BaseRecal_2.out
    }


    call PrintReads {
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
        in_bam = bqsr_bam,
        bqsr = BaseRecal_2.out
    }


}

        #TODO reincorporate known_sites arg
        # Array[File] known_sites
        # -knownSites ${sep=" -knownSites " known_sites}

task BaseRecalibrator {
    String gatk
    File ref
    File dict
    File amb
    File ann
    File bwt
    File fai
    File pac
    File sa
    File bam
    File index
    String ? bqsr
    String out_file
    command {
        set -euo pipefail
        /opt/src/algutil/monitor_start.py
        samtools index ${bam}
        java -Xmx4G  -jar ${gatk} -T BaseRecalibrator -nct 8 -nt 1 -R ${ref} -I  ${bam} -o ${out_file} ${" -BQSR " + bqsr}
        /opt/src/algutil/monitor_stop.py
    }
    output {
        File out = "${out_file}"
        File monitor_start="monitor_start.log"
        File monitor_stop="monitor_stop.log"
        File dstat="dstat.log"
    } runtime {
        task_name: "BaseRecalibrator"
        docker : "gcr.io/btl-dockers/btl_gatk:1"
    }
    parameter_meta {
        gatk: "The absolute path to the gatk executable jar."
        ref: "fasta file of reference genome"
        bam_list: "The text file containing absolute paths to all bam files, one per line."
        known_sites: "An array of databases of known polymorphic sites."
        bqsr: "Full path to bqsr file."
    }
}

task AnalyzeCovariates {
    String gatk
    File ref
    File dict
    File amb
    File ann
    File bwt
    File fai
    File pac
    File sa
    String sample_name
    String before
    String after
    command {
        set -euo pipefail
        /opt/src/algutil/monitor_start.py

        #use R-3.1
        java -Xmx8G -jar ${gatk} \
        -T AnalyzeCovariates \
        -R ${ref} \
        -before ${before} \
        -after ${after} \
        -plots ${sample_name}.recalibration_plots.pdf
        /opt/src/algutil/monitor_start.py

    } output {
        String done = "done"
        File out = "${sample_name}.recalibration_plots.pdf"
        File monitor_start="monitor_start.log"
        File monitor_stop="monitor_stop.log"
        File dstat="dstat.log"
    } runtime {
        task_name: "AnalyzeCovariates"
        docker : "gcr.io/btl-dockers/btl_gatk:1"
    }
    parameter_meta {
        gatk: "The absolute path to the gatk executable jar."
        ref: "fasta file of reference genome"
        sample_dir: "The sample-specific directory inside output_dir for each sample."
        sample_name: "The name of the sample as indicated by the 1st column of the gatk.samples_file json input."
        known_sites: "An array of databases of known polymorphic sites."
        before: "The table output from the first BaseRecalibrator step."
        after: "The table output from the second BaseRecalibrator step."
    }
}

task PrintReads {
    String gatk
    String sample_name
    String ref
    File in_bam
    File dict
    File amb
    File ann
    File bwt
    File fai
    File pac
    File sa
    String bqsr
    String out =  "${sample_name}.bqsr.bam"
    command {
        set -euo pipefail
        /opt/src/algutil/monitor_start.py
        samtools index ${in_bam}
        java -Xmx4G -jar ${gatk} -T PrintReads -nct 8 -nt 1 -R ${ref} -I ${in_bam} -BQSR ${bqsr} -o ${out}
        samtools index ${out}
        /opt/src/algutil/monitor_stop.py

    }
    output {
        File bam = out
        # File index = "$[sample_name}.bqsr.bam.bai"
        String done = "done"
        File monitor_start="monitor_start.log"
        File monitor_stop="monitor_stop.log"
        File dstat="dstat.log"
    } runtime {
        task_name: "PrintReads"
        docker : "gcr.io/btl-dockers/btl_gatk:1"
    }
    parameter_meta {
        gatk: "The absolute path to the gatk executable jar."
        ref: "fasta file of reference genome"
        in_bam: "The input bam for PrintReads."
        bqsr: "Full path to bqsr file."
        out: "The bam file with bqsr applied to it."
    }
}