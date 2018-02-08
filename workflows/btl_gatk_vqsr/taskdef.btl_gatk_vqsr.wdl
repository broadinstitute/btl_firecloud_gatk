workflow gatk_vqsr{
    # https://software.broadinstitute.org/gatk/documentation/article.php?id=1259
    call gatk_vqsr_task
}


task gatk_vqsr_task {
    String gatk = "/humgen/gsa-hpprojects/GATK/bin/GenomeAnalysisTK-3.7-93-ge9d8068/GenomeAnalysisTK.jar"
    File reference_tgz
    File ? intervals

    File genotype_caller_vcf
    Array[String] snp_resource
    Array[String] indel_resource
    Array[String] snp_annotation
    Array[String] indel_annotation
    Int ? snp_max_gaussians
    Int ? indel_max_gaussians
    Int ? mq_cap_snp
    Int ? mq_cap_indel
    Float ts_filter_snp
    Float ts_filter_indel
    String ? extra_vr_params



    String output_disk_gb 
    String boot_disk_gb = "10"
    String ram_gb = "10"
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

run('echo STARTING VariantRecalibrator-SNP')
run('date')
run('''\
		java -Xmx8G -jar ${gatk} \
			-T VariantRecalibrator \
			-R ref.fasta \
			${default="" "--intervals " + intervals} \
			-input ${genotype_caller_vcf} \
			-mode "snp" \
            -resource:${sep=" -resource:" snp_resource} \
			-recalFile snp.recal \
			-tranchesFile snp.tranches \
			-rscriptFile snp.plots.R \
			-an ${sep=" -an " snp_annotation} \
			--maxGaussians ${snp_max_gaussians} \
			--MQCapForLogitJitterTransform ${mq_cap_snp}
			${default="\n" extra_vr_params}
''')

run('echo STARTING ApplyRecalibration-SNP')
run('date')
run('''\
		java -Xmx8G -jar ${gatk} \
            -T ApplyRecalibration \
            -R ref.fasta \
            -input ${genotype_caller_vcf} \
            --ts_filter_level ${ts_filter_snp} \
            -tranchesFile snp.tranches \
            -recalFile snp.recal \
            -mode "snp" \
            -o snp.recalibrated.filtered.vcf
''')

run('echo STARTING VariantRecalibrator-INDEL')
run('date')
run('''\
		java -Xmx8G -jar ${gatk} \
			-T VariantRecalibrator \
			-R ref.fasta \
			${default="" "--intervals " + intervals} \
			-input snp.recalibrated.filtered.vcf \
			-mode "indel" \
            -resource:${sep=" -resource:" indel_resource} \
			-recalFile indel.recal \
			-tranchesFile indel.tranches \
			-rscriptFile indel.plots.R \
			-an ${sep=" -an " indel_annotation} \
			--maxGaussians ${indel_max_gaussians} \
			--MQCapForLogitJitterTransform ${mq_cap_indel}
			${default="\n" extra_vr_params}
''')

run('echo STARTING ApplyRecalibration-INDEL')
run('date')
run('''\
		java -Xmx8G -jar ${gatk} \
            -T ApplyRecalibration \
            -R ref.fasta \
            -input snp.recalibrated.filtered.vcf \
            --ts_filter_level ${ts_filter_indel} \
            -tranchesFile indel.tranches \
            -recalFile indel.recal \
            -mode "indel" \
            -o snp_indel.recalibrated.filtered.vcf
''')

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

        # create bundle conditional on failure of the Python section
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
        File out_bam = "${out_bam}"
        File out_bam_index = "${out_bam_index}"
        File recalibration_plots = "${recalibration_plots_fn}"
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

    }

}

