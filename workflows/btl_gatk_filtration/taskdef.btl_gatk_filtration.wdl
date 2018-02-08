workflow gatk_filtration {
# http://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set



    File ? vqsr_vcf
    File ? genotype_vcf

    File sv_vcf = select_first([vqsr_vcf, genotype_vcf])


    String filtration_type #combined or variant

    if (filtration_type == "combined") {
        call gatk_combined_filtration_task
    }
    if (filtration_type == "variant") {
        call gatk_variant_filtration_task
    }
}




task gatk_combined_filtration_task {
    String gatk = "/humgen/gsa-hpprojects/GATK/bin/GenomeAnalysisTK-3.7-93-ge9d8068/GenomeAnalysisTK.jar"
    String sample_name
    File reference_tgz
    String filter_expression

    String vcf_out = "${sample_name}.ALL.filtered.vcf"

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

run('echo STARTING VariantFiltration')
run('date')
run('''\
        java -Xmx8G -jar ${gatk} \
            -T VariantFiltration \
            -R ref.fasta \
            -V ${sv_vcf} \
            --filterExpression '${filter_expression}' \
            --filterName my_variant_filter \
            -o ${vcf_out}
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
        File vcf_out = "${vcf_out}"
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

task gatk_variant_filtration_task {
    String gatk = "/humgen/gsa-hpprojects/GATK/bin/GenomeAnalysisTK-3.7-93-ge9d8068/GenomeAnalysisTK.jar"
    String sample_name
    File reference_tgz
    String snp_filter_expression
    String indel_filter_expression

    String vcf_out = "${sample_name}.variant.filtered.vcf"


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

run('echo STARTING SelectVariants-SNP')
run('date')
run('''\
        java -Xmx8G -jar ${gatk} \
            -T SelectVariants \
            -R ${ref} \
            -V ${vcf_in} \
            -selectType SNP \
            -o selectSNP.vcf
''')


run('echo STARTING VariantFiltration-SNP')
run('date')
run('''\
        java -Xmx8G -jar ${gatk} \
            -T VariantFiltration \
            -R ref.fasta \
            -V selectSNP.vcf \
            --filterExpression '${snp_filter_expression}' \
            --filterName my_variant_filter \
            -o filtered_SNPs.vcf
''')

run('echo STARTING SelectVariants-INDEL')
run('date')
run('''\
        java -Xmx8G -jar ${gatk} \
            -T SelectVariants \
            -R ${ref} \
            -V ${vcf_in} \
            -selectType INDEL \
            -o selectINDEL.vcf
''')


run('echo STARTING VariantFiltration-INDEL')
run('date')
run('''\
        java -Xmx8G -jar ${gatk} \
            -T VariantFiltration \
            -R ref.fasta \
            -V selectINDEL.vcf \
            --filterExpression '${indel_filter_expression}' \
            --filterName my_variant_filter \
            -o filtered_INDELs.vcf
''')


run('echo STARTING VariantFiltration-INDEL')
run('date')
run('''\
        java -jar -Xmx8G ${gatk} \
            -T CombineVariants \
            -R ref.fasta \
            --variant filtered_SNPs.vcf \
            --variant filtered_INDELs.vcf \
            -o ${vcf_out} \
            -genotypeMergeOptions UNIQUIFY
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
        File vcf_out = "${vcf_out}"
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



