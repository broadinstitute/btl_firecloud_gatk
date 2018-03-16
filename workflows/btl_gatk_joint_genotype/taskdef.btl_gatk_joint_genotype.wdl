workflow gatk_joint_genotype {
    String? onprem_download_path
    Map[String, String]? handoff_files

    call gatk_joint_genotype_task
}


task gatk_joint_genotype_task {

    String ? extra_gg_params
    Array[File] HaplotypeCaller_gvcfs
    Boolean ? all_sites
    String gatk_path = "/humgen/gsa-hpprojects/GATK/bin/GenomeAnalysisTK-3.7-93-ge9d8068/GenomeAnalysisTK.jar"


    File reference_tgz

    String cohort_name
    String vcf_out_fn = "${cohort_name}.vcf"




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

# add intervals back in when actually scattering haplotype caller

run('''\
        java -Xmx8G -jar ${gatk_path} \
            -T GenotypeGVCFs \
            -R ref.fasta \
            -o ${vcf_out_fn} \
            -V ${sep=" -V " HaplotypeCaller_gvcfs} \
            ${true="-allSites" false="" all_sites} \
            ${default="\n" extra_gg_params}
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
        File vcf_out = "${vcf_out_fn}"
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






