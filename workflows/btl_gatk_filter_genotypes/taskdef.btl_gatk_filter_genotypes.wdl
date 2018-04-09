
workflow gatk_filter_genotypes {
    String? onprem_download_path
    Map[String, String]? handoff_files

    call gatk_filter_genotypes_task
}


task gatk_filter_genotypes_task {
    File vcf_in
    String cohort_name
    String vcf_out_fn = "${cohort_name}.genotype.filtered.vcf"

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

run('echo running filterGatkGenotypes.py...')
run('''python /opt/src/filterGatkGenotypes.py \
    --min_GQ 50 \
    --min_percent_alt_in_AD 0.8 \
    --min_total_DP 10 \
    ${vcf_in} > ${vcf_out_fn}
''')
"
    echo "$python_cmd"
    python -c "$python_cmd"
    export exit_code=$?
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
}


