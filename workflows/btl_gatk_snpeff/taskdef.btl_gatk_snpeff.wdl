
workflow gatk_snpeff {
    String? onprem_download_path
    Map[String, String]? handoff_files

    call gatk_snpeff_task
}

#TODO: How to handle the database issue: 1) Upload database to bucket location. 2) Pass location to WDL as a param
#TODO: Reconsider if we really need to pack/unpack the vcf file
#TODO: Reconsider if we really need to pass in 3 optional vcf parameters or we can just use one.
#TODO: Need to conditionally unpack tarball for database. Hard to do with CL python since it doesn't take expressions.
#TODO: Need to conditionally apply -dataDir flag. Hard to do with CL python since it doesn't take expressions.
#TODO: Would it make sense to upload the custom database to the compute VM itself and define in the config file there
# a local db store?

# Based on http://gatkforums.broadinstitute.org/gatk/discussion/50/adding-genomic-annotations-using-snpeff-and-variantannotator
task gatk_snpeff_task {
#    File ? vqsr_vcf_tb
#    File ? genotype_vcf_tb
#    File ? filtration_vcf_tb
#    File vcf_in_tb = select_first([filtration_vcf_tb, vqsr_vcf_tb, genotype_vcf_tb])
    File vcf_in
    String snpeff = "/cil/shed/apps/external/snpEff/snpEff-4.1g/snpEff.jar"
    String ? snpeff_db_tgz
    String ? snpeff_extra_params

    String cohort_name
    String vcf_out_fn = "${cohort_name}.snpeff.vcf"

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
import shutil
import subprocess

def run(cmd):
    print (cmd)
    subprocess.check_call(cmd, shell=True)

#run('echo STARTING tar xvf to unpack vcf')
#run('date')
#run('tar xvf ${vcf_in}')

run('tar xvf ${snpeff_db_tgz} -C data/')
run('sleep 9999999')
print "Copying Config..."
run('cp /cil/shed/apps/external/snpEff/snpEff-4.1g/snpEff.config .')

#shutil.copy2('/cil/shed/apps/external/snpEff/snpEff-4.1g/snpEff.config', '.')
print "Done copying config..."

print "Running snpeff..."
run(' \
    java -Xmx4G -jar ${snpeff} \
        -formatEff \
        -no-downstream \
        -no-intergenic \
        -no-upstream \
        -no-utr \
        ${"-dataDir data/ "} \
        ${vcf_in} \
        ${snpeff_extra_params} > ${vcf_out_fn} \
    ')

print "Done running snpeff."
"

        echo "$python_cmd"
        set +e
        python -c "$python_cmd"
        export _exit_code=$?
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
    runtime {
        docker : "gcr.io/btl-dockers/btl_gatk:1"
        memory: "${ram_gb}GB"
        cpu: "${cpu_cores}"
        disks: "local-disk ${output_disk_gb} HDD"
        bootDiskSizeGb: "${boot_disk_gb}"
        preemptible: "${preemptible}"
    }
    output {
        File vcf_out = vcf_out_fn
        File monitor_start="monitor_start.log"
        File monitor_stop="monitor_stop.log"
        File dstat="dstat.log"
        File debug_bundle="debug_bundle.tar.gz"
    }
    parameter_meta {
        vcf_in: "The input variants file."
        vcf_out: "The output variants file."
        snp_eff: "The path to the snpEff executable."
        snpeff_db: "The snpeff database to use."
        snpeff_extra_params: "A string for passing any additional parameters to snpeff."
    }
}
