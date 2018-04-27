from runcromwellremote import run_wdl

wdl = 'taskdef.btl_gatk_indexref.wdl'

inputs = {
    'gatk_indexref.gatk_indexref_task.ref_fasta':'gs://4b66fc8a/test_data/minion_illumina_hybrid_clean_MT.fasta',
    'gatk_indexref.gatk_indexref_task.ref_name':"minion_illumina_hybrid_clean_MT",
    'gatk_indexref.gatk_indexref_task.output_disk_gb':'13',
    'gatk_indexref.gatk_indexref_task.debug_dump_flag':'onfail',

    }

(status,outputs) = run_wdl(wdl,inputs)

print(outputs)