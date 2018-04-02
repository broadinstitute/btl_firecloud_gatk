import pytest
import tarfile
import glob
import subprocess

def test_verify_valid_comparison_dir(comparison_dir):

    def bam_read_count(fname):
        p = subprocess.Popen(['samtools view -F 0x40 ' + fname + ' | cut -f1 | sort | uniq | wc -l'], shell=True,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        lines = ""
        for line in p.stdout:
            lines = lines + line

        return int(lines.strip().split()[0])

    GATK_GOAL_PATH = "/cil/shed/resources/wdl/gatk/bqsr/output/CandidaAuris.bqsr.bam"
    COMPARISON_PATH= comparison_dir + "/CandidaAuris.bqsr.bam"

    assert(bam_read_count(GATK_GOAL_PATH) == bam_read_count(COMPARISON_PATH))