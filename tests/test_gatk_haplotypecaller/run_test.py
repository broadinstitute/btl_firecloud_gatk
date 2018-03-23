import pytest
import tarfile
import glob
import subprocess

def test_verify_valid_comparison_dir(comparison_dir):

    def word_count(fname):
        p = subprocess.Popen(['wc', '-l', fname], stdout=subprocess.PIPE,
                                              stderr=subprocess.PIPE)
        result, err = p.communicate()

        if p.returncode != 0:
            raise IOError(err)

        return int(result.strip().split()[0])

    GATK_GOAL_PATH = "/cil/shed/resources/wdl/gatk/haplotypecaller/output/CandidaAuris.gvcf"
    COMPARISON_PATH= comparison_dir + "/CandidaAuris.gvcf"

    assert(word_count(GATK_GOAL_PATH) == word_count(COMPARISON_PATH))