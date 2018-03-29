import pytest
import tarfile
import glob
import subprocess

def test_verify_valid_comparison_dir(comparison_dir):

    def bam_read_count(fname):
        p = subprocess.Popen(['samtools', 'view', '-F', '0x40', fname, '|', 'cut', '-f1', '|', 'sort', '|', 'uniq',
                              '|', 'wc', '-l'], stdout=subprocess.PIPE,
                                              stderr=subprocess.PIPE)
        result, err = p.communicate()

        if p.returncode != 0:
            raise IOError(err)

        return int(result.strip().split()[0])

    GATK_GOAL_PATH = "/cil/shed/resources/wdl/gatk/snpeff/output/CandidaAurisCohort.snpeff.vcf"
    COMPARISON_PATH= comparison_dir + "/CandidaAurisCohort.snpeff.vcf"

    tar_gz = glob.glob(comparison_dir + '/*.tgz')[0]
    tar = tarfile.open(tar_gz)
    tar.extractall(path=comparison_dir)
    tar.close()

    assert(bam_read_count(GATK_GOAL_PATH) == bam_read_count(COMPARISON_PATH))