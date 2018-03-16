import pytest
import glob
import os
import gffutils
import filecmp

def test_verify_valid_comparison_dir(comparison_dir):
    assert(len(comparison_dir) > 0)