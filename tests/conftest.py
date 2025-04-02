import tempfile
import shutil

import pytest
import os

@pytest.fixture
def temp_dir():
    tmpd = tempfile.mktemp()
    yield tmpd
    shutil.rmtree(tmpd)

@pytest.fixture
def temp_file(temp_dir):
    fd, fn = tempfile.mkstemp(dir=temp_dir)
    try:
        with open(fd, 'w') as fh:
            yield fn
    finally:
        if os.path.exists(fn):
            os.remove(fn)

@pytest.fixture
def inputs_dna(tmp_path):
    input_file = 'test_data/test_for_python_DNA_aligned.fasta'
    output_file = tmp_path / "test_output_dna.tsv"

    yield str(input_file), str(output_file)

@pytest.fixture
def inputs_protein(tmp_path):
    input_file = 'test_data/test_for_python_protein_aligned.fasta'
    output_file = tmp_path / "test_output_protein.tsv"
    similarity_groups = "GA,ST,NDEQ,KRH,VILM,YWF,CP"
    my_similarity_groups = [tuple(group.strip()) for group in similarity_groups.split(',')]

    yield str(input_file), str(output_file), my_similarity_groups