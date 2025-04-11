import tempfile
import shutil
import pytest
import os
import matplotlib.pyplot as plt
import pandas as pd


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

@pytest.fixture
def no_show(monkeypatch):
    monkeypatch.delattr(plt, 'show')

@pytest.fixture
def sample_data():
    data = {
        'Seq_1': ['A', 'A', 'B'],
        'Seq_2': ['B', 'C', 'C'],
        'Per_Id': [90, 85, 80],
        'Per_Sim': [0.9, 0.85, 0.8]
    }
    return pd.DataFrame(data)

@pytest.fixture
def protein_data():
    data = {'Seq_1': ['A', 'B', 'A'],
            'Seq_2': ['B', 'C', 'C'],
            'Per_Id': [0.8, 0.9, 0.7],
            'Per_Sim': [0.6, 0.5, 0.4]}
    df = pd.DataFrame(data)
    return df

