import aps_dop.APS_DoP as prog
import pandas as pd
import pytest
import sys
import os
import glob
from unittest.mock import patch
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import subprocess

def test_run_aps_dop(inputs_dna):
    """
    Tests the execution of the APS_DoP.py script as a subprocess for DNA alignment.
    It checks if the script runs successfully and returns a zero exit code.
    """
    test_dir = os.path.dirname(os.path.abspath(__file__))
    aps_dop_path = os.path.join(test_dir, "..", "aps_dop", "APS_DoP.py")
    aps_dop_path = os.path.abspath(aps_dop_path)
    input_file, output_file = inputs_dna

    try:
        result = subprocess.run(
            ["python", aps_dop_path, "-i", input_file, "-o", output_file, "-m", "DNA"],
            capture_output=True,
            text=True,
            check=True,
        )
        assert result.returncode == 0
    except subprocess.CalledProcessError as e:
        pytest.fail(f"Subprocess failed with error: {e.stderr}")

def test_dna(inputs_dna):
    """
    Tests the DNA alignment processing function.
    It checks if the output file is created with the expected content and format.
    """
    input_file, output_file = inputs_dna
    prog.process_alignment_DNA(input_file, output_file)

    with open(output_file, 'r') as f:
        lines = f.readlines()

    assert len(lines) == 4
    assert "Seq_1\tSeq_2\tPer_Id\tLength\n" == lines[0]

    data_line = lines[1].strip().split("\t")
    assert len(data_line) == 4

def test_protein(inputs_protein):
    """
   Tests the protein alignment processing function.
   It checks if the output file is created with the expected content and format,
   including the similarity score.
    """
    input_file, output_file, similarity_groups = inputs_protein
    prog.process_alignment_protein(input_file, output_file, similarity_groups)

    with open(output_file, 'r') as f:
        lines = f.readlines()

    assert len(lines) == 67
    assert "Seq_1\tSeq_2\tPer_Id\tPer_Sim\tLength\n" == lines[0]

    data_line = lines[1].strip().split("\t")
    assert len(data_line) == 5

def test_load_data_dna(tmp_path):
    """
    Tests the data loading function for DNA data.
    It checks if the function correctly loads the data from a TSV file into a Pandas DataFrame
    and identifies the data type as 'dna'.
    """
    file_path = tmp_path / "test_dna.tsv"
    data = {'Seq_1': ['A', 'B'], 'Seq_2': ['C', 'D'], 'Per_Id': [90, 80], 'Length': [100, 200]}
    df = pd.DataFrame(data)
    df.to_csv(file_path, sep='\t', index=False)

    loaded_df, data_type = prog.load_data(str(file_path))

    assert data_type == 'dna'
    pd.testing.assert_frame_equal(loaded_df, df)

def test_process_alignment_dna_file_not_found(capsys):
    """
    Tests the DNA alignment processing function's behavior when the input file is not found.
    It checks if the function raises a FileNotFoundError exception.
    """
    with pytest.raises(FileNotFoundError):
        prog.process_alignment_DNA("nonexistent_file.fasta", "output.tsv")

def test_process_alignment_dna_invalid_fasta(tmp_path):
    """
    Tests the DNA alignment processing function's behavior when the input file is an invalid FASTA file.
    It checks if the function raises an Exception.
    """
    invalid_fasta_file = tmp_path / "invalid.fasta"
    with open(invalid_fasta_file, "w") as f:
        f.write("This is not a valid FASTA file")

    output_file = tmp_path / "output.tsv"

    with pytest.raises(Exception):
        prog.process_alignment_DNA(str(invalid_fasta_file), str(output_file))

def test_calculate_identity_similarity_protein_empty_sequences():
    """
    Tests the function for calculating identity and similarity for protein sequences when the sequences are empty.
    It checks if both identity and similarity are 0.
    """
    identity, similarity = prog.calculate_identity_similarity_protein("", "", [("A", "G")])
    assert identity == 0
    assert similarity == 0

def test_calculate_identity_similarity_protein_different_lengths():
    """
    Tests the function for calculating identity and similarity for protein sequences when sequences are different length.
    """
    seq1 = "ATGC"
    seq2 = "ATG"

    identity, similarity = prog.calculate_identity_similarity_protein(seq1, seq2, [("A", "G")])
    print(identity, similarity)

def test_process_alignment_protein_invalid_fasta(tmp_path):
    """
    Tests the protein alignment processing function's behavior when the input file is an invalid FASTA file.
    It checks if the function raises an Exception.
    """
    invalid_fasta_file = tmp_path / "invalid.fasta"
    with open(invalid_fasta_file, "w") as f:
        f.write("This is not a valid FASTA file")
    output_file = tmp_path / "output.tsv"
    similarity_groups = [("A", "G")]
    with pytest.raises(Exception):
        prog.process_alignment_protein(str(invalid_fasta_file), str(output_file), similarity_groups)

def test_parse_arguments_dna(monkeypatch):
    """
    Tests the argument parsing function for DNA mode.
    It checks if the function correctly parses the command-line arguments for DNA alignment.
    """
    monkeypatch.setattr(sys, 'argv', ['script_name.py', '-i', 'input.fasta', '-m', 'DNA', '-o', 'output.tsv'])
    args = prog.parse_arguments()
    assert args.input == 'input.fasta'
    assert args.mode == 'DNA'
    assert args.output == 'output.tsv'
    assert args.aminoacid is None

def test_parse_arguments_protein(monkeypatch):
    """
   Tests the argument parsing function for protein mode.
   It checks if the function correctly parses the command-line arguments for protein alignment,
   including the amino acid similarity groups.
   """
    monkeypatch.setattr(sys, 'argv', ['script_name.py', '-i', 'input.fasta', '-m', 'protein', '-o', 'output.tsv', '-amk', 'GA,ST'])
    args = prog.parse_arguments()
    assert args.input == 'input.fasta'
    assert args.mode == 'protein'
    assert args.output == 'output.tsv'
    assert args.aminoacid == 'GA,ST'

def test_parse_arguments_protein_missing_aminoacid(monkeypatch):
    """
   Tests the argument parsing function for protein mode when the amino acid similarity groups are missing.
   It checks if the function raises a SystemExit exception with the correct error code.
   """
    monkeypatch.setattr(sys, 'argv', ['script_name.py', '-i', 'input.fasta', '-m', 'protein', '-o', 'output.tsv'])
    with pytest.raises(SystemExit) as excinfo:
        prog.parse_arguments()
    assert excinfo.value.code == 2

def test_create_matrices_dna(sample_data):
    """
    Tests the function for creating identity matrices for DNA data.
    It checks if the function correctly creates the identity matrix from the input data.
    """
    identity_matrix, _ = prog.create_matrices(sample_data, data_type='dna')

    expected_identity_matrix = pd.DataFrame(100, index=['A', 'B', 'C'], columns=['A', 'B', 'C'], dtype='float64')
    expected_identity_matrix.loc['A', 'B'] = 90
    expected_identity_matrix.loc['B', 'A'] = 90
    expected_identity_matrix.loc['A', 'C'] = 85
    expected_identity_matrix.loc['C', 'A'] = 85
    expected_identity_matrix.loc['B', 'C'] = 80
    expected_identity_matrix.loc['C', 'B'] = 80

    pd.testing.assert_frame_equal(identity_matrix, expected_identity_matrix)

def test_create_matrices_protein(protein_data):
    """
    Tests the function for creating identity and similarity matrices for protein data.
    It checks if the function correctly creates both matrices from the input data.
    """
    identity_matrix, similarity_matrix = prog.create_matrices(protein_data, 'protein')
    assert isinstance(identity_matrix, pd.DataFrame)
    assert isinstance(similarity_matrix, pd.DataFrame)
    assert identity_matrix.loc['A', 'B'] == 0.8
    assert similarity_matrix.loc['A', 'B'] == 0.6

def test_plot_heatmaps_dna(tmp_path):
    """
    Tests the function for plotting heatmaps for DNA data.
    It checks if the function correctly creates the identity heatmap and does not create a similarity heatmap.
    """
    results_dir = tmp_path / "results"
    results_dir.mkdir()
    os.chdir(tmp_path)

    id_matrix = pd.DataFrame([[100, 90], [90, 100]], index=['seq1', 'seq2'], columns=['seq1', 'seq2'])

    prog.plot_heatmaps(id_matrix, None, 'dna')

    identity_files = glob.glob(os.path.join("results", "identity_*.png"))
    assert len(identity_files) == 1

    similarity_files = glob.glob(os.path.join("results", "similarity_*.png"))
    assert len(similarity_files) == 0
    os.chdir(os.path.dirname(__file__))

def test_plot_heatmaps_protein(tmp_path, no_show):
    """
    Tests the function for plotting heatmaps for protein data.
    It checks if the function correctly creates both identity and similarity heatmaps.
    """
    results_dir = tmp_path / "results"
    results_dir.mkdir()
    os.chdir(tmp_path)

    id_matrix = pd.DataFrame([[100, 90], [90, 100]], index=['seq1', 'seq2'], columns=['seq1', 'seq2'])
    sim_matrix = pd.DataFrame([[100, 95], [95, 100]], index=['seq1', 'seq2'], columns=['seq1', 'seq2'])

    prog.plot_heatmaps(id_matrix, sim_matrix, 'protein')

    identity_files = glob.glob(os.path.join("results", "identity_*.png"))
    similarity_files = glob.glob(os.path.join("results", "similarity_*.png"))
    assert len(identity_files) == 1
    assert len(similarity_files) == 1

    os.chdir(os.path.dirname(__file__))

def test_plot_heatmaps_creates_results_directory_if_not_exists(tmp_path, no_show):
    """
    Tests that the plot_heatmaps function creates the results directory if it doesn't already exist.
    """
    os.chdir(tmp_path)

    id_matrix = pd.DataFrame([[100, 90], [90, 100]], index=['seq1', 'seq2'], columns=['seq1', 'seq2'])

    prog.plot_heatmaps(id_matrix, None, 'dna')

    assert os.path.exists("results")

    os.chdir(os.path.dirname(__file__))

def test_main_dna(monkeypatch, capsys, tmp_path):
    """
    Tests the main function for DNA mode.
    It mocks the underlying functions and checks if they are called with the correct arguments.
    """
    input_file = tmp_path / "input.fasta"
    input_file.write_text(">seq1\nATGC\n>seq2\nTGCA")
    output_file = tmp_path / "output.tsv"
    monkeypatch.setattr(sys, 'argv', ['script_name.py', '-i', str(input_file), '-m', 'DNA', '-o', str(output_file)])

    with patch('aps_dop.APS_DoP.process_alignment_DNA') as mock_process_dna, \
         patch('aps_dop.APS_DoP.load_data') as mock_load_data, \
         patch('aps_dop.APS_DoP.create_matrices') as mock_create_matrices, \
         patch('aps_dop.APS_DoP.plot_heatmaps') as mock_plot_heatmaps:

        mock_load_data.return_value = (pd.DataFrame({'Seq_1': [], 'Seq_2': [], 'Per_Id': [], 'Length': []}), 'dna')
        mock_create_matrices.return_value = (pd.DataFrame(), None)

        prog.main()

        mock_process_dna.assert_called_once_with(str(input_file), str(output_file))
        mock_load_data.assert_called_once_with(str(output_file))
        mock_create_matrices.assert_called_once_with(mock_load_data.return_value[0], 'dna')
        mock_plot_heatmaps.assert_called_once_with(mock_create_matrices.return_value[0], mock_create_matrices.return_value[1], 'dna')

def test_main_protein(monkeypatch, capsys, tmp_path):
    """
   Tests the main function for protein mode.
   It mocks the underlying functions and checks if they are called with the correct arguments,
   including the amino acid similarity groups.
   """
    input_file = tmp_path / "input.fasta"
    input_file.write_text(">seq1\nATGC\n>seq2\nTGCA")
    output_file = tmp_path / "output.tsv"
    monkeypatch.setattr(sys, 'argv', ['script_name.py', '-i', str(input_file), '-m', 'protein', '-o', str(output_file), '-amk', 'GA,ST'])

    with patch('aps_dop.APS_DoP.process_alignment_protein') as mock_process_protein, \
         patch('aps_dop.APS_DoP.load_data') as mock_load_data, \
         patch('aps_dop.APS_DoP.create_matrices') as mock_create_matrices, \
         patch('aps_dop.APS_DoP.plot_heatmaps') as mock_plot_heatmaps:

        mock_load_data.return_value = (pd.DataFrame({'Seq_1': [], 'Seq_2': [], 'Per_Id': [], 'Per_Sim': [], 'Length': []}), 'protein')
        mock_create_matrices.return_value = (pd.DataFrame(), pd.DataFrame())

        prog.main()

        mock_process_protein.assert_called_once()
        mock_load_data.assert_called_once_with(str(output_file))
        mock_create_matrices.assert_called_once_with(mock_load_data.return_value[0], 'protein')
        mock_plot_heatmaps.assert_called_once_with(mock_create_matrices.return_value[0], mock_create_matrices.return_value[1], 'protein')

def test_main_protein_missing_aminoacid(monkeypatch, capsys, tmp_path):
    """
   Tests the main function for protein mode when the amino acid similarity groups are missing.
   It checks if the function raises a SystemExit exception and prints the correct error message.
   """
    monkeypatch.setattr(sys, 'argv', ['script_name.py', '-i', 'input.fasta', '-m', 'protein', '-o', 'output.tsv'])
    with pytest.raises(SystemExit) as excinfo:
        prog.main()
    assert excinfo.value.code == 2
    captured = capsys.readouterr()
    assert "error: --aminoacid is required when mode is 'protein'" in captured.err

def test_main_dna_exception(monkeypatch, capsys, tmp_path):
    """
    Tests the main function for DNA mode when an exception occurs during DNA alignment processing.
    It checks if the function catches the exception and raises a new exception with a descriptive message.
    """
    input_file = tmp_path / "input.fasta"
    input_file.write_text(">seq1\nATGC\n>seq2\nTGCA")
    output_file = tmp_path / "output.tsv"
    monkeypatch.setattr(sys, 'argv', ['script_name.py', '-i', str(input_file), '-m', 'DNA', '-o', str(output_file)])

    with patch('aps_dop.APS_DoP.process_alignment_DNA', side_effect=Exception("DNA processing failed")):
        with pytest.raises(Exception, match="Error processing DNA alignment: DNA processing failed"):
            prog.main()

def test_main_data_processing_exception(monkeypatch, capsys, tmp_path):
    """
    Tests the main function when an exception occurs during data loading.
    It asserts that an exception with a descriptive message is raised.
   """
    input_file = tmp_path / "input.fasta"
    input_file.write_text(">seq1\nATGC\n>seq2\nTGCA")
    output_file = tmp_path / "output.tsv"
    monkeypatch.setattr(sys, 'argv', ['script_name.py', '-i', str(input_file), '-m', 'DNA', '-o', str(output_file)])

    with patch('aps_dop.APS_DoP.process_alignment_DNA') as mock_process_dna, \
         patch('aps_dop.APS_DoP.load_data', side_effect=Exception("Data processing failed")), \
         patch('aps_dop.APS_DoP.create_matrices') as mock_create_matrices, \
         patch('aps_dop.APS_DoP.plot_heatmaps') as mock_plot_heatmaps:
        with pytest.raises(Exception, match="Error during data processing or plotting: Data processing failed"):
            prog.main()

def test_main_protein_exception(monkeypatch, capsys, tmp_path):
    """
    Tests the main function for protein mode when an exception occurs during protein alignment processing.
    It checks if the function catches the exception and raises a new exception with a descriptive message.
    """
    input_file = tmp_path / "input.fasta"
    input_file.write_text(">seq1\nATGC\n>seq2\nTGCA")
    output_file = tmp_path / "output.tsv"
    monkeypatch.setattr(sys, 'argv', ['script_name.py', '-i', str(input_file), '-m', 'protein', '-o', str(output_file), '--aminoacid', 'A, B'])

    with patch('aps_dop.APS_DoP.process_alignment_protein', side_effect=Exception("Protein processing failed")):
        with pytest.raises(Exception, match="Error processing protein alignment: Protein processing failed"):
            prog.main()

def test_file_not_found():
    """
    Tests that FileNotFoundError is raised when the specified file doesn't exist.
    """
    with pytest.raises(FileNotFoundError) as excinfo:
        prog.process_alignment_protein("nonexistent_file.fasta", "output.txt", {})
    assert str(excinfo.value) == "File 'nonexistent_file.fasta' not found."

def test_less_than_two_sequences(capsys):
    """
   Tests that the program returns an error message when the specified file contain less then two sequences.
   """
    with patch("Bio.AlignIO.read", return_value=[SeqRecord(Seq("ACGT"), id="seq1")]):
        prog.process_alignment_protein("dummy_file.fasta", "output.txt", {})
    captured = capsys.readouterr()
    assert captured.out.strip() == "Error: The alignment file must contain at least two sequences."

def test_load_data_file_not_found():
    """
    Tests that FileNotFoundError is raised when the specified file for load_data function doesn't exist.
    """
    input_file = "nonexistent_file.txt"
    with pytest.raises(FileNotFoundError) as excinfo:
        prog.load_data(input_file)
    assert str(excinfo.value) == f"File {input_file} doesn't found!"


