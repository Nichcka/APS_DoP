import aps_dop.APS_DoP as prog

def test_dna(inputs_dna):
    input_file, output_file = inputs_dna
    prog.process_alignment_DNA(input_file, output_file)

    with open(output_file, 'r') as f:
        lines = f.readlines()

    assert len(lines) == 4
    assert "Seq_1\tSeq_2\tPer_Id\tLength\n" == lines[0]

    data_line = lines[1].strip().split("\t")
    assert len(data_line) == 4

def test_protein(inputs_protein):
    input_file, output_file, similarity_groups = inputs_protein
    prog.process_alignment_protein(input_file, output_file, similarity_groups)

    with open(output_file, 'r') as f:
        lines = f.readlines()

    assert len(lines) == 67
    assert "Seq_1\tSeq_2\tPer_Id\tPer_Sim\tLength\n" == lines[0]

    data_line = lines[1].strip().split("\t")
    assert len(data_line) == 5




