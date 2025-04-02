from Bio import AlignIO
from itertools import combinations
import sys

def calculate_identity_DNA(alignment1, alignment2):
    """
    Calculates the percent identity between two alignment sequences.
    """
    matches = 0
    alignment_length = len(alignment1)
    for i in range(alignment_length):
        if alignment1[i] == alignment2[i] and alignment1[i] != '-':
            matches += 1

    return (matches / alignment_length) * 100 if alignment_length > 0 else 0


def process_alignment_DNA(alignment_file, output_file):
    """
    Processes a FASTA alignment file and outputs a table of pairwise comparisons.
    Saves the results to the specified file.
    """

    try:
        alignment = AlignIO.read(alignment_file, "fasta")
        sequences = list(alignment)
    except FileNotFoundError:
        print(f"Error: File '{alignment_file}' not found.")
        return
    except Exception as e:
        print(f"Error reading file '{alignment_file}': {e}")
        return

    with open(output_file, "w") as outfile:  # Open output file for writing
        outfile.write("Seq_1\tSeq_2\tPer_Id\tLength\n")  # Write header

        pairs = combinations(sequences, 2)

        for seq1, seq2 in pairs:
            # Calculate the percentage of identity
            percent_identity = calculate_identity_DNA(str(seq1.seq), str(seq2.seq))

            # Output the result in table format and write it to a file
            outfile.write(f"{seq1.id}\t{seq2.id}\t{percent_identity:.2f}\t{len(seq1.seq)}\n")

    print(f"The results are saved to file: {output_file}")

def calculate_identity_similarity_protein(alignment1, alignment2, similarity_groups):
    """
    Calculates percent identity and similarity between two aligned protein sequences,
    using only the provided similarity groups.
    """
    matches = 0
    similar = 0
    alignment_length = len(alignment1)

    if alignment_length == 0:
        return 0, 0

    for i in range(alignment_length):

        if i >= len(alignment1) or i >= len(alignment2):
            continue

        if alignment1[i] == alignment2[i] and alignment1[i] != '-':
            matches += 1
            similar += 1  # Identical is always similar

        elif alignment1[i] != '-' and alignment2[i] != '-':

            for group in similarity_groups:
                if alignment1[i] in group and alignment2[i] in group:
                    similar += 1
                    break  # Only count one group

    identity = (matches / alignment_length) * 100
    similarity = (similar / alignment_length) * 100

    return identity, similarity


def process_alignment_protein(alignment_file, output_file, similarity_groups):
    """Processes a FASTA alignment file and outputs a table of pairwise comparisons."""

    try:
        alignment = AlignIO.read(alignment_file, "fasta")
        sequences = list(alignment)

    except FileNotFoundError:
        print(f"Error: File '{alignment_file}' not found.")
        return
    except Exception as e:
        print(f"Error reading file '{alignment_file}': {e}")
        return

    if len(sequences) < 2:
        print("Error: The alignment file must contain at least two sequences.")
        return

    with open(output_file, "w") as outfile:
        outfile.write("Seq_1\tSeq_2\tPer_Id\tPer_Sim\tLength\n")  # Write header

        pairs = combinations(sequences, 2)

        for seq1, seq2 in pairs:
            identity, similarity = calculate_identity_similarity_protein(str(seq1.seq), str(seq2.seq), similarity_groups)
            outfile.write(f"{seq1.id}\t{seq2.id}\t{identity:.2f}\t{similarity:.2f}\t{len(seq1.seq)}\n")

    print(f"The results are saved to file: {output_file}")   

