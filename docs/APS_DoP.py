from Bio import AlignIO
from itertools import combinations
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
from datetime import datetime

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
    print(f"Processing DNA alignment from {alignment_file} to {output_file}")
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

def load_data(input_file):
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Файл {input_file} не найден!")
    df = pd.read_csv(input_file, sep='\t')
    data_type = 'protein' if 'Per_Sim' in df.columns else 'dna'
    print(f"Тип данных: {data_type.upper()}")
    return df, data_type

def create_matrices(df, data_type):
    all_ids = sorted(set(df['Seq_1']).union(set(df['Seq_2'])))
    if data_type == 'dna':
        identity_matrix = pd.DataFrame(100, index=all_ids, columns=all_ids)
        for _, row in df.iterrows():
            id1, id2 = row['Seq_1'], row['Seq_2']
            identity_matrix.loc[id1, id2] = identity_matrix.loc[id2, id1] = row['Per_Id']
        return identity_matrix, None
    else:
        identity_matrix = pd.DataFrame(100, index=all_ids, columns=all_ids)
        similarity_matrix = pd.DataFrame(100, index=all_ids, columns=all_ids)
        for _, row in df.iterrows():
            id1, id2 = row['Seq_1'], row['Seq_2']
            identity_matrix.loc[id1, id2] = row['Per_Id']
            similarity_matrix.loc[id1, id2] = row['Per_Sim']
        return identity_matrix, similarity_matrix


def plot_heatmaps(id_matrix, sim_matrix, data_type):
    os.makedirs("results", exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    # Identity Graph
    matrix_size = id_matrix.shape[0]
    figsize = (max(10, matrix_size * 0.8), max(8, matrix_size * 0.6))  # Adjust factors
    plt.figure(figsize=figsize)
    sns.heatmap(id_matrix, annot=True, fmt=".1f", cmap="YlOrRd", vmin=0, vmax=100)
    plt.title("Pairwise Identity (%)")
    identity_file = f"results/identity_{timestamp}.png"
    plt.savefig(identity_file, bbox_inches='tight', dpi=300)
    plt.close()
    print(f"Сохранен: {identity_file}")

    # Similarity plot (for proteins only)
    if data_type == 'protein':
        matrix_size = id_matrix.shape[0]
        figsize = (max(10, matrix_size * 0.8), max(8, matrix_size * 0.6))  # Adjust factors
        plt.figure(figsize=figsize)
        sns.heatmap(sim_matrix, annot=True, fmt=".1f", cmap="YlGnBu", vmin=0, vmax=100)
        plt.title("Pairwise Similarity (%)")
        similarity_file = f"results/similarity_{timestamp}.png"
        plt.savefig(similarity_file, bbox_inches='tight', dpi=300)
        plt.close()
        print(f"Сохранен: {similarity_file}")

import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, help='Input file path')
    parser.add_argument('-m', '--mode', required=True, choices=['DNA', 'protein'], help='Mode: DNA or protein')
    parser.add_argument('-o', '--output', default='result_table.tsv', help='Output file path (default: result_table.tsv)')
    parser.add_argument('-amk', '--aminoacid', help='Amino acid (required for protein mode)')
    args = parser.parse_args()

    if args.mode == 'DNA':
        try:
            process_alignment_DNA(args.input, args.output)
        except Exception as e:
            print(f"Error processing DNA alignment: {e}")

    elif args.mode == 'protein':
        if args.aminoacid is None:
            parser.error("--aminoacid is required when mode is 'protein'")
        clear_aminoacid = [tuple(group.strip()) for group in args.aminoacid.split(',')]
        try:
            process_alignment_protein(args.input, args.output, clear_aminoacid)
        except Exception as e:
            print(f"Error processing protein alignment: {e}")

    df, data_type = load_data(args.output)
    id_matrix, sim_matrix = create_matrices(df, data_type)
    plot_heatmaps(id_matrix, sim_matrix, data_type)
