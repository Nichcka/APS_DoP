## APS_DoP (Alignment Pairwise Statistics DNA or Protein)

### version 1.0

APS_DoP is a command-line tool for analyzing multiple sequence alignments of DNA or protein sequences. It provides detailed information on sequence identity and similarity, useful for studying evolutionary relationships, domain conservation, and other bioinformatics tasks.

## Installation

```
git clone https://github.com/Nichcka/APS_DoP
cd APS_DoP
pip install poetry
poetry install
```

Run tests to ensure correct program functionality: `pytest`

## Key Features
Identity and Similarity Analysis: APS_DoP calculates the percentage of identity and similarity for each sequence pair in the provided multiple sequence alignment. Similarity is computed based on a user-defined set of residue groups (e.g., "GAVLI, FYW, CM, ST, KRH, DENQ, P" for amino acids), where residues within each group are considered equivalent.

Pairwise Comparisons: The tool performs pairwise comparisons of all sequences in the alignment, providing a comprehensive overview of the relationships between them.

Tabular Output: APS_DoP generates a tab-separated values (.tsv) table containing the following information for each sequence pair:

- Column 1: Identifier of sequence 1 
- Column 2: Identifier of sequence 2 
- Column 3: Percent identity between sequences 1 and 2, calculated from the alignment 
- Column 4: Percent similarity between sequences, calculated based on the alignment and the defined residue groups (residues within a group are considered equivalent)
- Column 5: Total alignment length of the two sequences

Heatmaps (Plots): The tool generates identity and similarity heatmaps, visualizing the pairwise sequence comparisons. This allows for quick identification of clusters of similar sequences and the detection of patterns.

## Input
A multi-FASTA file containing a multiple sequence alignment of DNA or protein sequences, and comma-separated groups of similar residues (e.g., "GAVLI, FYW, CM, ST, KRH, DENQ, P" for amino acids). 

## Quick usage
for DNA:

``python aps_dop/APS_DoP.py -i input.fasta -m DNA -o output.tsv``

for protein:

``python aps_dop/APS_DoP.py -i input.fasta  -m protein -o output.tsv -amk GA,ST,NDEQ,KRH,VILM,YWF,CP``

## Contributers
A project for the "Python Programming" course at St. Petersburg State University

Code contributors:
- Vinichenko Veronika
- Vlasevskaya Anastasia
- Tsapulina Ekaterina
