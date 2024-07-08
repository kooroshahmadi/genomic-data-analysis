
# genomic_data_analysis.py

"""
Genomic Data Analysis Project

This project performs basic genomic data analysis using Biopython. It includes tasks such as sequence alignment and mutation detection.

Requirements:
- Biopython
- Clustal Omega (for sequence alignment)

To install Biopython, use:
    pip install biopython
"""

from Bio import SeqIO
from Bio.Align.Applications import ClustalOmegaCommandline

def read_fasta(file_path):
    """
    Reads a FASTA file and returns the sequences.
    
    Args:
        file_path (str): Path to the FASTA file.
    
    Returns:
        list: List of sequences.
    """
    sequences = list(SeqIO.parse(file_path, "fasta"))
    return sequences

def align_sequences(input_file, output_file):
    """
    Aligns sequences using Clustal Omega.
    
    Args:
        input_file (str): Path to the input FASTA file.
        output_file (str): Path to the output aligned file.
    """
    clustalomega_cline = ClustalOmegaCommandline(infile=input_file, force=True, outfile=output_file, verbose=True, auto=True)
    clustalomega_cline()

def detect_mutations(seq1, seq2):
    """
    Detects mutations between two sequences.
    
    Args:
        seq1 (str): First sequence.
        seq2 (str): Second sequence.
    
    Returns:
        list: List of mutation positions.
    """
    mutations = []
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            mutations.append((i, seq1[i], seq2[i]))
    return mutations

if __name__ == "__main__":
    # Example usage
    input_fasta = "example_sequences.fasta"
    output_aligned = "aligned_sequences.aln"
    
    # Read sequences
    sequences = read_fasta(input_fasta)
    print(f"Read {len(sequences)} sequences from {input_fasta}")
    
    # Align sequences
    align_sequences(input_fasta, output_aligned)
    print(f"Aligned sequences saved to {output_aligned}")
    
    # Detect mutations (example for two sequences)
    seq1 = str(sequences[0].seq)
    seq2 = str(sequences[1].seq)
    mutations = detect_mutations(seq1, seq2)
    print(f"Detected mutations: {mutations}")
