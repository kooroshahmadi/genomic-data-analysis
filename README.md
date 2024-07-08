# genomic-data-analysis
This project performs basic genomic data analysis using Biopython.  It includes tasks such as sequence alignment and mutation detection.

### Explanation

#### Header and Imports
"""
Genomic Data Analysis Project

This project performs basic genomic data analysis using Biopython. 
It includes tasks such as sequence alignment and mutation detection.

Requirements:
- Biopython
- Clustal Omega (for sequence alignment)
"""

from Bio import SeqIO
from Bio.Align.Applications import ClustalOmegaCommandline

- The first section is a multi-line comment (docstring) describing the purpose of the script, its tasks, and its requirements.
- `from Bio import SeqIO`: Imports the `SeqIO` module from Biopython, which is used for reading and writing sequence file formats.
- `from Bio.Align.Applications import ClustalOmegaCommandline`: Imports the Clustal Omega command-line wrapper from Biopython for sequence alignment.

#### Function: `read_fasta`
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

- `def read_fasta(file_path)`: Defines a function to read a FASTA file.
- `sequences = list(SeqIO.parse(file_path, "fasta"))`: Uses `SeqIO.parse` to read the FASTA file and converts the resulting iterator to a list of sequences.
- `return sequences`: Returns the list of sequences.

#### Function: `align_sequences`
def align_sequences(input_file, output_file):
    """
    Aligns sequences using Clustal Omega.

    Args:
        input_file (str): Path to the input FASTA file.
        output_file (str): Path to the output aligned file.
    """
    clustalomega_cline = ClustalOmegaCommandline(
        infile=input_file, force=True, outfile=output_file, verbose=True, auto=True
    )
    clustalomega_cline()

- `def align_sequences(input_file, output_file)`: Defines a function to align sequences using Clustal Omega.
- `clustalomega_cline = ClustalOmegaCommandline(...)`: Creates a Clustal Omega command-line object with specified input and output files.
- `clustalomega_cline()`: Executes the Clustal Omega command-line tool.

#### Function: `detect_mutations`
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
        if seq1[i] != seq2[i]):
            mutations.append((i, seq1[i], seq2[i]))
    return mutations

- `def detect_mutations(seq1, seq2)`: Defines a function to detect mutations between two sequences.
- `mutations = []`: Initializes an empty list to store mutation positions.
- `for i in range(len(seq1))`: Loops through each position in the first sequence.
- `if seq1[i] != seq2[i]`: Checks if the characters at position `i` are different between the two sequences.
- `mutations.append((i, seq1[i], seq2[i]))`: Adds the position and differing characters to the mutations list.
- `return mutations`: Returns the list of mutation positions.

#### Main Block
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

- `if __name__ == "__main__":`: Ensures the code within this block runs only if the script is executed directly, not when imported as a module.
- `input_fasta = "example_sequences.fasta"`: Sets the input FASTA file path.
- `output_aligned = "aligned_sequences.aln"`: Sets the output aligned file path.
- `sequences = read_fasta(input_fasta)`: Reads the sequences from the input FASTA file.
- `print(f"Read {len(sequences)} sequences from {input_fasta}")`: Prints the number of sequences read.
- `align_sequences(input_fasta, output_aligned)`: Aligns the sequences using Clustal Omega.
- `print(f"Aligned sequences saved to {output_aligned}")`: Prints a message indicating the output file location.
- `seq1 = str(sequences[0].seq)`: Converts the first sequence to a string.
- `seq2 = str(sequences[1].seq)`: Converts the second sequence to a string.
- `mutations = detect_mutations(seq1, seq2)`: Detects mutations between the two sequences.
- `print(f"Detected mutations: {mutations}")`: Prints the detected mutations.
