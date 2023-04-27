# 4 - Alineament de seqüències
################################
"""
- Utilitzeu BioPython per alinear una de les seqüències amb la resta de seqüències.
- Hi ha d'haver al menys 10 comparacions.
- Guardeu-vos la puntuació de l'alineament i ordeneu les comparacions segons la similitud


1. Obtener localización del gen
2. Obtener seq del gen
3. Alinear con Bio Align
4. Obtener resultados, se pueden guardar en formato Fasta.
"""

# Imports
# --------------------------------------------------------------------------
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner, PairwiseAlignments, PairwiseAlignment, substitution_matrices
from pathlib import Path
from Bio.SeqIO import SeqRecord
from pprint import pprint
from Bio.SeqFeature import SeqFeature
# from parse_genbank import read_list_genbanks, genbank_list


# DIR
# --------------------------------------------------------------------------
DATA_DIR: Path = Path(__file__).parents[2]/"data"

mirabilis_file = DATA_DIR/"welwitschia_mirabilis.gb"
trifurca_file = DATA_DIR/"ephedra_trifurca.gb"
montanum_file = DATA_DIR/"gnetum_montanum.gb"


# Get sequence
# --------------------------------------------------------------------------
def get_gene_sequence(genbank_file, gene_location):

    # For each record in genbank file
    for record in SeqIO.parse(genbank_file, "genbank"):

        # Get gene sequence with the location
        gene_sequence = record.seq[gene_location.start:gene_location.end]

        # Return the gene sequence
        return gene_sequence


# Push the gene seq in the list
# --------------------------------------------------------------------------
def gene_seq(gene_name, sequence_file, mlist):

    # For each record in genbank file
    for record in SeqIO.parse(sequence_file, "genbank"):

        # For each feature:
        for feature in record.features:
            # if feature.type == "CDS":

            # If the feature type is CDS and the gene name is the given as a parameter:
            if feature.type == "CDS" and feature.qualifiers.get('gene')[0] == gene_name:

                # Get gene location
                gene_location = feature.location

                # And get the gene sequence
                gene_sequence = get_gene_sequence(sequence_file, gene_location)

                # Next, append the sequence in the list
                mlist.append(gene_sequence)

                # print(gene_location)
                # return gene_sequence

# Function to insert the genes seq in a list
def insert_gene_seq(gene_name):
    mylist = []

    gene_seq(gene_name, mirabilis_file, mylist)
    gene_seq(gene_name, trifurca_file, mylist)
    gene_seq(gene_name, montanum_file, mylist)

    return mylist


# Function to align the gene
# --------------------------------------------------------------------------
def align_genes(seq1, seq2):

    # Create Global Aligner
    aligner: PairwiseAligner = PairwiseAligner()

    # Get BLOSUM62 matrix
    blosum62_matrix: Array = substitution_matrices.load('BLOSUM62')

    # Put the matrix in the aligner.
    aligner.substitution_matrix = blosum62_matrix
    
    # Get first global alignment
    alignments: PairwiseAlignments = aligner.align(seq1, seq2)
    alignment = alignments[0]

    # Print it, its score, etc.
    print(alignment)
    print(f"Score: {alignment.score}")
    print(f"Num of aligned seqs:  {len(alignment)}")
    # Shape
    print(f"Alignment shape:       {alignment.shape}")
    print(f"Num of aligned seqs:    {alignment.shape[0]}")
    print(f"Num of aligned letters: {alignment.shape[1]}")
    print(f"Matched indexes in target seq:  {alignment.aligned[0]}")
    print(f"Matched indexes in query seq:   {alignment.aligned[1]}")
    return alignment


# Main
# -----------------------------------------------------------------------
if __name__ == "__main__":

    ##################### Get seqs in list ####################

    # Genes
    gene1 = "matK"
    gene2 = "rbcL"
    gene3 = "atpF"
    gene4 = "psaB"

    # Get gene seq of psbA
    gene1_seq_list = insert_gene_seq(gene1)

    # Get gene seq of rbcL
    gene2_seq_list = insert_gene_seq(gene2)

    # Get gene seq of atpF
    gene3_seq_list = insert_gene_seq(gene3)

    # Get gene seq of psaB
    gene4_seq_list = insert_gene_seq(gene4)

    # Checks
    # print(gene2_seq_list)
    # print(len(gene1_seq_list[2]))


    # Constants
    MIRABILIS = 0
    TRIFURCA  = 1
    MONTANUM  = 2

    # print(len(list1[2]))

    ##################### ALIGNMENT ####################
    # psbA
    # align_genes(gene1_seq_list[MIRABILIS], gene1_seq_list[TRIFURCA])
    # align_genes(gene1_seq_list[MIRABILIS], gene1_seq_list[MONTANUM])


    # rbcL
    # align_genes(gene2_seq_list[MIRABILIS], gene2_seq_list[TRIFURCA])
    # align_genes(gene2_seq_list[MIRABILIS], gene2_seq_list[MONTANUM])

    # # atpF
    # align_genes(gene3_seq_list[MIRABILIS], gene3_seq_list[TRIFURCA])
    # align_genes(gene3_seq_list[MIRABILIS], gene3_seq_list[MONTANUM])

    # psaB
    # align_genes(gene4_seq_list[MIRABILIS], gene4_seq_list[TRIFURCA])
    # align_genes(gene4_seq_list[MIRABILIS], gene4_seq_list[MONTANUM])