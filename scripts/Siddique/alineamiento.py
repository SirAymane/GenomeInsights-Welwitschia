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
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner, PairwiseAlignments, PairwiseAlignment, substitution_matrices
from pathlib import Path
from Bio.SeqIO import SeqRecord
from pprint import pprint

from parse_genbank import read_list_genbanks, genbank_list


# Get welwitschia_mirabilis
records = SeqIO.read('/bio/exercice/m14_biopython_pt1/welwitschia_mirabilis.gb', 'genbank')
for record in SeqIO.parse("/bio/exercice/m14_biopython_pt1/welwitschia_mirabilis.gb", "genbank"):
  #Imprimir la posición de la secuencia
  print(record.features[0].location)


#Leer el archivo de GenBank
for record in SeqIO.parse("/bio/exercice/m14_biopython_pt1/welwitschia_mirabilis.gb", "genbank"):
  #Comprobar si el registro es para el gen psbA
  if record.features[0].qualifiers["gene"] == 'psbA':
    #Imprimir la posición de la secuencia
    # Print locations
    print(record.features[0].location)
    # Print seq
    print(record.seq)



# Show all the genes in the record
def show_genes(records):
    for feature in records.features:
        if feature.type == 'gene':
            print(feature.qualifiers['gene'][0])
            # if (feature.qualifiers['gene']=='psbA'):
            #     print('ok')


# Get seq by the index of the gen
def get_gen_seq(param1, param2):
    seq = record.seq[param1:param2]
    return seq
    


seqpsbA = get_gen_seq(1, 899)
#print(seqpsbA)


# Gen psbA gen seq
seq1 = (record.seq[1:899])
#print(seq1)
#print(len(seq1))


# Gen rbcL gen seq
seq2 = (record.seq[42893:44320])
#print(len(seq2))
#print(seq2)



# Align globally
def align_globally(seq1: Seq, seq2: Seq):
    # Create Aligner
    aligner: PairwiseAligner = PairwiseAligner()

    # Get score only
    score: float = aligner.score(seq1, seq2)
    print(f"Global alignment score: {score}")

    # Get global alignments. Global is the default.
    alignments: PairwiseAlignments = aligner.align(seq1, seq2)

    # Alignments looks like an iterator, but it is not.
    # We can ask the length and has some extra methods. See dir().
    print(f"Number of alignments: {len(alignments)}")

    # Each alignment is a PairwiseAlignment object. It has additional methods.
    alignment: PairwiseAlignment
    for alignment in alignments:
        print(alignment)
    

# Local alignment. Score and number of alignments differ from global alignment.
# ---------------------------------------------------------------------
def align_locally(seq1: Seq, seq2: Seq):

    # Create Aligner
    aligner: PairwiseAligner = PairwiseAligner()
    aligner.mode = 'local'

    # Get score only
    score: float = aligner.score(seq1, seq2)
    print(f"Local alignment score: {score}")

    # Get global alignments. Global is the default.
    alignments: PairwiseAlignments = aligner.align(seq1, seq2)

    # Alignments looks like an iterator, but it is not.
    print(f"Number of alignments: {len(alignments)}")

    # Each alignment is a PairwiseAlignment object. It has additional methods.
    alignment: PairwiseAlignment
    for alignment in alignments:
        print(alignment)


def inspect_global_alignment(param1, param2):

    # Create Aligner
    aligner: PairwiseAligner = PairwiseAligner()

    # Seqs
    seq1: Seq = param1
    seq2: Seq = param2

    # Get global alignments.
    alignments: PairwiseAlignments = aligner.align(seq1, seq2)

    # Get first alignment
    alignment: PairwiseAlignment = alignments[0]

    # Alignment attributes
    print(f"Score: {alignment.score}")
    print(f"Target seq: {alignment.target}")
    print(f"Query seq:  {alignment.query}")
    print(f"Num of aligned seqs:  {len(alignment)}")
    print()

    # Print alignment
    print(alignment)

    # Shape
    print(f"Alignment shape:       {alignment.shape}")
    print(f"Num of aligned seqs:    {alignment.shape[0]}")
    print(f"Num of aligned letters: {alignment.shape[1]}")
    print()

    # Aligned
    print(f"Aligned attribute:      {alignment.aligned}")
    print(f"Matched indexes in target seq:  {alignment.aligned[0]}")
    print(f"Matched indexes in query seq:   {alignment.aligned[1]}")



# Main
# ---------------------------------------------------------------------
this_module: str = __name__
main_module: str = "__main__"

if this_module == main_module:

    #prueba = read_list_genbanks([Ephedra_trifurca_genbank, Gnetum_montanum_genbank, Welwitschia_mirabilis_genbank])

    genbanks: dict[str,SeqRecord] = read_list_genbanks(genbank_list)
 
    #pprint(genbanks)

    # show_genes(record)
    #align_globally(seq1,seq2)
    #align_locally()
# ---------------------------------------------------------------------