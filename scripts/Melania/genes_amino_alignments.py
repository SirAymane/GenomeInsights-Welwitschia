
"""
Genes search and aminoacids alignments

"""

from Bio import SeqIO
from Bio.Seq        import Seq
from Bio.SeqRecord  import SeqRecord
from Bio.SeqFeature import SeqFeature
from Bio.Align      import PairwiseAligner, PairwiseAlignments, PairwiseAlignment, substitution_matrices

from Bio.Align.substitution_matrices import Array





def get_translation(gen_name: str, genbanks: dict): 

    trans_list: list[Seq] = []

    for key in genbanks.keys():
        genbank = genbanks[key]
        #print(genbank)
        for feature in genbank.features:
            if feature.type == "CDS":
                gene: str = feature.qualifiers['gene'][0]
                if gene == gen_name:
                    trans_list.append(Seq(feature.qualifiers['translation'][0]))
    return trans_list


def align_aminos(seq_list: list[Seq], seq1_pos: int, seq2_pos: int):

    # Create Global Aligner
    aligner: PairwiseAligner = PairwiseAligner()

    # Get BLOSUM62 matrix
    blosum62_matrix: Array = substitution_matrices.load('BLOSUM62')


    # Put the matrix in the aligner.
    aligner.substitution_matrix = blosum62_matrix
    
    # Get first global alignment
    alignments: PairwiseAlignments = aligner.align(seq_list[seq1_pos], seq_list[seq2_pos])
    alignment = alignments[0]

    # Print it, its score, etc.
    print(alignment)
    print(f"Score: {alignment.score}")

    print(f"Matched indexes in target seq:  {alignment.aligned[0]}")
    print(f"Matched indexes in query seq:   {alignment.aligned[1]}")

    return alignment




# Main
# ---------------------------------------------------------------------
this_module: str = __name__
main_module: str = "__main__"

if this_module == main_module:

    genbanks: dict[str,SeqRecord] = read_list_genbanks(genbank_list) #welwitschia_mirabilis.gb, ephedra_trifurca.gb, gnetum_montanum.gb
    gen1: str = "matK" # Codifica para una proteína implicada en el fotosistema II de las plantas. La proteína MatK es un componente estructural de la subunidad de reacción del fotosistema II, que es responsable de la transferencia de electrones durante la fotosíntesis. Además de su papel en la fotosíntesis, matK se ha utilizado como un marcador molecular útil para estudios filogenéticos de plantas debido a su conservación a través de una amplia gama de especies y su relativamente alta tasa de evolución. 
    gen2: str = "rbcL" # Codifica la proteína de la subunidad de la cianobacteriana Rubisco, que es la enzima clave en la fotosíntesis.
    gen3: str = "atpF" # Codifica una subunidad de la ATP sintasa, una enzima que produce ATP a partir de la energía química liberada en la fotosíntesis.
    gen4: str = "clpP" # Codifica una proteína que interviene en el proceso de degradación de proteínas en el cloroplasto.

    trans_gen1 = get_translation(gen1, genbanks)
    trans_gen2 = get_translation(gen2, genbanks)
    trans_gen3 = get_translation(gen3, genbanks)
    trans_gen4 = get_translation(gen4, genbanks)

    # Alineamiento matK 1: welwitschia_mirabilis - ephedra trifurca
    align_aminos(trans_gen1, 0, 1)
    # Alineamiento matK 2: welwitschia_mirabilis - gnetum_montanum
    align_aminos(trans_gen1, 0, 2)
    # Alineamiento matK 3: ephedra trifurca - gnetum_montanum
    align_aminos(trans_gen1, 1, 2)

    # Alineamiento rbcL 1: welwitschia_mirabilis - ephedra trifurca
    align_aminos(trans_gen2, 0, 1)
    # Alineamiento rbcL 2: welwitschia_mirabilis - gnetum_montanum
    align_aminos(trans_gen2, 0, 2)
    # Alineamiento rbcL 3: ephedra trifurca - gnetum_montanum
    align_aminos(trans_gen2, 1, 2)

    # Alineamiento atpF 1: welwitschia_mirabilis - ephedra trifurca
    align_aminos(trans_gen3, 0, 1)
    # Alineamiento atpF 2: welwitschia_mirabilis - gnetum_montanum
    align_aminos(trans_gen3, 0, 2)
    # Alineamiento atpF 3: ephedra trifurca - gnetum_montanum
    align_aminos(trans_gen3, 1, 2)

    # Alineamiento psaB 1: welwitschia_mirabilis - ephedra trifurca
    align_aminos(trans_gen4, 0, 1)
    # Alineamiento psaB 2: welwitschia_mirabilis - gnetum_montanum
    align_aminos(trans_gen4, 0, 2)
    # Alineamiento psaB 3: ephedra trifurca - gnetum_montanum
    align_aminos(trans_gen4, 1, 2)

