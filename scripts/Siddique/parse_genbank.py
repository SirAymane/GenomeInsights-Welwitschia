from Bio                import SeqIO            # SeqIO is a module
from Bio.Seq            import Seq              # Seq is a class
from Bio.SeqRecord      import SeqRecord        # SeqRecord is a class
from Bio.SeqFeature     import SeqFeature

from pathlib import Path
from pprint import pp


# Constants
DATA_DIR: Path = Path(__file__).parents[2]/"data"

ephedra_trifurca_genbank         = DATA_DIR/'ephedra_trifurca.gb'
gnetum_montanum_genbank          = DATA_DIR/'gnetum_montanum.gb'
welwitschia_mirabilis_genbank    = DATA_DIR/'welwitschia_mirabilis.gb'

genbank_list: list[Path] = [ephedra_trifurca_genbank,gnetum_montanum_genbank,welwitschia_mirabilis_genbank]



# #############################################################################
# - PT1 Biopython - Aymane
# #############################################################################


def read_list_genbanks(genbank_list : list) -> dict[str,SeqRecord]: 



    genbank_dict: dict = {}

    # 1. Show what information you have in the first part of the genbank
    for genbank_path in genbank_list:
        # 1. Show what information you have in the first part of the genbank
        list_record: SeqRecord = SeqIO.read(genbank_path, 'genbank')
        genbank_dict [genbank_path.name] = list_record
    
    return genbank_dict

    


# Main
# -----------------------------------------------------------------------------
this_module = __name__
main_module = '__main__'

if (this_module == main_module):

    genbanks: dict[str,SeqRecord] = read_list_genbanks(genbank_list)

    for genbank in genbanks.keys():
        pp(genbanks[genbank].annotations.keys())
        # Show all annotations
        pp(genbanks[genbank].annotations)
        # Show main reference
        pp(genbanks[genbank].annotations['references'][0])
    