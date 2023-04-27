from Bio import Entrez
from io import TextIOWrapper
from pathlib import Path

Entrez.email = "welwitschia@mirabilis.edu"

species: list[dict] = [
    {
        "file_name": "welwitschia_mirabilis.gb",
        "accession_id": "NC_010654"
    },
    {
        "file_name": "ephedra_trifurca.gb",
        "accession_id": "NC_065630"
    },
    {
        "file_name": "gnetum_montanum.gb",
        "accession_id": "NC_021438"
    }
]

def save_gb(gb: str, file_path: Path) -> None:
    """
    Save a genbank as text file

    Args:
        gb: genbank as text
        file_path (Path): destination path

    Returns:
        None
    """
        
    with open(file_path, "w") as file:
        file.write(gb)
        print(f"{file_path.name} has been exported!")


def get_gb_by_accession(species: list[dict]) -> None:
    """
    Given a list of species (accession id), returns the
    genbank text of NCBI.
    The function only queries for the genbank if the file
    is not present locally. Otherwise no query is make

    Args:
        species (list[str]): accession names of species to search

    Returns:
        None
    """

    data_path: Path = Path(__file__).parents[2]/"data"

    for specie in species:

        file_path: Path = data_path/specie["file_name"]

        if not file_path.exists():

            handle: TextIOWrapper = Entrez.efetch(db="nucleotide", id=specie["accession_id"], rettype="gb", retmode="text")

            gb: str = ""
            for line in handle:
                gb += line

            save_gb(gb,file_path)

        else:
            print(f"{file_path.name} already exists!")



if __name__ == "__main__":

    get_gb_by_accession(species)
