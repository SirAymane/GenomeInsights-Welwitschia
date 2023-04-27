import re
from pathlib import Path
from Bio.SeqIO import SeqRecord
from Bio.SeqFeature import Reference
from pprint import pprint
from parse_genbank import genbank_list, read_list_genbanks
from pprint import pp
from Bio import SeqIO


# #############################################################################
# - PT1 Biopython - Aymane
# #############################################################################


genbanks: dict[str, SeqRecord] = read_list_genbanks(genbank_list)


features = []


def get_features(gb: dict[str, SeqRecord]) -> list[str]:
    """ Applies a regen to a specific genbank

    Args:
        gb (dict[str, SeqRecord]): Takes a genbank as an argument, SeqRecord

    Returns:
        list[str]: returns a list of fields to be displayed in a single column
    """

    features: dict[str, list[str]] = {}
    feature_list: list[str] = []

    for genbank_file in genbank_list:
        
        # Using read.text to find a match with the string.
        genbank_text: str = Path(genbank_file).read_text()
        genbank: SeqRecord = gb[genbank_text]
        
        reference: Reference = genbank.annotations["references"]

        if re.match("^LOCUS", genbank):
            locus = re.sub("^LOCUS\s*", "", genbank)
            print(locus)
        elif re.match("^ACCESSION", genbank):
            accession = re.sub("^ACCESSION\s*", "", genbank)
            print(accession)
        elif re.match("^REFERENCE", genbank):
            reference = re.sub("^REFERENCE\s*", "",
                               genbank, flags=re.MULTILINE)
            print(reference)
            for feature in features:
                featurename = re.search(r"^{5}(\S+)", feature).group(1)
                print(feature)
                if re.search("keywords", locus, re.IGNORECASE) or re.search("keywords", reference, re.IGNORECASE) or re.search("keywords", feature, re.IGNORECASE):
                    print(accession)

    # Takes a list of items and unites them to be showed in a table
    feature_list.append(f"{locus}, {accession}, {reference}".strip())
    return feature_list


if __name__ == "__main__":

    pass
    # genbanks = read_list_genbanks(genbank_list)

    # c = get_features(genbanks)
    # print(c)
