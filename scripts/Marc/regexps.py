import re 
from pathlib import Path 
from Bio.SeqIO import SeqRecord
from Bio.SeqFeature import Reference
from pprint import pprint
#from parse_genbank import genbank_list, read_list_genbanks
from pprint import pp


def get_apa_citation(gb: dict[str,SeqRecord]) -> dict[str,list[str]]:
    """
    Make the APA citation with available information
    in the genbank

    Args:
        gb (dict[str,SeqRecord]) a dictionary of genbanks

    Returns:
        citations (dict[str,list[str]]):    a dictionary of lists of strings (APA citations),
                                            where keys are the genbank name
    """

    # Initialize result
    citations: dict[str,list[str]] = {}

    # Date and Journal regexp pattern (2 groups)
    # Sometimes, the pattern follows date+journal (straight),
    # Sometimes, the pattern follows journal+date (reverse)
    # If straight is empty, reverse are applied 
    straight_date_and_journal_regexp: str = r"\(([0-3][0-9]-[A-Z]{3}-[0-9]{4}|[0-9]{4})\) (.*)?"
    reverse_date_and_journal_regexp: str = r"(.*) \(([0-3][0-9]-[A-Z]{3}-[0-9]{4}|[0-9]{4})\)"

    for specie in gb.keys():
        genbank: SeqRecord = gb[specie]
        references: Reference = genbank.annotations["references"]

        cites: list[str] = []
        for reference in references:

            date: str = ""
            journal: str = ""
            
            authors = reference.authors

            matches = re.findall(straight_date_and_journal_regexp,reference.journal)
            if len(matches) == 0:
                matches = re.findall(reverse_date_and_journal_regexp,reference.journal)
                date = matches[0][1]
                journal = matches[0][0]
            else:
                date = matches[0][0]
                journal = matches[0][1]

            cites.append(f"{authors} {date}, {reference.title}, {journal}".strip())

        citations[specie] = cites

    return citations


if __name__ == "__main__":
    pass

    # genbanks = read_list_genbanks(genbank_list)

    # c = get_apa_citation(genbanks)
    # print(c["ephedra_trifurca.gb"][0])
