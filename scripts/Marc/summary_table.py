from pandas import DataFrame





if __name__ == "__name__":

    main_table: dict[str,list[any]]= {
        "genbank": ["welwitschia_mirabilis.gb","ephedra_trifurca.gb","gnetum_montanum.gb"],
        "Accession number": [],
        "Descriptive name": [],
        "Regexp APA": []
    }

    alignment_table: dict[str,list[any]] = {
        "Reference genbank": [],
        "Aligned to": [],
        "Alignment type": [],
        "Gene": [],
        "Score": []
    }


