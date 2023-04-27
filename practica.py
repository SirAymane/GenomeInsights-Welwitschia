
from Bio.SeqIO import SeqRecord
from pandas import DataFrame
from Bio.Align import Alignment

from scripts.Marc.get_gb import species, get_gb_by_accession
from scripts.Marc.regexps import get_apa_citation
from scripts.Aymane.parse_genbank import genbank_list, read_list_genbanks
from scripts.Melania.genes_amino_alignments import get_translation, align_aminos
from scripts.Siddique.alignment_final import get_gene_sequence, gene_seq, insert_gene_seq, align_genes

from pprint import pp


def get_main_table(genbanks: dict[str,SeqRecord], print_table: bool=True, export_excel: bool = True) -> DataFrame:

    """
    Wrapps up and returns the main findings of the genbanks introduced.
    For every genbank in dict accession number, descriptive name an rellevant
    regexp information are presented.

    Args:
        genbanks (dict[str,SeqRecord]): dictionary of SeqRecords where key is the genbank name id
        print_table (bool) optional: want to print the table? True by default
        export_excel (bool) optional: want to export as excel file? True by default
    
    """

    #Structure
    main_table: dict[str,list[any]]= {
        "genbank": ["welwitschia_mirabilis.gb","ephedra_trifurca.gb","gnetum_montanum.gb"],
        "Accession number": [],
        "Descriptive name": [],
        "Regexp APA": []
    }

    #Get accession numbers and names
    accession_numbers: list[str] = []
    names: list[str] = []
    for key in genbanks:
        gb: SeqRecord = genbanks[key]
        accession_number: str = gb.name
        name: str = gb.description

        accession_numbers.append(accession_number)
        names.append(name)

    # Add to main table
    main_table["Accession number"] = accession_numbers
    main_table["Descriptive name"] = names

    #Get APA citations by Regexp
    apa_citations: dict[str,list[str]] = get_apa_citation(genbanks)
    main_table["Regexp APA"] = apa_citations.values()

    #Do main table as DataFrame
    main_table_df: DataFrame = DataFrame.from_dict(main_table)

    if print_table:
        print(main_table_df)

    if export_excel:
        main_table_df.to_excel("main_table.xlsx", index=None)

    return main_table_df


def get_alignment_table(aligments: list[Alignment], print_table: bool=True, export_excel: bool = True, filename = "aligmnet_table.xlsx") -> DataFrame:

    """
    Returns a table of the alignments sorted by max puntuation

    Args:
        aligment: list of Aligment

    Returns:
        alignment_table (DataFrame)
    
    """

    #Structure
    main_table: dict[str,list[any]]= {
        "Alignment": [],
        "Score": []
    }

    aligment_n: int = 1
    for aligment in aligments:
        print(aligment)
        score: int = aligment.score

        main_table["Alignment"].append(aligment_n)
        main_table["Score"].append(score)

        aligment_n += 1

    main_table_df: DataFrame = DataFrame.from_dict(main_table)
    main_table_df = main_table_df.sort_values(by=["Score"], ascending=False)


    if print_table:
        print(main_table_df)

    if export_excel:
        main_table_df.to_excel(filename, index=None)

    return main_table_df




if __name__ == "__main__":

    """
    (1) Introducción: hipótesis
    
    Existe un tipo de planta curiosa que es considerada un fósil viviente dada su longevidad. 
    Welwitschia Mirabilis es una especie monotípica dentro de la família Welwitschia (única en su género). 
    Se estima que algunos especimenes llegan a superar los 1.000 años.
    
    El siguiente trabajo tiene por objetivo analizar el genóma de esta planta. Para su análisis, se han 
    obtenido los genómas de 2 especies cercanas más: Ephedra Trifurca y Gnetum Montanum.
    Todas ellas forman parte del orden de las Gnetales.

    A lo largo del trabajo, se han obtenido sus genómas y se han parseado y analizando biopython. 
    El objetivo es destacar diferencias y similitudes entre Welwitschia Mirabilis y los otros dos especimenes.

    Básicamente, el experimento planteado es el siguiente:
    disponemos de 3 genomas completos de la familia Welwitschia, donde el género
    Marabilis es único en su forma y composición. Vamos a realizar un primer estado de la cuestión
    sobre su bibliografía, sus datos y, además, vamos a comparar ciertos genes coincidentes
    para ver como difieren en su alineamiento. 

    La hipótesis es que debido a sus cualidades de supervivencia y su exposición a climas extremos, 
    la especia Welwitschia Mirabilis podria presentar diferencias en algunos de sus genes respecto a otras plantas.


    ###############################################
    
    (2) Descarregar Genbanks utilitzant Bio.Entrez

    En esta parte se pasa a obtener los genbanks de las 3 especies. La lógica de obtención es la siguiente:
    (1) se crea una lista de especies, 
    (2) se mira si dicho gb ya existe, 
    (3) si no existe se realiza la consulta y se guarda un fichero gb en formato texto.
    """

    get_gb_by_accession(species)

    """
    Una vez obtenidos los ficheros, los cargamos en memoria para trabajarlos posteriormente
    La variable 'genbanks' resultante es un diccionario, donde cada llave es el nombre
    del fichero genbank y su contenido el genbank en clase SeqRecord
    """

    genbanks: dict[str,SeqRecord] = read_list_genbanks(genbank_list)
    #print(genbanks)

    ###############################################

    """
    (3) Expresiones regulares
    En este apartado se procede a extraer información desestructurada
    y, habitualmente, variable mediante expresiones regulares.
    El ejercicio consiste en, a través de expresiones regulares, construir las citas
    APA de las referencias de los genbanks
    """

    apa_citations: dict[str,list[str]] = get_apa_citation(genbanks)
    #print(apa_citations["welwitschia_mirabilis.gb"][0])

    


    """
    (4) Alineamientos
    """
    ## Alineament d'aminoàcids
    genbanks: dict[str,SeqRecord] = read_list_genbanks(genbank_list) #welwitschia_mirabilis.gb, ephedra_trifurca.gb, gnetum_montanum.gb
    gen1: str = "matK" # Codifica para una proteína implicada en el fotosistema II de las plantas. La proteína MatK es un componente estructural de la subunidad de reacción del fotosistema II, que es responsable de la transferencia de electrones durante la fotosíntesis. Además de su papel en la fotosíntesis, matK se ha utilizado como un marcador molecular útil para estudios filogenéticos de plantas debido a su conservación a través de una amplia gama de especies y su relativamente alta tasa de evolución. 
    gen2: str = "rbcL" # Codifica la proteína de la subunidad de la cianobacteriana Rubisco, que es la enzima clave en la fotosíntesis.
    gen3: str = "atpF" # Codifica una subunidad de la ATP sintasa, una enzima que produce ATP a partir de la energía química liberada en la fotosíntesis.
    gen4: str = "psaB" # El gen psaB codifica una subunidad proteica que es un componente esencial del Fotosistema I en plantas y algunos otros organismos.

    # Obtenir les traduccions de cada gen dels tres genbanks
    trans_gen1 = get_translation(gen1, genbanks)
    trans_gen2 = get_translation(gen2, genbanks)
    trans_gen3 = get_translation(gen3, genbanks)
    trans_gen4 = get_translation(gen4, genbanks)

    # Alineaments primer gen
    # Alineamiento matK 1: welwitschia_mirabilis - ephedra trifurca
    al1 = align_aminos(trans_gen1, 0, 1)
    # Alineamiento matK 2: welwitschia_mirabilis - gnetum_montanum
    al2 = align_aminos(trans_gen1, 0, 2)
    # Alineamiento matK 3: ephedra trifurca - gnetum_montanum
    al3 = align_aminos(trans_gen1, 1, 2)

    # Alineaments segon gen
    # Alineamiento rbcL 1: welwitschia_mirabilis - ephedra trifurca
    al4 = align_aminos(trans_gen2, 0, 1)
    # Alineamiento rbcL 2: welwitschia_mirabilis - gnetum_montanum
    al5 = align_aminos(trans_gen2, 0, 2)
    # Alineamiento rbcL 3: ephedra trifurca - gnetum_montanum
    al6 = align_aminos(trans_gen2, 1, 2)

    # Alineaments tercer gen
    # Alineamiento atpF 1: welwitschia_mirabilis - ephedra trifurca
    al7 = align_aminos(trans_gen3, 0, 1)
    # Alineamiento atpF 2: welwitschia_mirabilis - gnetum_montanum
    al8 = align_aminos(trans_gen3, 0, 2)
     # Alineamiento atpF 3: ephedra trifurca - gnetum_montanum
    al9 = align_aminos(trans_gen3, 1, 2)

    # Alineaments quart gen
    # Alineamiento psaB 1: welwitschia_mirabilis - ephedra trifurca
    al10 = align_aminos(trans_gen4, 0, 1)
    # Alineamiento psaB 2: welwitschia_mirabilis - gnetum_montanum
    al11 = align_aminos(trans_gen4, 0, 2)
    # Alineamiento psaB 3: ephedra trifurca - gnetum_montanum
    al12 = align_aminos(trans_gen4, 1, 2)

    alignment_list: list[Alignment] = [
        al1,
        al2,
        al3,
        al4,
        al5,
        al6,
        al7,
        al8,
        al9,
        al10,
        al11,
        al12
    ]


    # Alineament nucleotids
    ##################### Get seqs in list ####################

    # Genes
    gene1 = "matK"
    gene2 = "rbcL"
    gene3 = "atpF"
    gene4 = "psaB"

    # Get gene seq of matK
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
    ali1 = align_genes(gene1_seq_list[MIRABILIS], gene1_seq_list[TRIFURCA])
    ali2 = align_genes(gene1_seq_list[MIRABILIS], gene1_seq_list[MONTANUM])


    # rbcL
    ali3 = align_genes(gene2_seq_list[MIRABILIS], gene2_seq_list[TRIFURCA])
    ali4 = align_genes(gene2_seq_list[MIRABILIS], gene2_seq_list[MONTANUM])

    # # atpF
    ali5 = align_genes(gene3_seq_list[MIRABILIS], gene3_seq_list[TRIFURCA])
    ali6 = align_genes(gene3_seq_list[MIRABILIS], gene3_seq_list[MONTANUM])

    # psaB
    ali7 = align_genes(gene4_seq_list[MIRABILIS], gene4_seq_list[TRIFURCA])
    ali8 = align_genes(gene4_seq_list[MIRABILIS], gene4_seq_list[MONTANUM])


    alignment_list2: list[Alignment] = [
        ali1,
        ali2,
        ali3,
        ali4,
        ali5,
        ali6,
        ali7,
        ali8,
    ]

    """
    (5) Resultados
    """
    main_table: DataFrame = get_main_table(genbanks)
    aligment_table: DataFrame = get_alignment_table(alignment_list)
    aligment_table2: DataFrame = get_alignment_table(alignment_list2, filename="alignment_table2.xlsx")

    """
    Los resultados muestran que la hipótesis de que la planta welwitschia_mirabilis podría presentar
    diferencias en alguno de sus genes principales de su cloroplastos respecto a otras especies debido 
    a sus cualidades de sobrevivir en un ambiente extremadamente hostil durante cientos de años no se 
    ha cumplido ya que observamos mínimas diferencias en cualquiera de los alineamientos.

    Para llegar a esta conclusión, es especialmente significativo que en los alineamientos con una 
    puntuación más baja, también se da una puntuación igualmente baja entre las dos especies de control, 
    lo que descartaría una diferencia significativa de la planta estudiada.
    """