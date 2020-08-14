

# imports
from string import ascii_letters

ILLEGAL_CHARACTERS = ["'",'"','=',';',':','[',']','>','<','|','\\',"/",'*','-','_','.',',','?',')','(','^','#','!','`','~','+','{','}','@','$','%','&']


class MultiGeneBlastException(Exception):
    """
    Create MultiGeneBlastExceptions for error that are expected from MultiGeneBlast
    """
    pass


def remove_illegal_characters(string):
    """
    Remove any character from ILLEGAL_CHARACTERS from string

    :param string: a string
    :return: the string without illegal characters. If no illegal characters
    are found the originall string is returned
    """
    clean_string = ""
    for letter in string:
        if letter not in ILLEGAL_CHARACTERS:
            clean_string += letter
    return clean_string

def is_valid_accession(accession):
    """
    Test if accession number is probably real GenBank/RefSeq acc nr

    :param accession: a string that is a probable Genbank/Refseq accession
    :return: a boolean
    """
    nrnumbers = 0
    nrletters = 0
    for value in accession:
        if value.isdigit():
            nrnumbers += 1
        elif value.isalpha():
            nrletters += 1
    return nrnumbers >= 3 or nrletters >= 1

def complement(seq):
    """
    Create the complement strand of a sequence. Use replacing for speed

    :param seq: a DNA sequence as a string
    :return: the complement of that string
    """
    seq = seq.replace("a", "x").replace("A", "X")
    seq = seq.replace("t", "a").replace("T", "A")
    seq = seq.replace("x", "t").replace("X", "T")

    seq = seq.replace("c", "x").replace("C", "X")
    seq = seq.replace("g", "c").replace("G", "C")
    seq = seq.replace("x", "g").replace("X", "G")
    return seq

def translate(sequence):
    """
    Translate a nucleotide sequence into an amino acid sequence

    :param sequence: a string that is a nucleotide sequence
    :return: a string that is an amino acid sequence
    """
    #Translation table standard genetic code; according to http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
    transldict = { 'ttt': 'F', 'tct': 'S', 'tat': 'Y', 'tgt': 'C',
                   'ttc': 'F', 'tcc': 'S', 'tac': 'Y', 'tgc': 'C',
                   'tta': 'L', 'tca': 'S', 'taa': '*', 'tga': '*',
                   'ttg': 'L', 'tcg': 'S', 'tag': '*', 'tgg': 'W',
                   'ctt': 'L', 'cct': 'P', 'cat': 'H', 'cgt': 'R',
                   'ctc': 'L', 'ccc': 'P', 'cac': 'H', 'cgc': 'R',
                   'cta': 'L', 'cca': 'P', 'caa': 'Q', 'cga': 'R',
                   'ctg': 'L', 'ccg': 'P', 'cag': 'Q', 'cgg': 'R',
                   'att': 'I', 'act': 'T', 'aat': 'N', 'agt': 'S',
                   'atc': 'I', 'acc': 'T', 'aac': 'N', 'agc': 'S',
                   'ata': 'I', 'aca': 'T', 'aaa': 'K', 'aga': 'R',
                   'atg': 'M', 'acg': 'T', 'aag': 'K', 'agg': 'R',
                   'gtt': 'V', 'gct': 'A', 'gat': 'D', 'ggt': 'G',
                   'gtc': 'V', 'gcc': 'A', 'gac': 'D', 'ggc': 'G',
                   'gta': 'V', 'gca': 'A', 'gaa': 'E', 'gga': 'G',
                   'gtg': 'V', 'gcg': 'A', 'gag': 'E', 'ggg': 'G'}
    protseq = ""
    for triplet_start in range(0, len(sequence), 3):
        triplet = sequence[triplet_start : triplet_start + 3].lower()
        if "n" in triplet or triplet not in transldict:
            protseq += "X"
        else:
            protseq += transldict[triplet]
    if  len(protseq) > 0 and protseq[-1] == "*":
        protseq = protseq[:-1]
    return protseq
