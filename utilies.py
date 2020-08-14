

# imports
from string import ascii_letters
from utilies import MultiGeneBlastException

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

def testaccession(accession):
    #TODO look into this function at some time
    #Test if accession number is probably real GenBank/RefSeq acc nr
    numbers = list(range(0,10))
    letters = []
    for i in ascii_letters:
        letters.append(i)
    nrnumbers = 0
    nrletters = 0
    for i in accession:
        if i in letters:
            nrletters += 1
        try:
            j = int(i)
            if j in numbers:
                nrnumbers += 1
        except:
            pass
    return not (nrnumbers < 3 or nrletters < 1)

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
    #TODO take a better look at this function, for now it seems fine :)
    #Translation table standard genetic code; according to http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
    transldict = { 'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',
                   'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C',
                   'TTA': 'L', 'TCA': 'S', 'TAA': '*', 'TGA': '*',
                   'TTG': 'L', 'TCG': 'S', 'TAG': '*', 'TGG': 'W',
                   'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R',
                   'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',
                   'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',
                   'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',
                   'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S',
                   'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',
                   'ATA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',
                   'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',
                   'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G',
                   'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',
                   'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',
                   'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G',
                   'ttt': 'F', 'tct': 'S', 'tat': 'Y', 'tgt': 'C',
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
    triplets = []
    triplet = ""
    a = 0
    for i in sequence:
        if a < 2:
            a += 1
            triplet = triplet + i
        elif a == 2:
            triplet = triplet + i
            triplets.append(triplet)
            triplet = ""
            a = 0
    protseq = ""
    aanr = 0
    for i in triplets:
        aanr += 1
        if aanr == 1:
            protseq = protseq + "M"
        else:
            if "n" in i or "N" in i or i not in list(transldict.keys()):
                protseq = protseq + "X"
            else:
                protseq = protseq + transldict[i]
    if  len(protseq) > 0 and protseq[-1] == "*":
        protseq = protseq[:-1]
    return protseq