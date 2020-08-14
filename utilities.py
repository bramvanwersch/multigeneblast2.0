

# imports
from string import ascii_letters
import logging
import time
import os, sys
import datetime

ILLEGAL_CHARACTERS = ["'",'"','=',';',':','[',']','>','<','|','\\',"/",'*','-','_','.',',','?',')','(','^','#','!','`','~','+','{','}','@','$','%','&']


def setup_logger(outdir, starttime, level=logging.DEBUG):
    """
    Function for setting up a logger that will write output to file as well as
    to sdout.

    :param outdir: output directory specified by the user. this is where the
    log file should end up.
    :param starttime: the time the program was started. The logger is created
    slightly later.
    :param level: the logging level, default is DEBUG
    """

    #if the directory exists simply ignore it, that can be expected
    dir_exists = False
    try:
        os.mkdir(outdir)
    except FileExistsError:
        dir_exists = True

    log_file_loc = "{}{}{}".format(outdir, os.sep, 'run.log')
    #make sure that a potentially existing logfile is emptied
    if os.path.exists(log_file_loc):
        open(log_file_loc, "w").close()

    #TODO think about making the level a user definable parameter
    logging.basicConfig(filename=log_file_loc, level=level, format='%(levelname)s: %(asctime)s - %(message)s')

    #configure a handler that formats the logged events properly and prints the events to file as well as stdout
    handler = logging.StreamHandler(sys.stdout)
    formatter = MyFormatter('%(currentTime)s (%(passedTime)s sec) - %(message)s', starttime=starttime)
    handler.setFormatter(formatter)
    logging.getLogger().addHandler(handler)

    logging.debug('Logger created')

    if dir_exists:
        logging.warning("The output directory already exists. Files may be overwritten.")


class MyFormatter(logging.Formatter):
    """
    Formatter subclass that saves the start time so the time can be displayed
    since starting to run MultiGeneBlast
    """
    def __init__(self, fmt, starttime=time.time()):
        logging.Formatter.__init__(self, fmt)
        self._start_time = starttime

    def format(self, record):
        """
        Overwrite of the format function that prints the passed time and adds
        current time to the existing format

        :See: logging.Formatter.format()
        """
        #difference = datetime.datetime.now() - self._start_time
        record.passedTime = "{:.3f}".format(time.time() - self._start_time)
        record.currentTime = datetime.datetime.now().time()
        return super(MyFormatter, self).format(record)


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
