

# imports
from string import ascii_letters
import logging
import subprocess
import time
import os, sys
import datetime
import multiprocessing

from constants import TEMP

ILLEGAL_CHARACTERS = ["'",'"','=',';',':','[',']','>','<','|','\\',"/",'*','-','.',',','?',')','(','^','#','!','`','~','+','{','}','@','$','%','&']


def setup_temp_folder():
    """
    Setup a the working directory for MultiGeneBlast. This directory is emptied
    when setting it up ensuring no confilicting files are present when running
    any of the scripts
    """
    mgb_temp_folder = TEMP
    try:
        shutil.rmtree(mgb_temp_folder)
    except:
        pass
    try:
        os.mkdir(mgb_temp_folder)
    except:
        pass
    os.chdir(mgb_temp_folder)

def setup_logger(outdir, starttime, level="basic"):
    """
    Function for setting up a logger that will write output to file as well as
    to sdout.
    :param outdir: output directory specified by the user. this is where the
    log file should end up.
    :param starttime: the time the program was started. The logger is created
    slightly later.
    :param level: the logging level, default is basic, meaning all logging.info
    messages
    """

    if level == "basic":
        level = logging.INFO
    elif level == "all":
        level = logging.DEBUG
    else:
        level = logging.WARNING

    #if the directory exists simply ignore it, that can be expected
    dir_exists = False
    try:
        os.mkdir(outdir)
    except FileExistsError:
        pass

    log_file_loc = "{}{}{}".format(outdir, os.sep, 'run.log')
    #make sure that a potentially existing logfile is emptied
    if os.path.exists(log_file_loc):
        open(log_file_loc, "w").close()

    logging.basicConfig(filename=log_file_loc, level=logging.DEBUG,
                        format='%(levelname)s: %(asctime)s - %(message)s')

    # configure a handler that formats the logged events properly and prints the events to file as well as stdout
    handler = logging.StreamHandler(sys.stdout)
    formatter = MyFormatter('%(levelname)s %(currentTime)s (%(passedTime)s sec) - %(message)s', starttime=starttime)
    handler.setFormatter(formatter)
    handler.setLevel(level)
    logging.getLogger().addHandler(handler)

    logging.debug('Logger created')


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

def run_commandline_command(command, max_retries = 5):
    """
    Run a command line command that can be repeadet when a error is returned.
    This function is meant to run the BLAST+ command line tools

    :param command: a string that can be deployed on the command line as a
    command
    :param max_retries: The maximum amount of times the program should retry the
    command when an error is returned
    :raises MultiGeneBlastError: when the command returned max_retries amount of
    errors
    """
    new_env = os.environ.copy()
    command_stdout = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=new_env)
    command_stdout = command_stdout.stdout.read()
    retries = 0
    while "error" in str(command_stdout.lower()):
        logging.debug("The following command {} returned the following error {}. Retrying: {}/{}".format(command, command_stdout, retries, max_retries))
        if max_retries <= retries:
            logging.critical("Command {} keeps returing an error. Exiting...".format(command))
            raise MultiGeneBlastException()
        command_stdout = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=new_env)
        command_stdout = command_stdout.stdout.read()
        retries += 1


class MultiGeneBlastException(Exception):
    """
    Create MultiGeneBlastExceptions for error that are expected from MultiGeneBlast
    """
    pass


def determine_cpu_nr(cores):
    """
    Determine the number of CPU's needed based on the nr provided by the user.

    :param cores: a string that represents the amount of cores
    :return: an integer that is between 1 and maximum cores

    the cores can be any integer or 'all'. If 'all' is provided the maximum
    amount of cores is selected. If a number higher then the maximum number of
    cores is requested the maximum number is returned.
    """
    if cores.lower() == "all":
        try:
            nrcpus = multiprocessing.cpu_count()
        except(IOError,OSError,NotImplementedError):
            nrcpus = 1
    else:
        cores = int(cores)
        try:
            nrcpus = multiprocessing.cpu_count()
        except(IOError,OSError,NotImplementedError):
            nrcpus = 1
        if cores < nrcpus:
            nrcpus = cores
    return nrcpus

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

def fasta_to_dict(file_name, check_headers = True):
    """
    Creates a dictionary containing the name of the fasta sequence as key
    and the sequence as value.

    :param file_name: A string that represents a file path towards a fasta
    file containing one or more sequences
    :param check_headers: Boolean if the headers should be filtered for illegal
    characters. Default is True
    :return: a dictionary with sequence names as keys and the sequence
    itself as values.
    """
    try:
        with open(file_name, "r") as f:
            text = f.read()
    except Exception:
        logging.critical("Could not read file {}. Exiting...".format(file_name))
        raise MultiGeneBlastException("The file {} does not exist anymore or cannot be read.".format(file))
    sequences = {}
    # do not include the first empty match that results from the split
    fasta_entries = text.split(">")[1:]
    if len(fasta_entries) == 0:
        logging.critical("Invalid fasta format for file '{}'. Exiting...".format(file_name))
        raise MultiGeneBlastException("Fasta file '{}' does not contain any sequences.".format(file_name))
    for entry in fasta_entries:
        lines = entry.split("\n")

        #make sure that names do not contain illegal characters. This makes sure no strange results are retrieved from
        #blast+ tools
        name = remove_illegal_characters(lines[0])
        #make sure no trailing newlines
        sequence = "".join(lines[1:]).strip()
        #skip incomplete entries
        if len(sequence) == 0:
            logging.warning("Invalid fasta format for entry '{}' in file '{}'. Skipping...".format(name, file_name))
        elif name in sequences:
            logging.warning("Double fasta entry '{}' in file '{}'. Skipping...".format(name, file_name))
        else:
            sequences[name] = sequence
    return sequences

def is_dna(sequence):
    """
    Determine is a sequence is most likely a DNA strand or not.
    Note is unreliable for sequences that are to short < 10bp

    :param sequence: a sequence as string
    :return: a boolean
    """
    sequence = sequence.lower()
    if sequence.count("a") + sequence.count("g") + sequence.count("t") + sequence.count("c") < 0.5 * len(sequence):
        return False
    return True