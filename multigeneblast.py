#!/usr/bin/env python
## Copyright (c) 2012 Marnix H. Medema
## Department of Microbial Physiology / Groningen Bioinformatics Centre
## University of Groningen
## License: GNU General Public License v3 or later
## A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

##Imported modules

import os
from os import system
import sys
import datetime

sys.path.append('.\\pysvg')

import time
import multiprocessing
from multiprocessing import Process, freeze_support
import random
import fileinput
import subprocess
import argparse
import logging
from collections import OrderedDict

#constants
FASTA_EXTENSIONS = ("fasta","fas","fa","fna")
EMBL_EXTENSIONS = ("embl","emb")
GENBANK_EXTENSIONS = ("gbk","gb","genbank")
HITS_PER_PAGE = 50
SVG_CORE_EXTENSION = 5000

ILLEGAL_CHARACTERS = ["'",'"','=',';',':','[',']','>','<','|','\\',"/",'*','-','_','.',',','?',')','(','^','#','!','`','~','+','{','}','@','$','%','&']
global GUI
global OUTBOX
global FRAME
global CURRENTDIR
global MGBPATH
global APPDATA
global TEMP
global DBPATH


#Find path to mgb files if run from another directory
pathfolders = os.environ['PATH'].split(os.pathsep)
pathfolders.reverse()
pathfolders.append(os.getcwd())
pathfolders.reverse()
CURRENTDIR = os.getcwd()
MGBPATH = ""
for folder in pathfolders:
    try:
        if "read_input_gui.py" in os.listdir(folder) and "guilib.py" in os.listdir(folder) and "empty.xhtml" in os.listdir(folder) and "multigeneblast.py" in os.listdir(folder) and "mgb_gui.py" in os.listdir(folder):
            MGBPATH = folder
            break
    except:
        pass
try:
    if  MGBPATH == "" and os.sep in sys.argv[0] and "read_input_gui.py" in os.listdir(sys.argv[0].rpartition(os.sep)[0]) and "guilib.py" in os.listdir(sys.argv[0].rpartition(os.sep)[0]):
        MGBPATH = sys.argv[0].rpartition(os.sep)[0]
        os.chdir(MGBPATH)
except:
    pass
if MGBPATH == "":
    print("Error: Please add the MultiGeneBlast installation directory to your $PATH environment variable before running the executable from another folder.")
    sys.exit(1)
#Find path to Application Data
if sys.platform == ('win32'):
    APPDATA = os.environ['ALLUSERSPROFILE'] + os.sep + 'Application Data'
elif sys.platform == ('darwin'):
    APPDATA = os.path.expanduser("~") + "/Library/Application Support"
else:
    try:
        if os.path.exists(os.getcwd() + os.sep + "multigeneblast_data"):
            APPDATA = os.getcwd() + os.sep + "multigeneblast_data"
        else:
            os.mkdir(os.getcwd() + os.sep + "multigeneblast_data")
            APPDATA = os.getcwd() + os.sep + "multigeneblast_data"
    except:
        try:
            if os.path.exists(os.environ['HOME'] + os.sep + "multigeneblast_data"):
                APPDATA = os.getcwd() + os.sep + "multigeneblast_data"
            else:
                os.mkdir(os.environ['HOME'] + os.sep + "multigeneblast_data")
                APPDATA = os.environ['HOME'] + os.sep + "multigeneblast_data"
        except:
            print("No permission to write to installation folder. Please change user or save somewhere else.")
            sys.exit()
if sys.platform == ('darwin') or sys.platform == ('win32'):
    try:
        os.mkdir(APPDATA + os.sep + 'MultiGeneBlast')
        APPDATA = APPDATA + os.sep + 'MultiGeneBlast'
    except:
        if os.path.exists(APPDATA + os.sep + 'MultiGeneBlast'):
            APPDATA = APPDATA + os.sep + 'MultiGeneBlast'
#Find path to temporary files
if sys.platform == ('win32'):
    TEMP = os.environ['TEMP']
elif sys.platform == ('darwin'):
    TEMP = os.environ['TMPDIR']
else:
    try:
        os.mkdir(os.environ['HOME'] + os.sep + ".mgbtemp")
        TEMP = os.environ['HOME'] + os.sep + ".mgbtemp"
    except:
        TEMP = APPDATA
#Set other environment variables
os.environ['EXEC'] = MGBPATH + os.sep + "exec"
os.environ['PATH'] = APPDATA + os.pathsep + os.environ['EXEC'] + os.pathsep + os.environ['PATH']


from pysvg.filter import *
from pysvg.gradient import *
from pysvg.linking import *
from pysvg.script import *
from pysvg.shape import *
from pysvg.structure import *
from pysvg.style import *
from pysvg.text import *
from pysvg.builders import *
from string import ascii_letters
import urllib.request, urllib.error, urllib.parse
from urllib.request import Request, urlopen
from urllib.error import URLError, HTTPError
import http.client
from http.client import BadStatusLine,HTTPException
import urllib.request, urllib.parse, urllib.error
import tarfile
import pickle as pickle
from tkinter import *
from tkinter.messagebox import askyesno, showerror
import shutil

from genbank_parsing import GenbankFile

### GENERAL UTILITY###

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

def write_fasta(protein_list, file):
    """
    Write a fasta file using a list of Protein objects to the Temp directory of
    the operating system

    :param protein_list: a list of Protein objects
    :param file: the output file
    """
    #the working directory is set to a temporary directory
    try:
       with open(file,"w") as out_file:
            for prot in protein_list:
                out_file.write(prot.fasta_text())
    except Exception:
        logging.critical("No fasta file created for sequences with file name {}".format(file))
        raise MultiGeneBlastException("Cannot open file {}".format(file))
    logging.debug("Saved file {} at {}.".format(file, os.getcwd()))

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

        #make sure that names do not contain illegal characters.
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

####PARSE OPTIONS####

class Options:
    """
    Options object for saving all the user defined options aswell as some
    variables that can be directly calculated from these options like the
    architecture_mode value
    """
    def __init__(self, arguments):
        """
        :param arguments: an Argparser object containing all the options
        specified by the user.
        """
        #the input file
        self.infile = arguments.i
        #the output directory
        self.outdir = arguments.o
        #the database to use
        #TODO make sure this gets configured based on the database input
        self.db_mode = "local"
        self.dbtype, self.db = self.__configure_db_and_type(arguments.db)

        self.architecture_mode = self.__check_architecture_mode()
        #values that define the query
        self.startpos = arguments.f
        self.endpos = arguments.t
        self.ingenes = arguments.g

        self.cores = arguments.c
        self.minseqcov = arguments.msc
        self.minpercid = arguments.mpi
        self.hitspergene = arguments.hpg
        #TODO make sure that the *2 is correct
        self.distancekb = int(arguments.dkb * 1000)
        self.muscle = arguments.m
        self.pages = arguments.op
        self.syntenyweight = arguments.sw

        #TODO move this value to a more logical place
        self.screenwidth = 1024
        self.gui = False

    def __check_architecture_mode(self):
        """
        Check if MultiGeneBlast is running in architecture or homolgy search mode

        :return: Boolean that is True if the infile ends in a fasta extension
        meaning MultiGeneBlast is running in architecture mode.
        """
        if any(self.infile.lower().endswith(ext) for ext in FASTA_EXTENSIONS):
            return True
        return False

    def __configure_db_and_type(self, database_path):
        to_path, db_file = os.path.split(database_path)
        #very important to configure. Otherwise BLAST+ cannot find database files
        os.environ['BLASTDB'] = to_path

        # make sure the databae file is of the correct file type
        dbname, ext = db_file.split(".")
        if ext == "pal":
            return "prot", dbname
        elif ext == "nal":
            return "nucl", dbname

def get_arguments():
    """
    Parse the command line arguments using argparser.

    :return: an Option object that holds the options specified by the user
    """
    parser = argparse.ArgumentParser(description='Run multigeneblast on a '
                                'specified database using a specified query.',
                                     epilog="-from, -to and -genes are "
                                            "mutually exclusive")
    parser.add_argument("-i", "-in", help="Query file name: GBK/EMBL file for "
                                "homology search,FASTA file with multiple"
                                " protein sequences for architecture search",
                        required=True, type=check_in_file, metavar="file path")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-f", "-from", help="Start position of query region", type=int
                       , metavar="Int")
    #nargs allows one or more arguments
    group.add_argument("-g", "-genes", help="Accession codes of genes constituting "
                                      "query multigene module", nargs='+'
                       , metavar="coma-separted gene identifiers")

    parser.add_argument("-t", "-to", help="End position of query region", type=int,
                        metavar="Int")

    parser.add_argument("-o", "-out", help="Output folder path in which results will"
                                     " be stored", type=check_out_folder, metavar="file path")
    parser.add_argument("-db", "-database", help="Blast database to be queried",
                        required=True, type=check_db_folder, metavar="file path")

    parser.add_argument("-c", "-cores", help="Number of parallel CPUs to use for "
                                       "threading (default: all)",
                        type=determine_cpu_nr, default="all", metavar="Int")
    parser.add_argument("-hpg", "-hitspergene", help="Number of Blast hits per query "
                                             "gene to be taken into account "
                                             "(default: 250), max = 10000",
                        type=int, choices=range(50,10001), default=250,
                        metavar="[50 - 10000]")
    parser.add_argument("-msc", "-minseqcov", help="Minimal %% coverage of a Blast hit"
                                           " on hit protein sequence to be "
                                           "taken into account (default: 25)",
                        type=int, choices=range(0,101), default=25,
                        metavar="[0 - 100]")
    parser.add_argument("-mpi", "-minpercid", help="Minimal %% identity of a Blast hit"
                                           " on hit protein sequence to be "
                                           "taken into account (default: 30)",
                        type=int, choices=range(0,101), default=30,
                        metavar="[0 - 100]")
    parser.add_argument("-dkb", "-distancekb", help="Maximum kb distance between two"
                                            " blast hits to be counted as"
                                            " belonging to the same locus"
                                            " (default: 10)", default=20,
                        type=int, choices=range(0,100),
                        metavar="[0 - 100]")
    parser.add_argument("-sw", "-syntenywheight", help="Weight of synteny conservation"
                                               " in hit sorting score: (default"
                                               ": 0.5)", type=restricted_float,
                        default=0.5, metavar="Float")
    parser.add_argument("-m", "-muscle", help="Generate a Muscle multiple sequence"
                                        " alignments of all hits of each input"
                                        " gene (default: n)", default="n",
                        choices=["y","n"])
    parser.add_argument("-op", "-outpages", help="Maximum number of output pages"
                                          " (with 50 hits each) to be generated"
                                          " (default: 5)", type=int, default=5,
                        choices=range(1,41), metavar="[1 - 40]")

    name_space = parser.parse_args()

    #some final checks for certain arguments, raise error when appropriate
    #when defining a from without to or a to without a from
    if (name_space.f is not None and name_space.t is None) or \
            (name_space.f is None and name_space.t is not  None):
        parser.error("When specifying -from also specify -to.")

    #if the -in argument is a fasta file make sure that no -from, -to or genes are defined
    #otherwise make sure that they are defined
    if (any(name_space.i.endswith(ext) for ext in  FASTA_EXTENSIONS) and
        (name_space.t is not None or name_space.f is not None or name_space.g is not None)):
        parser.error("When providing a fasta file (running architecture mode) -f, -t and"
                     " -g arguments cannot be specified")
    elif (not any(name_space.i.endswith(ext) for ext in  FASTA_EXTENSIONS) and
        (name_space.t is None and name_space.f is None and name_space.g is None)):
        parser.error("When providing a genbank or embl file -f, -t or -g should be specified")

    user_options = Options(name_space)
    return user_options

def check_in_file(path):
    """
    Check the provided input file

    :param path: a string that is an absolute or relative path
    :raises ArgumentTypeError: when the path provided does not exist or the 
    or the specified file has the wrong extension.
    :return: the original path
    """

    try:
        #make sure relatively defined paths are also correct
        my_path = os.path.abspath(os.path.dirname(__file__))
        path = os.path.join(my_path, path)

        assert os.path.exists(path)
        #the extension of the file
        ext = os.path.split(path)[1].split(".")[1]
        assert ext.lower() in [*GENBANK_EXTENSIONS, *EMBL_EXTENSIONS, *FASTA_EXTENSIONS]
        return path
    except (AssertionError, IndexError):
        raise argparse.ArgumentTypeError("Please supply input file with valid"
                                         " GBK / EMBL extension (homology "
                                         "search) or FASTA extension "
                                         "(architecture search).")

def check_out_folder(path):
    """
    Check the provided output folder

    :param path: a string that is an absolute or relative path
    :raises ArgumentTypeError: when the path provided does not exist or the
    specified name for the folder contains illegal characters
    :return: the original path
    """
    try:
        #make sure relatively defined paths are also correct
        my_path = os.path.abspath(os.path.dirname(__file__))
        path = os.path.join(my_path, path)

        to_path, folder_name = os.path.split(path)
        assert os.path.exists(to_path)

        #assert that the folder name does not contain illegal characters
        assert folder_name.replace("_", "").isalnum()
        return path
    except AssertionError:
        raise argparse.ArgumentTypeError("Output folder does exist or cannot be"
                                         " created.")

def check_db_folder(path):
    """
    Check the provided database folder to make sure that all the required files
    are present

    :param path: a string that is an absolute or relative path
    :raises ArgumentTypeError: when the path provided does not exist or the
    database folder does not contain all neccesairy files
    :return: the original path
    """
    try:
        #make sure relatively defined paths are also correct
        my_path = os.path.abspath(os.path.dirname(__file__))
        path = os.path.join(my_path, path)

        assert os.path.exists(path)
        to_path, db_file = os.path.split(path)

        #make sure the databae file is of the correct file type
        db_name, ext = db_file.split(".")
        assert ext in ("pal", "nal")
        #make sure the db folder contains all files required
        db_folder = os.listdir(to_path)
        #TODO make sure that not more files need to be checked
        if ext == "pal":
            assert "{}.{}".format(db_name, "phr") in db_folder
        else:
            assert "{}.{}".format(db_name, "nhr") in db_folder
        return path
    except (AssertionError, IndexError):
        raise argparse.ArgumentTypeError("The provided path is incorrect or "
                                "not all neccesairy data base files exist.")

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

def restricted_float(x):
    """
    Make sure the float for syntenywheight is inbetween 0.0 and 2.0

    :param x: a string that is a float
    :raises ArgumentTypeError: when the value is not a float or not inbetween 0
    and 2
    :return: the float representation of the string
    """
    try:
        x = float(x)
    except ValueError:
        raise argparse.ArgumentTypeError("{} not a floating-point literal".format(x))

    if x < 0.0 or x > 2.0:
        raise argparse.ArgumentTypeError("{} not in range [0.0, 2.0]".format(x))
    return x

####READ THE INPUT FILE###

def read_query_file(user_options):
    """
    Read the query file supplied by the user. This can be a fasta, genbank or
    embl file. A fasta version of the file is saved in query.fasta for blasting

    :param user_options: an Option object containing the user defined options
    :return: a dictionary containing keys for the fasta headers the protein
    sequences are saved with in the query.fasta file and Protein objects as
    values
    """
    logging.info("Starting to parse input file...")
    if user_options.architecture_mode:
        query_proteins = query_proteins_from_fasta(user_options.infile)
        write_fasta(query_proteins.values(), "query.fasta")
        return query_proteins
    elif any(user_options.infile.lower().endswith(ext) for ext in GENBANK_EXTENSIONS):
        gb_file = GenbankFile(user_options.infile, protein_range=[user_options.startpos, user_options.endpos], allowed_proteins=user_options.ingenes)
        query_proteins = gb_file.proteins
    else:
        genbank_file = embl_to_genbank(user_options.infile)
        gb_file = GenbankFile(user_options.infile, file_text=genbank_file, protein_range=[user_options.startpos, user_options.endpos], allowed_proteins=user_options.ingenes)
        query_proteins = gb_file.proteins

    if len(query_proteins) == 0:
        logging.critical("No proteins found in the provided query. Exiting...")
        raise MultiGeneBlastException("No proteins found in the provided query")

    write_fasta(query_proteins.values(), "query.fasta")
    logging.info("{} proteins have been extracted from the query.".format(len(query_proteins)))
    return query_proteins

def query_proteins_from_fasta(fasta_file):
    """
    Convert a fasta file containing query proteins for architecture search into
    a dictionary of Protein objects.

    :param fasta_file: a file path to a fasta file
    :return: a dictionary containing keys for the fasta headers the protein
    sequences are saved with in the query.fasta file and Protein objects as
    values
    """
    logging.debug("Started parsing architecture input file")

    #read the fasta file with headers that do not contain forbidden characters
    fasta_entries = fasta_to_dict(fasta_file)

    #make for each fasta sequence a protein object and save them in a dictionary
    total_lenght = 0
    query_proteins = {}
    for entry_name in fasta_entries:
        sequence = fasta_entries[entry_name]
        entry_protein = Protein(sequence, total_lenght, entry_name, "+", start_header="input|c1")
        query_proteins[protein.name] = entry_protein

        #TODO figure out why the +100. I assume for drawing
        total_lenght += entry_protein.nt_lenght + 100
    logging.debug("Finished parsing architecture input file")
    return query_proteins

def embl_to_genbank(emblfile):
    """
    Convert an embl file in such a way that the GenbankFile object can read it
    and convert it appropriately

    :param emblfile: a path to an embl file.
    :return: a string that can be read by a GenbankFile object
    """
    logging.debug("Converting to make embl file {} readable for the GenbankFile object".format(emblfile))
    try:
        with open(emblfile, "r") as f:
            file_text = f.read()
    except Exception as e:
        logging.critical("Invalid embl file {}. Exiting...".format(emblfile))
        raise MultiGeneBlastException("Invalid embl file {}.".format(emblfile))

    # make sure to remove potential old occurances of \r. Acts like a \n
    file_text = file_text.replace("\r", "\n")

    # do a basic check to see if the embl file is valid
    if "FT   CDS " not in file_text or ("\nSQ" not in file_text):
        logging.critical("Embl file {} is not properly formatted or contains no sequences. Exiting...".format(emblfile))
        raise MultiGeneBlastException("Embl file {} is not properly formatted or contains no sequences".format(emblfile))
    text_lines = file_text.split("\n")

    #change certain line starts
    line_count = 0
    while line_count < len(text_lines):
        line = text_lines[line_count]
        if line.startswith("FT"):
            text_lines[line_count] = line.replace("FT", "  ", 1)
        if line.startswith("SQ"):
            text_lines[line_count] = "ORIGIN"
        elif line.startswith("AC"):
            text_lines[line_count] = line.replace("AC   ", "ACCESSION   ", 1)
        elif line.startswith("DE"):
            #change all the definition lines to make them readable
            text_lines[line_count] = line.replace("DE   ", "DEFINITION  ", 1)
            line_count += 1
            def_line = text_lines[line_count]
            while not def_line.startswith("XX"):
                text_lines[line_count] = def_line.replace("DE   ", "            ", 1)
                line_count += 1
                def_line = text_lines[line_count]
        line_count += 1
    logging.debug("Embl file {} made readable for genbank file object.".format(emblfile))
    return "\n".join(text_lines)


##Functions necessary for this script
def get_sequence(fasta):
    """get the description and trimmed dna sequence"""
    in_file = open(fasta, 'r')
    content = in_file.readlines()
    in_file.close()
    content2 = []
    for i in content:
        if i != "":
            content2.append(i)
    content = content2
    while content[0] == "" or content[0] == "\n":
        content = content[1:]
    header = content[0]
    content = content[1:]
    content = [x.rstrip() for x in content]
    seq = "".join(content)
    if ">" not in header or ">" in seq:
        print("FASTA file not properly formatted; should be single sequence starting with '>' and sequence name.", file=sys.stderr)
        sys.exit(1)
    return seq

def parse_dna_from_embl(embl_string):
    "Parse DNA sequence from EMBL input"
    seq_array = []
    lines = embl_string.split('\n')
    for line in lines:
        if line.lower().find('sequence') > -1:
            continue
        line = line.strip()
        line = line.rstrip('0123456789')
        line = line.rstrip('/')
        line = line.strip()
        seq_array.append(line)

    return "".join(seq_array)

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

def internal_blast(user_options, query_proteins):
    """
    Run an internal blast of the query against itself.

    :param user_options: an Option object containing user specified options
    :param query_proteins: a dictionary of Protein objects saved by gene_name
    :return: a dictionary of blast results grouped around results that have
    returned blast hits against one another that are at least of a certain
    coverage and identity as defined by the user.
    """
    logging.info("Finding internal homologs...")
    internalhomologygroupsdict = {}
    clusternumber = 1

    #Make Blast db for using the sequences saved in query.fasta in the previous step
    #TODO when running this again after a crash with a different database, a new database
    #TODO continue: not created. This is problematic because it can lead to inconsistent nameing
    make_blast_db_command = "{}\\exec_new\\makeblastdb.exe -in query.fasta -out query_db -dbtype prot".format(CURRENTDIR)
    logging.debug("Started making internal blast database...")
    try:
        run_commandline_command(make_blast_db_command, max_retries=5)
    except MultiGeneBlastException:
        raise MultiGeneBlastException("Data base can not be created. The input query is most likely invalid.")
    logging.debug("Finished making internal blast database.")

    #Run and parse BLAST search
    blast_search_command = "{}\\exec_new\\blastp.exe  -db query_db -query query.fasta -outfmt 6" \
                           " -max_target_seqs 1000 -evalue 1e-05 -out internal_input.out" \
                           " -num_threads {}".format(CURRENTDIR, user_options.cores)

    logging.debug("Running internal blastp...")
    try:
        run_commandline_command(blast_search_command, max_retries=5)
    except MultiGeneBlastException:
        raise MultiGeneBlastException("Running blast when finding internal homologs returns multiple errors.")
    logging.debug("Finished running internal blastp.")

    try:
        with open("internal_input.out","r") as f:
            blastoutput = f.read()
    except Exception:
        logging.critical("Something went wrong reading the blast output file. Exiting...")
        raise MultiGeneBlastException("Something went wrong reading the blast output file.")

    #extract blast hits into a dictionary
    blast_dict = blast_parse(user_options, query_proteins, blastoutput)

    #because query proteins are sorted this is fine
    query_cluster = Cluster(list(query_proteins.values())[0].start, list(query_proteins.values())[-1].stop, "query")
    groups = []
    for results in blast_dict.values():
        for blast_result in results:
            query_cluster.add_protein(query_proteins[blast_result.subject], blast_result)
    return query_cluster

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

#TODO see where these are used
def sortdictkeysbyvalues(dict):
    items = [(value, key) for key, value in list(dict.items())]
    items.sort()
    return [key for value, key in items]

def sortdictvaluesbykeys(dict):
    items = [(key, value) for key, value in list(dict.items())]
    items.sort()
    return [value for key, value in items]

def sortdictkeysbyvaluesrev(dict):
    items = [(value, key) for key, value in list(dict.items())]
    items.sort()
    items.reverse()
    return [key for value, key in items]

def sortdictkeysbyvaluesrevv(dict):
    items = [(value, key) for key, value in list(dict.items())]
    items.sort()
    items.reverse()
    return [value for value, key in items]

class BlastResult:
    """
    Simply makes it easier to see what values are used for comparissons
    See also http://www.metagenomics.wiki/tools/blast/blastn-output-format-6
    """
    def __init__(self, tabs, query_protein):
        """
        :param tabs: is a list values in a line of blast output split on \t
        :param query_protein: a Protein object that corresponds to the query
        protein that was put in.
        """
        self.line = tabs

        #protein object that is the query
        self.query_protein = query_protein

        self.query = tabs[0]
        self.subject = tabs[1]
        self.percent_identity = float(tabs[2])
        self.allignment_lenght = int(tabs[3])
        #number of mismatches
        self.mismatch = int(tabs[4])
        #number of gap openings
        self.gapopen = int(tabs[5])

        #allignment in query
        self.query_start = int(tabs[6])
        self.query_stop = int(tabs[7])

        #allignment in subject
        self.subject_start = int(tabs[8])
        self.subject_stop = int(tabs[9])

        #expected value
        self.evalue = float(tabs[10])
        self.bit_score = float(tabs[11])

        # calculated values
        #percentage of the query that is covered by the subject
        self.percent_coverage = self.allignment_lenght / self.query_protein.aa_lenght * 100

    def summary(self):
        return "\t".join(self.line)


def blast_parse(user_options, query_proteins, blast_output):
    """
    Parse the output of a blast search using BlastResult objects.

    :param user_options: a Option object containing the options supplied by the user
    :param query_proteins: a dictionary that links gene_names against Protein objects
    :param blast_output: the output produced by NCBI BLAST+ blastp algorithm version
    2.2.18
    :return: a dictionary that contains a key for every query that had a siginificant
    blast result with the values being BlastResult objects that have the key as
    query and are above user defined treshholds.
    """
    logging.debug("Started parsing blast output...")
    #remove the last split. It is an empty line
    blast_lines = blast_output.split("\n")[:-1]

    #Filter for unique blast comparissons
    query_subject_combinations = []
    unique_blast_results = {}
    for line in blast_lines:
        tabs = line.split("\t")
        query = tabs[0]
        subject = tabs[1]
        query_subject_combination = query + subject
        if not (query_subject_combination in query_subject_combinations):
            query_subject_combinations.append(query_subject_combination)
            blast_result = BlastResult(tabs, query_proteins[query])
            unique_blast_results[blast_result.query + blast_result.subject] = blast_result

    #Filters blastlines to get rid of hits that do not meet criteria and save
    #hits that do in a dictionary that links a protein name to BlastHit objects.
    internal_homolog_dict = {}
    count = 0
    #TODO remove the or statement, now here for consistency sake
    for result in unique_blast_results.values():
        if int(result.percent_identity) > user_options.minpercid and\
                (int(result.percent_coverage) > user_options.minseqcov or result.allignment_lenght > 40):
            count += 1
            if result.query not in internal_homolog_dict:
                internal_homolog_dict[result.query] = [unique_blast_results[result.query + result.subject]]
            else:
                internal_homolog_dict[result.query].append(unique_blast_results[result.query + result.subject])
    logging.debug("Finished parsing blast output.")
    return internal_homolog_dict

def db_blast(query_proteins, user_options):
    """
    Run blast agianst the user defined database

    :param query_proteins: a dictionary of Protein objects
    :param user_options: an Option object with user options
    :return: the blast output
    """
    logging.info("Running NCBI BLAST+ searches on the provided database {}..".format(user_options.db))
    if user_options.dbtype == "prot":
        command_start = "{}\\exec_new\\blastp.exe".format(CURRENTDIR)
    else:
        command_start = "{}\\exec_new\\tblastn.exe".format(CURRENTDIR)

    complete_command = "{} -db {} -query {}\\query.fasta -outfmt 6 -max_target_seqs" \
                       " {} -evalue 1e-05 -out {}\\input.out -num_threads {}"\
        .format(command_start, user_options.db, os.getcwd(), user_options.hitspergene,
                os.getcwd(), user_options.cores)

    logging.debug("Started blasting against the provided database...")
    db_blast_process = Process(target=run_commandline_command, args=[complete_command])
    db_blast_process.start()
    while db_blast_process.is_alive():
        continue
    #when the process exits with an error
    if db_blast_process.exitcode != 0:
        raise MultiGeneBlastException("Blasting against the provided database returns multiple errors.")
    logging.debug("Finished blasting")

    try:
        with open("input.out","r") as f:
            blastoutput = f.read()
    except:
        logging.critical("Cannot open the file that contains the database blast output. Exiting...")
        raise MultiGeneBlastException("Cannot open the file that contains the database blast output.")
    return blastoutput

def parse_db_blast(user_options, query_proteins, blast_output):
    """
    Wrapper function for blast_parse where an exit can be catched if no blast
    results are returned

    :param user_options: a Option object containing the options supplied by the user
    :param query_proteins: a dictionary that links gene_names against Protein objects
    :param blast_output: the output produced by NCBI BLAST+ blastp algorithm version
    2.2.18
    :return: a dictionary that contains a key for every query that had a siginificant
    blast result with the values being BlastResult objects that have the key as
    query and are above user defined treshholds.
    """
    blast_dict = blast_parse(user_options, query_proteins, blast_output)
    if len(blast_dict) == 0:
        logging.info("No blast hits encountered against the provided database."
                     " Maybe try and lower min_percentage_identity ({}) or "
                     "minimum_sequence_coverage ({})"
                     " Extiting...".format(user_options.minpercid,
                                           user_options.minseqcov))
        sys.exit(0)
    return blast_dict

def generate_rgbscheme(nr):
    #TODO look furter into function if there is time look into antismash for methods here
    usablenumbers = [1,2,4,8,12,18,24,32,48,64,10000]
    lengthsdict = {1:[1,1,1],2:[1,1,2],4:[1,2,2],8:[2,2,2],12:[2,2,3],18:[2,3,3],24:[3,3,3],32:[3,3,4],48:[3,4,4],64:[4,4,4]}
    shortestdistance = 10000
    for i in usablenumbers:
        distance = i - nr
        if distance >= 0:
            if distance < shortestdistance:
                shortestdistance = distance
                closestnr = i
    toohigh = "n"
    if closestnr == 10000:
        toohigh = "y"
        closestnr = 64
    xyznumbers = lengthsdict[closestnr]
    x = xyznumbers[0]
    y = xyznumbers[1]
    z = xyznumbers[2]
    xpoints = []
    xpoint = (255/z)/2
    for i in range(x):
        xpoints.append(xpoint)
        xpoint += (255/x)
    ypoints = []
    ypoint = (255/z)/2
    for i in range(y):
        ypoints.append(ypoint)
        ypoint += (255/y)
    zpoints = []
    zpoint = (255/z)/2
    for i in range(z):
        zpoints.append(zpoint)
        zpoint += (255/z)
    colorlist = []
    for i in xpoints:
        for j in ypoints:
            for k in zpoints:
                rgb = "rgb(" + str(i) + "," + str(j) + "," + str(k) + ")"
                colorlist.append(rgb)
    if toohigh == "y":
        colorlist = colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist
    if closestnr == 24:
        colorlist = colorlist[:15] + colorlist[18:]
    if closestnr == 32:
        colorlist = colorlist[:21] + colorlist[24:]
    colorlist2 = []
    if closestnr == 1:
        colorlist2.append("red")
    if closestnr == 2:
        colorlist2.append("red")
        colorlist2.append("green")
    if closestnr == 4:
        colorlist2.append("red")
        colorlist2.append("green")
        colorlist2.append("blue")
        colorlist2.append("yellow")
    if closestnr == 8:
        neworder=[4,1,2,5,6,7,3,0]
        colorlist2 = [colorlist[i] for i in neworder]
    if closestnr == 12:
        neworder=[6,3,5,9,7,2,11,4,8,1,10,0]
        colorlist2 = [colorlist[i] for i in neworder]
    if closestnr == 18:
        neworder=[9,6,2,14,15,8,12,10,3,5,7,11,4,1,16,13,0]
        colorlist2 = [colorlist[i] for i in neworder]
    if closestnr == 24:
        neworder=[15,12,9,6,5,0,21,1,16,14,8,17,2,23,22,3,13,7,10,4,18,20,19,11]
        colorlist2 = [colorlist[i] for i in neworder]
    if closestnr == 32:
        neworder = [21,19,27,6,8,1,14,7,20,13,9,30,4,23,18,12,5,29,24,17,11,31,2,28,22,15,26,3,20,16,10,25]
        colorlist2 = [colorlist[i] for i in neworder]
    if closestnr > 32:
        random.shuffle(colorlist)
        colorlist2 = colorlist
    colorlist = colorlist2
    return colorlist

def _gene_arrow(start,end,strand,color,base,height):
    #TODO look to clean up this function if there is time left
    halfheight = height/2
    if start > end:
        start2 = end
        end2 = start
        start = start2
        end = end2
    oh = ShapeBuilder()
    if (end - start) < halfheight:
        if (strand == "+"):
            pointsAsTuples=[(start,base),
                            (end,base - halfheight),
                            (start,base - height),
                            (start,base)
                            ]
        if (strand == "-"):
            pointsAsTuples=[(start,base - halfheight),
                            (end,base - height),
                            (end,base),
                            (start,base - halfheight)
                            ]
    else:
        if (strand == "+"):
            arrowstart = end-halfheight
            pointsAsTuples=[(start,base),
                            (arrowstart,base),
                            (end,base-halfheight),
                            (arrowstart,base - height),
                            (start,base - height),
                            (start,base)
                            ]
        if (strand == "-"):
            arrowstart = start + halfheight
            pointsAsTuples=[(start,base - halfheight),
                            (arrowstart,base - height),
                            (end,base - height),
                            (end,base),
                            (arrowstart,base),
                            (start,base - halfheight)
                            ]
    pg=oh.createPolygon(points=oh.convertTupleArrayToPoints(pointsAsTuples),strokewidth=1, stroke='black', fill=color)
    return pg

def relativepositions(starts, ends, largestclustersize, screenwidth):
    rel_starts = []
    rel_ends = []
    #Assign relative start and end sites for visualization
    lowest_start = int(starts[0])
    leftboundary = lowest_start
    for i in starts:
        i = float(float(int(i) - int(leftboundary)) / largestclustersize) * float(screenwidth * 0.75)
        i = int(i)
        rel_starts.append(i)
    for i in ends:
        i = float(float(int(i) - int(leftboundary)) / largestclustersize) * float(screenwidth * 0.75)
        i = int(i)
        rel_ends.append(i)
    return [rel_starts,rel_ends]

def startendsitescheck(starts,ends):
    #Check whether start sites are always lower than end sites, reverse if necessary
    starts2 = []
    ends2 = []
    a = 0
    for i in starts:
        if int(i) > int(ends[a]):
            starts2.append(ends[a])
            ends2.append(i)
        else:
            starts2.append(i)
            ends2.append(ends[a])
        a += 1
    ends = ends2
    starts = starts2
    return [starts,ends]

def calculate_colorgroups(blast_dict, query_cluster):
    #TODO it is the case that genes that do not have any blast hits are ignored in the query

    #first get groups internal homologs in the query cluster
    groups = []
    for key in query_cluster.blast_hit_proteins:
        # create a set to automatically filter for duplicates
        group = {key}
        for blast_result in query_cluster.blast_hit_proteins[key]:
            group.add(blast_result.query)

        # make sure that a no member of the new group is present in any of the
        # current groups. If that is the case add to current groups, else make
        # new group
        added = False
        for present_group in groups:
            for value in present_group:
                if value in group:
                    present_group.update(group)
                    added = True
                    break
        if not added:
            groups.append(group)

    color_groups = []
    for st in groups:
        color_group = set()
        for protein_name in st:
            if not protein_name in blast_dict:
                continue
            color_group.add(protein_name)
            for blast_result in blast_dict[protein_name]:
                color_group.add(blast_result.subject)
        if len(color_group) > 0:
            color_groups.append(list(color_group))
    # Generate RGB scheme
    rgb_color_scheme = generate_rgbscheme(len(color_groups))
    gene_color_dict = {}
    for gr_indx, group in enumerate(color_groups):
        for protein_name in group:
            gene_color_dict[protein_name] = rgb_color_scheme[gr_indx]
    return gene_color_dict

def equal_to_query_strand_orientation(cluster):
    equal_score = 0
    for protein_name, result_set in cluster.blast_hit_proteins.items():
        for blast_result in result_set:
            # compare the query protein strand with that of the cluster
            if blast_result.query_protein.strand == cluster.get_protein(protein_name).strand:
                equal_score += 1
            else:
                equal_score -= 1
    return equal_score >= 0

def draw_gene_cluster(proteins, dpn, index, screenwidth, color_dict, line = True, reverse=False):
    groups = []
    builder = ShapeBuilder()
    group = g()
    if line:
        group.addElement(builder.createLine(10, 35 + HITS_PER_PAGE * index, screenwidth, 35 + HITS_PER_PAGE * index, strokewidth=1, stroke="grey"))
    groups.append(group)
    rel_cluster_start = int(proteins[0].start * dpn)
    rel_cluster_stop = int(proteins[-1].stop * dpn)

    center_distance = (screenwidth - (rel_cluster_stop - rel_cluster_start)) / 2
    if reverse:
        center_distance = - (screenwidth - (rel_cluster_start - rel_cluster_stop)) / 2

    for prot in proteins:
        group = g()
        if prot.name in color_dict:
            color = color_dict[prot.name]
        else:
            color = "#FFFFFF"
        rel_start = abs(center_distance + (prot.start * dpn - rel_cluster_start))
        rel_stop = abs(center_distance + (prot.stop * dpn - rel_cluster_start))
        strand = prot.strand
        if reverse:
            if strand == "+":
                strand = "-"
            else:
                strand = "+"
        group.addElement(_gene_arrow(rel_stop, rel_start, strand, color, 40 + HITS_PER_PAGE * index, 10))

        group.set_id("all_" + str(index) + "_" + "%s" % prot.name)
        groups.append(group)
    return groups

def clusterblastresults(query_cluster, clusters, gene_color_dict, user_options):

    screenwidth = user_options.screenwidth * 0.8

    svg_image = svg(x = 0, y = 0, width = screenwidth, height = 2770)
    viewbox = "0 0 " + str(screenwidth) + " " + str(2680)
    svg_image.set_viewBox(viewbox)
    svg_image.set_preserveAspectRatio("none")

    #first element is the reference for the rest of the clusters
    all_clusters = clusters + [query_cluster]

    #biggest core size
    biggest_size = 0
    for cluster in all_clusters:
        start, stop = cluster.core_start_stop()
        size = (stop - start) + 2 * SVG_CORE_EXTENSION
        if size > biggest_size:
            biggest_size = size

    #a factor that can be applied to get a scaled distance relative to the screen
    #size
    distance_per_nucleotide = screenwidth / biggest_size

    #first draw the query
    line = True
    if user_options.architecture_mode:
        line = False
    groups = draw_gene_cluster(list(query_cluster.proteins.values()), distance_per_nucleotide, 0, screenwidth, gene_color_dict, line=line)
    for group in groups:
        svg_image.addElement(group)

    query_cluster_size = query_cluster.core_start_stop()[1] - query_cluster.core_start_stop()[0]

    for index, cluster in enumerate(clusters):
        #make sur the cluster is sorted. This should be the case but calling
        #this method does no harm when the cluster is sorted already
        cluster.sort()

        #take an a good size sample of proteins
        reduced_proteins = reduced_cluster_proteins(cluster, query_cluster_size)

        # check if the orientation of a cluster is inverted compared to the query
        # if so reverse the proteins.
        reverse = False
        if not equal_to_query_strand_orientation(cluster):
            reduced_proteins.reverse()
            reverse = True

        groups = draw_gene_cluster(reduced_proteins, distance_per_nucleotide, index + 1, screenwidth, gene_color_dict, reverse=reverse)
        for group in groups:
            svg_image.addElement(group)
    return svg_image

def reduced_cluster_proteins(cluster, query_cluster_size):
    # reduce the amount of proteins in the cluster by ignoreing proteins that
    # are to far from blast hits
    cluster_proteins = list(cluster.proteins.values())
    core_start, core_stop = cluster.core_start_stop()
    core_size = core_stop - core_start
    if query_cluster_size * 0.8 > core_size:
        extra_distance = (query_cluster_size * 0.8 - core_size) / 2
        core_start = max(core_start - extra_distance, cluster.start_stop()[0])
        core_stop = min(core_stop + extra_distance, cluster.start_stop()[1])
    min_start = core_start - SVG_CORE_EXTENSION
    max_end = core_stop + SVG_CORE_EXTENSION

    reduced_proteins = []
    for prot in cluster_proteins:
        if prot.start > min_start and prot.stop < max_end:
            reduced_proteins.append(prot)
    return reduced_proteins

def parse_absolute_paths(infile):
    #Parse absolute paths if found
    originalfilename = infile
    if "/" in infile or "\\" in infile:
        lastpos = max([infile.rfind("\\"),infile.rfind("/")])
        originpath = infile[:(lastpos + 1)]
        infile = infile[(lastpos + 1):]
        #if os.getcwd() != originalfilename[:lastpos] and os.getcwd() != originalfilename[:lastpos].replace("/","\\"):
        #  shutil.copyfile(originalfilename, infile)
    return infile


def load_genecluster_info(dbname, allgenomes):
    #Load gene cluster info to memory
    DBPATH = os.environ['BLASTDB']
    clusters = {}
    allgenomes_tags = [genomename[:6] for genomename in allgenomes]
    for i in fileinput.input(DBPATH + os.sep + dbname + "_all_descrs.txt"):
        tabs = i.split("\t")
        if len(tabs) > 0 and (tabs[0] in allgenomes or tabs[0] in allgenomes_tags):
            accession = tabs[0]
            clusterdescription = tabs[1]
            clusters[accession] = clusterdescription
    nucdescriptions = clusters
    frame_update()
    return nucdescriptions, clusters


def load_databases(query_proteins, blast_dict, user_options):
    #TODO improve this to only include scaffolds with a hit. Needs improvement in general in combination with improving the database
    logging.info("Loading GenBank positional info into memory...")
    db_path = os.environ["BLASTDB"]

    picklefile = open("{}\\{}.pickle".format(db_path, user_options.db), "rb")
    p_gbk = pickle.load(picklefile)
    picklefile.close()

    return p_gbk
    # #Load GenBank positional info into memory
    # if user_options.dbtype == "prot":
    #     hit_proteins = load_dbproteins_info(query_proteins, blast_dict, user_options.db)
    #     proteininfo = load_other_genes(allgenomes, proteininfo, dbname, blastdict)
    # else:
    #     allgenomes, nucdict, proteininfo = load_ndb_info(querylist, blastdict, dbname)
    # nucdescriptions, clusters = load_genecluster_info(dbname, allgenomes)
    # return nucdescriptions, nucdict, proteininfo

def find_gene_clusters(blast_dict, user_options, database):
    """
    Find clusters of genes in the blast_dict

    :param blast_dict: a dictionary that contains all query genes with a list
    of BlastResults objects that where against these genes.
    :param user_options: Option object of user specified options
    :param database: a database loaded from a pickled class.
    :return: A list of Cluster objects that are clusters found in the blast_dict
    plus all proteins within those cluster regions
    """
    logging.info("Started finding clusters in the blast results...")
    hits_per_contig, query_per_blast_hit = sort_hits_per_scaffold(blast_dict, database)
    clusters_per_contig = find_blast_clusters(hits_per_contig, query_per_blast_hit, user_options.distancekb)
    add_additional_proteins(clusters_per_contig, database)
    clusters = []

    #add clusters to one big list and sort them
    for contig in clusters_per_contig:
        for cluster in clusters_per_contig[contig]:
            cluster.sort()
            clusters.append(cluster)
    return clusters

def sort_hits_per_scaffold(blast_dict, database):
    """
    Extract all the different scaffolds from the blast_dict and sort the
    blast_hits on each contig

    :param blast_dict: a dictionary that contains all query genes with a list
    of BlastResults objects that where against these genes.
    :param database: a database loaded from a pickled class.
    :return: A dictionary of blast hits with scaffolds as keys sorted on
    location and a dictionary that linkes the subject of a BlastResult to the
    result itself.
    """
    logging.debug("Sorting blast-hits per scaffold...")
    hits_per_contig = {}
    query_per_blast_hit = {}
    for hits in blast_dict.values():
        for hit in hits:
            #extract the subject proteins from the Blastresults using the database
            subject_protein = database.proteins[hit.subject]
            if hit.subject in query_per_blast_hit:
                query_per_blast_hit[hit.subject].append(hit)
            else:
                query_per_blast_hit[hit.subject] = [hit]
            if subject_protein.genbank_file not in hits_per_contig:
                hits_per_contig[subject_protein.genbank_file] = set([subject_protein])
            else:
                hits_per_contig[subject_protein.genbank_file].add(subject_protein)
    #sort each contig using start and stop locations
    for scaffold in hits_per_contig:
        hits_per_contig[scaffold] = list(hits_per_contig[scaffold])
        hits_per_contig[scaffold].sort(key=lambda x: (x.start, x.stop))
    return hits_per_contig, query_per_blast_hit

def find_blast_clusters(hits_per_contig, query_per_blast_hit, extra_distance):
    """
    Find clusters for each scaffold.

    :param hits_per_contig: A dictionary of blast hits with scaffolds as keys
    sorted on location.
    :param query_per_blast_hit: a dictionary that linkes the subject of a
    BlastResult to a list of results with that subject
    :param extra_distance: the distance that is allowed between 2 blast hits
    for them to be considered part of the same cluster
    :return: a dictionary that contains contigs as keys and lists of Cluster
    objects as values.
    """
    #find all the clusters
    logging.debug("Finding clusters using blasthits...")

    clusters_per_contig = {}
    for scaffold in hits_per_contig:

        scaffold_proteins = hits_per_contig[scaffold]
        start = scaffold_proteins[0].start
        blast_hits = [scaffold_proteins[0]]
        index = 0
        clusters = []

        while index + 1 < len(scaffold_proteins):
            protein = scaffold_proteins[index]
            next_protein = scaffold_proteins[index + 1]
            blast_hits.append(protein)
            if protein.stop + extra_distance < next_protein.start:
                c = Cluster(start - extra_distance, protein.stop + extra_distance, protein.genbank_file)
                #add all the blast_hits to the cluster together with the query
                #gene they originate from
                for hit in blast_hits:
                    for result in query_per_blast_hit[hit.name]:
                        c.add_protein(hit, result)
                clusters.append(c)
                start = next_protein.start
                blast_hits = [next_protein]
            index += 1
        #make sure to add the final cluster
        c = Cluster(start - extra_distance, scaffold_proteins[index].stop + extra_distance, scaffold_proteins[index].genbank_file)
        for hit in blast_hits:
            for result in query_per_blast_hit[hit.name]:
                c.add_protein(hit, result)
        clusters.append(c)
        tot = 0

        #all cluster objects are from the same contig
        clusters_per_contig[clusters[0].contig] = clusters
    return clusters_per_contig

def add_additional_proteins(clusters_per_contig, database):
    """
    Add proteins from the database that are within the bounds of the cluster
    but no blzst hit, to the clusters.

    :param clusters_per_contig: a dictionary that contains contigs as keys and
    lists of Cluster objects as values.
    :param database: a database loaded from a pickled class.
    """
    logging.debug("Adding additional proteins...")
    for contig in clusters_per_contig:
        contig_proteins = database.contigs[contig].proteins.values()
        # make a copy of all the proteins of a contig and sort them
        cluster_index = 0
        clusters = clusters_per_contig[contig]

        #because all clusters are sorted aswell as the proteins one loop suffices
        #to assign all proteins
        for protein in contig_proteins:
            cluster = clusters[cluster_index]

            if (protein.start >= cluster.start and protein.start <= cluster.stop):
                cluster.add_protein(protein)
                #if a protein is added make sure to add the same protein to a
                #different cluster that can be overlapping in the extra region part
                if cluster_index + 1 < len(clusters):
                    next_cluster = clusters[cluster_index + 1]
                    if protein.stop >= next_cluster.start and protein.stop <= next_cluster.stop:
                        next_cluster.add_protein(protein)

            elif protein.start > cluster.start:
                cluster_index += 1
                #when all clusters have been covered
                if cluster_index >= len(clusters):
                    break
                cluster = clusters[cluster_index]
                if (protein.start >= cluster.start and protein.start <= cluster.stop):
                    cluster.add_protein(protein)


class Cluster:
    """
    Tracks all proteins in a cluster. The proteins should start within the
    cluster bounds. The cluster can be sorted but only when specifically asked
    by the user. The cluster is then guaranteed sorted unitl new proteins are
    added.
    """
    COUNTER = 0
    def __init__(self, start, stop, contig):
        """
        :param start: an integer that is the start of the cluster. Is allowed to
        be negative, this just means that the cluster starts from 0
        :param stop: an Integer that is the end of the cluster. Is allowed to be
        bigger then the contig. Simply means that the cluster stops at the end
        of the contig
        :param contig: the name of the contig this cluster is located on.
        """
        #track an increasing cluster number
        self.no = self.COUNTER
        Cluster.COUNTER += 1

        #start and stop are allowed to be negative or bigger then possible
        self.start = start
        self.stop = stop
        self.contig = contig
        #all proteins in the cluster
        self.__proteins = OrderedDict()
        self.__protein_indexes = OrderedDict()
        self.sorted = True

        #all protein names that are blast hits linked to the query protein that
        #that the blast hit was against
        self.__blast_proteins = OrderedDict()

        #a dictionary for tracking individual parts of the score
        self.__score = {"synteny" : 0, "accumulated_blast_score" : 0, "number unique hits" : 0}

    def add_protein(self, protein, query_blast_hit = False):
        """
        Controlled way of adding proteins to the cluster.

        :param protein: a Protein object
        :param query_blast_hit: a potential Protein object that is part of the
        query that defines to what query gene the 'protein' has a blast hit.
        """
        if query_blast_hit:
            if protein.name in self.__blast_proteins:
                self.__blast_proteins[protein.name].add(query_blast_hit)
            else:
                #create a set to auto filter duplicates
                self.__blast_proteins[protein.name] = {query_blast_hit}
        if not protein.name in self.__proteins:
            self.__proteins[protein.name] = protein
            self.__protein_indexes[protein.name] = len(self.__proteins) - 1
            self.__sorted = False

    def get_protein(self, name):
        return self.__proteins[name]

    @property
    def score(self):
        """
        Sum of the score dictionary

        :return: a float
        """
        return sum(self.__score.values())

    def start_stop(self):
        """
        Retrieve the start and stop coordinate of the cluster. This has tot
        sort the cluster if it is unsorted. Take care calling this method to
        often when adding proteins at the same time

        :return: two integers (start, stop)
        """
        self.sort()
        proteins = list(self.__proteins.values())
        return proteins[0].start, proteins[-1].stop

    def core_start_stop(self):
        """
        Retrieve the start and stop coordinate of the core of the cluster.
        The core is defines as the proteins that had blast hits against the
        query. This has to sort the cluster if it is unsorted. Take care
        calling this method to often when adding proteins at the same time

        :return: two integers (start, stop)
        """
        self.sort()
        proteins = list(self.__blast_proteins.keys())
        return self.__proteins[proteins[0]].start, self.__proteins[proteins[-1]].stop

    def set_score(self, synteny, accumulated_blast_score, number_unique_hits):
        """
        Set the score in a dictionary as 3 componenets

        :param synteny: the synteny score as float
        :param accumulated_blast_score: the accumulated blast bit score
        / 1.000.000
        :param number_unique_hits: integer that is the number of unique blast hits
        """
        self.__score["synteny"] = synteny
        self.__score["accumulated_blast_score"] = accumulated_blast_score
        self.__score["number unique hits"] = number_unique_hits

    def segmented_score(self):
        return self.__score

    @property
    def proteins(self):
        return self.__proteins

    @property
    def blast_hit_proteins(self):
        return self.__blast_proteins

    def sort(self):
        """
        Sort the __proteins dict and __protein_indexes dict if the cluster is
        not sorted.

        A cluster becomes unsorted if proteins are added. Calling sort on a
        sorted cluster does not do anything.
        """
        if not self.__sorted:
            self.__proteins = OrderedDict(sorted(self.__proteins.items(), key=lambda item: (item[1].start, item[1].stop)))
            self.__protein_indexes = OrderedDict((name, index) for index, name in enumerate(self.__proteins))
            self.__blast_proteins = OrderedDict(sorted(self.__blast_proteins.items(), key=lambda item: (self.__proteins[item[0]].start, self.__proteins[item[0]].stop)))
            self.__sorted = True

    def index(self, protein_name):
        """
        The index of a protein within the self.__proteins OrderedDictionary

        :param protein_name: a string that is the name of a protein
        :return: the index of the protein name in the self.__protein indexes
        dictionary
        """
        return self.__protein_indexes[protein_name]

    def short_summary(self):
        """
        A short summary of the cluster

        :return: a string
        """
        return "Cluster {}\t{}\t{:.7f}\t{}".format(self.no, self.contig, self.score, len(self.__blast_proteins))

    def summary(self):
        """
        A comprehensive string summary of a cluster, containing general
        information, and tables for genes in the cluster and the blast hits

        :return: a string that contains the information described above
        """
        full_summary = ""
        #general information
        full_summary += "- Source contig: {}\n".format(self.contig)
        full_summary += "- Number of proteins with BLAST hits to this cluster: {}\n".format(len(self.__blast_proteins))
        full_summary += "- MultiGeneBlast 2.0 score: {:.7f}\n  - synteny score: {}\n  " \
                        "- accumulated blast bit score / 1.000.000: {:.7f}\n  " \
                        "- Number of unique blast hits: {}\n\n".\
            format(self.score, self.__score["synteny"], self.__score["accumulated_blast_score"], self.__score["number unique hits"])

        #Table of all proteins in cluster
        full_summary += ">Table of genes in cluster:\n"
        full_summary += "Name\tstart\tstop\tstrand\tannotation\tlocus_tag\n"

        for protein in self.__proteins.values():
            full_summary += "{}\n".format(protein.summary())
        full_summary += "<\n"

        #Table of all the blast hits
        full_summary += "\n>Table of all blast hits:\n"
        full_summary += "query prot\tsubject prot\tperc ident\tallignment len" \
                        "\tmismatch\tgap open\tquery start\tquery stop\tsubject" \
                        " start\tsubject stop\te-value\tblast bit score\n"

        for result_set in self.__blast_proteins.values():
            for blast_result in result_set:
                full_summary += "{}\n".format(blast_result.summary())
        full_summary += "<\n"
        return full_summary


def score_clusters(clusters, query_proteins, user_options):
    """
    Score all clusters based on cumulative blast scores, the total number of
    unique queries that got a blast hit and synteny score when running homology
    search

    :param clusters: A list of Cluster objects
    :param query_proteins: A list of Protein objects that are present in the
    query
    :param user_options: An Option object that contains user specified options
    """
    logging.info("Scoring clusters...")
    index_query_proteins = OrderedDict((prot, i) for i, prot in enumerate(query_proteins))

    #calculate the cumulative blast scores over all blast hits of a cluster
    #TODO figure out why these scores differ slightly from mgb. I have looked but not found :(
    for cluster in clusters:
        cum_blast_score = 0
        hitnr = 0
        hit_positions_dict = {}
        for result_set in cluster.blast_hit_proteins.values():
            for blast_result in result_set:
                cum_blast_score += blast_result.bit_score / 1_000_000.0
                #get the best scoring blast result for each query entry
                if blast_result.query not in hit_positions_dict:
                    hitnr += 1
                    hit_positions_dict[blast_result.query] = blast_result
                else:
                    prev_br = hit_positions_dict[blast_result.query]
                    if prev_br.evalue > blast_result.evalue:
                        hit_positions_dict[blast_result.query] = blast_result

        synteny_score = 0
        if not user_options.architecture_mode:
            #make a tuple of index in the query and index in the cluster
            hit_positions = []
            for query_protein in hit_positions_dict:
                blast_result = hit_positions_dict[query_protein]
                hit_positions.append((index_query_proteins[blast_result.query],
                                     cluster.index(blast_result.subject)))
            synteny_score = score_synteny(hit_positions)

        #assign the score
        cluster.set_score(synteny_score * user_options.syntenyweight, cum_blast_score, hitnr)

    #sort the clusters on score
    clusters.sort(key=lambda x: x.score, reverse=True)

def score_synteny(hit_positions):
    """
    Calculate synteny score using a list of hit positions.

    :param hit_positions: A list of tuples containing the index of a protein in
    the query and the index of a protein in the cluster.
    :return: An integer that is a synteny score which is the number of adjacent
    genes in the same order as the query -1

    This function checks if the first indexes in consecutive tuples go up or
    down by 1 each time.
    """
    # modified from antismash
    score = 0
    for index, pos in enumerate(hit_positions[:-1]):
        query, subject = pos
        next_query, next_subject = hit_positions[index + 1]
        if (abs(query - next_query) < 2) and abs(query - next_query) == abs(subject - next_subject):
            score += 1
    return score

def write_txt_output(query_proteins, clusters, blast_output, user_options):
    """
    Write the output into a file containing information on all clusters.

    :param query_proteins: a list of Protein objects that where present in the
    provided query
    :param clusters: a list of Cluster objects
    :param blast_output: a dictionary that contains a key for every query that
    had a siginificant blast result with the values being BlastResult objects
    :param user_options: an Option object with user define doptions

    Note the formatting of the file is defined as follows:
    - '>>' signify a header of a new piece of information, the text directly after
    it is the descriptor
    - '>' signifies the start of a table followed by the descriptor of what is in
    the table. The first line of the table should always be the header.
    - '<' indicates the end of a table
    - '-' indicates a piece of information, taht can be gerarded as general
    facts of a certain header
    """
    logging.info("Writing .mgb output file into {}...".format(os.path.join(user_options.outdir, "clusterblast_output.mgb")))
    try:
        out_file = open("{}\clusterblast_output.mgb".format(user_options.outdir),"w")
    except (OSError, IOError):
        logging.critical("Writing .mgb file has failed. Can not opern file {}".format(os.path.join(user_options.outdir, "clusterblast_output.mgb")))
        raise MultiGeneBlastException("Can not open the output file {}.".format(os.path.join(user_options.outdir, "clusterblast_output.mgb")))

    #write for what database
    out_file.write("ClusterBlast scores for: {}\{}\n\n".format(os.environ["BLASTDB"], user_options.db))

    #write all the query proteins
    out_file.write(">>Query proteins:\n")
    out_file.write(">Table of query proteins:\n")
    out_file.write("Name\tstart\tstop\tstrand\tannotation\tlocus_tag\n")
    for protein in query_proteins.values():
        out_file.write("{}\n".format(protein.summary()))
    out_file.write("<\n")

    #write general information about the clusters
    out_file.write("\n>>Significant clusters:\n")
    out_file.write(">Table of the summary of significant clusters:\n")
    out_file.write("cluster no\tcontig id\tscore\tnumber of unique blast hits\n")
    for cluster in clusters:
        out_file.write("{}\n".format(cluster.short_summary()))
    out_file.write("<\n\n")

    #write more comprehensive information for each cluster
    for cluster in clusters:
        out_file.write(">> Details Cluster {}\n".format(cluster.no))
        out_file.write("{}\n".format(cluster.summary()))
    out_file.close()

def write_svg_files(clusters, user_options, query_cluster, blast_dict):
    logging.info("Writing visualization SVGs and XHTML...")
    try:
        os.mkdir("svg/")
    #if the folder exists ignore the error
    except(IOError,OSError):
        pass
    gene_color_dict = calculate_colorgroups(blast_dict, query_cluster)
    for index, cluster in enumerate(clusters):
        svg_clusters = [cluster]
        svg_image = clusterblastresults(query_cluster, svg_clusters, gene_color_dict, user_options)

        outfile = open("svg/clusterblast_{}.svg".format(cluster.no),"w")
        outfile.write(svg_image.getXML())
        outfile.close()
        if index % HITS_PER_PAGE == 0:
            page_nr = int(index / HITS_PER_PAGE)
            #when max pages are reached stop
            if page_nr > user_options.pages:
                break
            #Create svgs for multiple gene cluster alignment
            page_genes = clusters[page_nr * HITS_PER_PAGE : min((page_nr + 1) * HITS_PER_PAGE, len(clusters))]
            svg_image = clusterblastresults(query_cluster, page_genes , gene_color_dict, user_options)

            outfile = open("svg/clusterblast_page{}_all.svg".format(page_nr + 1),"w")
            outfile.write(svg_image.getXML())
            outfile.close()
            logging.info("Visual output page {}/{} created.".format(page_nr + 1, min(user_options.pages, int((len(clusters) - 1) / HITS_PER_PAGE) + 1)))

def runmuscle(args):
    os.system("muscle " + args)

def align_muscle(include_muscle, colorschemedict, seqdict):
    #Create Muscle alignments of colourgroups
    musclegroups = []
    if include_muscle == "y":
        log("Aligning homologous sequences with Muscle")
        try:
            os.mkdir("fasta")
        except(IOError,OSError):
            pass
        orthogroupsdup = list(colorschemedict.values())
        orthogroups = list(dict.fromkeys(orthogroupsdup).keys())
        for k in orthogroups:
            frame_update()
            accessions = []
            for l in list(colorschemedict.keys()):
                if colorschemedict[l] == k:
                    accessions.append(l)
            seqdict2 = {}
            for key in list(seqdict.keys()):
                seqdict2[key.split("|")[-1]] = seqdict[key]
            queryseqs = [">" + acc + "\n" + seqdict2[acc] + "\n" for acc in accessions if acc in seqdict2]
            accessions = [acc for acc in accessions if testaccession(acc) == "y"]
            if len(queryseqs) + len(accessions) < 2:
                continue
            musclegroups.append(k)
            url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&tool="multigeneblast"&ID=' + ",".join(accessions) + "&rettype=fasta&retmode=text"
            urltry = "n"
            tries = 0
            while urltry == "n":
                try:
                    time.sleep(3)
                    req = urllib.request.Request(url)
                    response = urllib.request.urlopen(req)
                    output = response.read()
                    if len(output) > 5:
                        urltry = "y"
                    if ">" not in output:
                        log("Downloading of FASTA sequences failed")
                        break
                except (IOError,http.client.BadStatusLine,URLError,http.client.HTTPException):
                    tries += 1
                    if tries == 5:
                        break
                    log("Waiting for connection... (4)")
                    time.sleep(60)
            outfile = open("fasta" + os.sep + "orthogroup" + str(k) + ".fasta","w")
            for seq in queryseqs:
                outfile.write(seq)
            outfile.write(output)
            outfile.close()
            args = "-quiet -in fasta" + os.sep + "orthogroup" + str(k) + ".fasta -out fasta" + os.sep + "orthogroup" + str(k) + "_muscle.fasta"
            muscleprocess = Process(target=runmuscle, args=[args])
            muscleprocess.start()
            while True:
                processrunning = "n"
                if muscleprocess.is_alive():
                    processrunning = "y"
                if processrunning == "y":
                    frame_update()
                elif processrunning == "y":
                    pass
                else:
                    break
    return musclegroups

def create_xhtml_template(queryclusterdata, page, pages):
    #Create HTML file with gene cluster info in hidden div tags
    htmlfile = open(MGBPATH + os.sep + "empty.xhtml","r")
    html = htmlfile.read()
    html = html.replace("\r","\n")
    htmlparts = html.split("<SPLIT HERE>")
    htmloutfile = open("displaypage" + str(page) + ".xhtml","w")
    htmloutfile.write(htmlparts[0] + '  displaycblastresults("' + str(page) + '","all")' + htmlparts[1])
    htmloutfile.write("var list=[" + ",".join([str(nr + (page - 1) * 50 + 1) for nr in range(queryclusterdata[1][0])]) + ",'all'];" + htmlparts[2])
    htmloutfile.write("var list=[" + ",".join([str(nr + (page - 1) * 50 + 1) for nr in range(queryclusterdata[1][0])]) + ",'all'];" + htmlparts[3])
    htmloutfile.write('<a class="bigtext"><br/><br/>&nbsp;Results pages: ')
    for pagenr in [pagenr + 1 for pagenr in range(pages)]:
        htmloutfile.write('<a href="displaypage' + str(pagenr) + '.xhtml" class="bigtext">' + str(pagenr) + '</a>')
        if pagenr != pages:
            htmloutfile.write(", ")
        else:
            htmloutfile.write("</a>")
    htmloutfile.write(htmlparts[4])
    return htmloutfile, htmlparts

def write_xhtml_output(htmloutfile, queryclusterdata, clusters, clusterblastpositiondata, nucname, page, screenwidth, blastdetails, mgb_scores, musclegroups, colorschemedict, dbtype):
    #Write ClusterBlast divs with pictures and description pop-up tags
    frame_update()
    htmloutfile.write('<div id="clusterblastview" class="clusterdescr">\n\n')
    #Add menu bar 3
    htmloutfile.write('<div id="bartext3" style="color:#FFFFFF; font-size:1em; position:absolute; z-index:2; top:3px; left:20px;"><b>MultiGeneBlast hits</b></div>')
    htmloutfile.write('<div id="descrbar3" style="position:absolute; z-index:1; top:0px;"><img src="images/bar.png" height="25" width="' + str(int(0.75*screenwidth)) + '"/></div>')
    htmloutfile.write('<div class="help" id="help3" style="position:absolute; z-index:1; top:2px; left:' + str(int(screenwidth * 0.75) - 30) + 'px;"><a href="http://multigeneblast.sourceforge.net/usage.html" target="_blank"><img border="0" src="images/help.png"/></a></div>')
    qclusternr = page
    nrhitclusters = queryclusterdata[1][0]
    hitclusterdata = queryclusterdata[1][1]
    htmloutfile.write('<div id="qcluster' + str(qclusternr) + '">\n<br/><br/>\n<div align="left">\n<form name="clusterform' + str(qclusternr) + '">\n<select name="selection' + str(qclusternr) + '" onchange="javascript:navigate(this);">\n')
    htmloutfile.write('<option value="">Select gene cluster alignment</option>\n')
    for i in range(nrhitclusters):
        cdescription = hitclusterdata[i + 1 + (page - 1) * 50][5][i].replace("&","&amp;")
        if len(cdescription) > 80:
            cdescription = cdescription[:77] + "..."
        htmloutfile.write('<option value="javascript:displaycblastresults(' + str(page) + ',' + str(i+1 + (page - 1) * 50) + ')">' + cdescription + '</option>\n')
    htmloutfile.write('</select>\n</form>\n\n</div>')
    htmloutfile.write('<div style="position:absolute; top:33px; left:' + str(screenwidth*0.625) + 'px;"><img src="images/button.gif" name="button' + str(qclusternr) + '" onclick="javascript:displaybutton(' + str(qclusternr) + ');"/></div>')
    for i in range(nrhitclusters):
        frame_update()
        hitclusterdata = queryclusterdata[1][1]
        queryclustergenes = hitclusterdata[list(hitclusterdata.keys())[0]][3]
        queryclustergenesdetails = hitclusterdata[list(hitclusterdata.keys())[0]][4]
        hitclusternumber =  i + 1 + (page - 1) * 50
        cluster_acc = hitclusterdata[hitclusternumber][6]
        cluster_blastdetails = blastdetails[cluster_acc]
        mgbscore = mgb_scores[cluster_acc][0]
        cumblastscore = mgb_scores[cluster_acc][1]
        hitclustergenes = clusters[cluster_acc][0]
        hitclustergenesdetails = hitclusterdata[hitclusternumber][2]
        relpositiondata = clusterblastpositiondata[str(page) + "_" + str(i + 1 + (page - 1) * 50)]
        qrel_starts = relpositiondata[0][0]
        qrel_ends = relpositiondata[0][1]
        hrel_starts = relpositiondata[1][hitclusternumber][0]
        hrel_ends = relpositiondata[1][hitclusternumber][1]
        strandsbalance = relpositiondata[2][hitclusternumber]
        hstarts = relpositiondata[3][hitclusternumber][0]
        hends = relpositiondata[3][hitclusternumber][1]
        invertedhstarts = [str(100000000 - int(l)) for l in hstarts]
        invertedhends = [str(100000000 - int(l)) for l in hends]
        if strandsbalance < 0:
            hitclustergenes.reverse()
        htmloutfile.write('<div id="hitcluster' + str(qclusternr) + '_' + str(i + 1 + (page - 1) * 50) + '">\n')
        #Load svg and embed it into XHTML
        svglines = open("svg" + os.sep + "clusterblast" + str(qclusternr) + '_' + str(i + 1 + (page - 1) * 50) + ".svg","r").read().split("\n")
        htmloutfile.write("\n" + svglines[0][:-1] + 'id="svg' + str(qclusternr) + '_' + str(i + 1 + (page - 1) * 50) + '" >' + "\n")
        for svgline in svglines[1:]:
            htmloutfile.write(svgline + "\n")
        #Insert gene cluster descriptions
        cgbkdescription = hitclusterdata[i + 1 + (page - 1) * 50][5][i].replace("&","&amp;").replace("\t"," ").partition(" ")[2].partition(" ")[2].split(", whole")[0].split(", complete")[0].split(", partial")[0]
        if len(cgbkdescription) > 90:
            cgbkdescription = cgbkdescription[:87] + "..."
        if testaccession(cluster_acc.rpartition("_")[0]) == "y":
            cdescription = '<a href="http://www.ncbi.nlm.nih.gov/nuccore/' + cluster_acc.rpartition("_")[0] + '" target="_blank"> ' + cluster_acc.rpartition("_")[0] + "</a>" + " : " + cgbkdescription + "&nbsp;&nbsp;&nbsp;&nbsp;Total score: " + mgbscore + "&nbsp;&nbsp;&nbsp;&nbsp; Cumulative Blast bit score: " + cumblastscore
        else:
            cdescription = cluster_acc.rpartition("_")[0] + " : " + cgbkdescription + "&nbsp;&nbsp;&nbsp;&nbsp;Total score: " + mgbscore + "&nbsp;&nbsp;&nbsp;&nbsp; Cumulative Blast bit score: " + cumblastscore
        if len(nucname) < 90:
            qdescription = "Query: " + nucname
        else:
            qdescription = "Query: " + nucname[0:87] + "..."
        htmloutfile.write('<div id="descriptionquery" style="text-align:left; position:absolute; top:60px; left:10px; font-size:10px; font-style:italic">' + qdescription + '</div>\n')
        htmloutfile.write('<div id="description' + str(qclusternr) + '" style="text-align:left; position:absolute; top:115px; left:10px; font-size:10px; font-style:italic">' + cdescription + '</div>\n')
        #Insert NCBI links
        htmloutfile.write('<div id="pub_pics" style="position:absolute; top:175px; left:' + str(int(screenwidth * 0.0)) + 'px; font-size:10px"> Hit cluster cross-links: \n')
        htmloutfile.write('&nbsp;&nbsp;<a href="http://www.ncbi.nlm.nih.gov/nuccore/' + cluster_acc.rpartition("_")[0] + '" target="_blank"><img align="absmiddle" border="0" src="images/genbank.gif"/></a>\n')
        htmloutfile.write('</div>\n\n')
        #Create gene pop-ups
        a = 0
        for j in queryclustergenes:
            j_accession = j
            htmloutfile.write('<div id="q' + str(qclusternr) + "_" + str(hitclusternumber) + "_" + str(a) + '_div" class="hidden popup" style="position:absolute; top:' + str(100) + 'px; left:' + str(int(float(qrel_starts[a])*0.875)) + 'px;">\n')
            htmloutfile.write(queryclustergenesdetails[j][3].replace("_"," ").replace("&","&amp;") + "\n")
            link = "http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins&amp;PROGRAM=blastp&amp;BLAST_PROGRAMS=blastp&amp;QUERY=" + j_accession + "&amp;LINK_LOC=protein&amp;PAGE_TYPE=BlastSearch"
            if j != queryclustergenesdetails[j][4] and testaccession(j) == "y":
                htmloutfile.write('<br/>Accession: <a href="http://www.ncbi.nlm.nih.gov/protein/' + j + '" target="_blank">' + j + "</a>\n")
            htmloutfile.write("<br/>Location: " + str(queryclustergenesdetails[j][0]) + "-" + str(queryclustergenesdetails[j][1]) + "\n")
            htmloutfile.write("</div>\n\n")
            htmloutfile.write('<div id="q' + str(qclusternr) + "_" + str(hitclusternumber) + "_" + str(a) + '_divtext" class="hidden genenames" style="position:absolute; top:' + str(75) + 'px; left:' + str(int(float((float(qrel_starts[a])+float(qrel_ends[a]))/2)*0.9375)) + 'px;">\n')
            if queryclustergenesdetails[j][4] != "" and queryclustergenesdetails[j][4] != "no_locus_tag":
                htmloutfile.write(queryclustergenesdetails[j][4])
            else:
                htmloutfile.write(j)
            htmloutfile.write("</div>\n\n")
            a+= 1
        a = 0
        for j in hitclustergenes:
            if ((hitclustergenesdetails[j][0] in hstarts or hitclustergenesdetails[j][0] in hends) and (hitclustergenesdetails[j][1] in hends or hitclustergenesdetails[j][1] in hstarts)) or ((hitclustergenesdetails[j][1] in invertedhstarts or hitclustergenesdetails[j][1] in invertedhends) and (hitclustergenesdetails[j][0] in invertedhends or hitclustergenesdetails[j][0] in invertedhstarts)):
                j_accession = j
                htmloutfile.write('<div id="h' + str(qclusternr) + "_" + str(hitclusternumber) + "_" + str(a) + '_div" class="hidden popup" style="position:absolute; top:' + str(151) + 'px; left:' + str(int(float(hrel_starts[a])*0.875)) + 'px;">\n')
                htmloutfile.write(hitclustergenesdetails[j][3].replace("_"," ").replace("&","&amp;") + "\n")
                link = "http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins&amp;PROGRAM=blastp&amp;BLAST_PROGRAMS=blastp&amp;QUERY=" + j_accession + "&amp;LINK_LOC=protein&amp;PAGE_TYPE=BlastSearch"
                if dbtype == "nucl":
                    htmloutfile.write('<br/>Accession: <a href="http://www.ncbi.nlm.nih.gov/nuccore/' + j.rpartition("_")[0] + '" target="_blank">' + j.rpartition("_")[0] + "</a>\n")
                elif j != hitclustergenesdetails[j][4] and testaccession(j) == "y":
                    htmloutfile.write('<br/>Accession: <a href="http://www.ncbi.nlm.nih.gov/protein/' + j + '" target="_blank">' + j + "</a>\n")
                htmloutfile.write("<br/>Location: " + str(hitclustergenesdetails[j][0]) + "-" + str(hitclustergenesdetails[j][1]) + "\n")
                if j in cluster_blastdetails:
                    for blasthit in cluster_blastdetails[j]:
                        htmloutfile.write("<br/><br/><b>BlastP hit with " + blasthit[0] + "</b>\n<br/>Percentage identity: " + blasthit[1] + " %\n")
                        htmloutfile.write("<br/>BlastP bit score: " + blasthit[2] + "\n<br/>Sequence coverage: " + blasthit[3].partition(".")[0] + " %\n")
                        htmloutfile.write("<br/>E-value: " + blasthit[4] + "\n<br/>")
                if testaccession(j) == "y" and dbtype != "nucl":
                    htmloutfile.write("<br/><a href=\"" + link + "\" target=\"_blank\"> NCBI BlastP on this gene </a>\n")
                if j in colorschemedict and colorschemedict[j] in musclegroups:
                    htmloutfile.write("<br/><a href=\"fasta" + os.sep + "orthogroup" + str(colorschemedict[j]) + "_muscle.fasta\" target=\"_blank\"> Muscle alignment of this gene with homologs </a>\n")
                htmloutfile.write("</div>\n\n")
                htmloutfile.write('<div id="h' + str(qclusternr) + "_" + str(hitclusternumber) + "_" + str(a) + '_divtext" class="hidden genenames" style="position:absolute; top:' + str(126) + 'px; left:' + str(int(float((float(hrel_starts[a])+float(hrel_ends[a]))/2)*0.9375)) + 'px;">\n')
                if hitclustergenesdetails[j][4] != "" and hitclustergenesdetails[j][4] != "no_locus_tag":
                    htmloutfile.write(hitclustergenesdetails[j][4])
                else:
                    htmloutfile.write(j)
                htmloutfile.write("</div>\n\n")
                a += 1
        htmloutfile.write('</div>\n')
    #Find new relative positions for display of all gene clusters in one picture
    relpositiondata = clusterblastpositiondata[str(page) + "_all"]
    if len(relpositiondata[0]) > 0:
        qrel_starts = relpositiondata[0][0]
        qrel_ends = relpositiondata[0][1]
        htmloutfile.write('<div id="hitcluster' + str(page) + '_all" style="display:none">\n')
        #Load svg and embed it into XHTML
        svglines = open("svg" + os.sep + "clusterblast" + str(qclusternr) + "_all.svg","r").read().split("\n")
        htmloutfile.write("\n" + svglines[0][:-1] + 'id="svg' + str(qclusternr) + '_all" >' + "\n")
        for svgline in svglines[1:]:
            htmloutfile.write(svgline + "\n")
        if len(nucname) < 90:
            qdescription = "Query: " + nucname
        else:
            qdescription = "Query: " + nucname[0:87] + "..."
        htmloutfile.write('<div id="descriptionquery" style="text-align:left; position:absolute; top:60px; left:10px; font-size:10px; font-style:italic">' + qdescription + '</div>\n')
        for i in range(nrhitclusters):
            frame_update()
            hitclusterdata = queryclusterdata[1][1]
            queryclustergenes = hitclusterdata[list(hitclusterdata.keys())[0]][3]
            queryclustergenesdetails = hitclusterdata[list(hitclusterdata.keys())[0]][4]
            hitclusternumber =  i + 1 + (page - 1) * 50
            hrel_starts = relpositiondata[1][hitclusternumber][0]
            hrel_ends = relpositiondata[1][hitclusternumber][1]
            cluster_acc = hitclusterdata[hitclusternumber][6]
            cluster_blastdetails = blastdetails[cluster_acc]
            mgbscore = mgb_scores[cluster_acc][0]
            cumblastscore = mgb_scores[cluster_acc][1]
            hitclustergenes = clusters[cluster_acc][0]
            hitclustergenesdetails = hitclusterdata[hitclusternumber][2]
            strandsbalance = relpositiondata[2][hitclusternumber]
            hstarts = relpositiondata[3][hitclusternumber][0]
            hends = relpositiondata[3][hitclusternumber][1]
            invertedhstarts = [str(100000000 - int(l)) for l in hstarts]
            invertedhends = [str(100000000 - int(l)) for l in hends]
            cgbkdescription = hitclusterdata[i + 1 + (page - 1) * 50][5][i].replace("&","&amp;").replace("\t"," ").partition(" ")[2].partition(" ")[2].split(", whole")[0].split(", complete")[0].split(", partial")[0]
            if len(cgbkdescription) > 90:
                cgbkdescription = cgbkdescription[:87] + "..."
            if testaccession(cluster_acc.rpartition("_")[0]) == "y":
                cdescription = str(i+1 + (page - 1) * 50) + ". : " + '<a href="http://www.ncbi.nlm.nih.gov/nuccore/' + cluster_acc.rpartition("_")[0] + '" target="_blank"> ' + cluster_acc.rpartition("_")[0] + "</a> " + cgbkdescription + "&nbsp;&nbsp;&nbsp;&nbsp; Total score: " + mgbscore + "&nbsp;&nbsp;&nbsp;&nbsp; Cumulative Blast bit score: " + cumblastscore
            else:
                cdescription = str(i+1 + (page - 1) * 50) + ". : " + cluster_acc.rpartition("_")[0] + " " + cgbkdescription + "&nbsp;&nbsp;&nbsp;&nbsp; Total score: " + mgbscore + "&nbsp;&nbsp;&nbsp;&nbsp; Cumulative Blast bit score: " + cumblastscore
            htmloutfile.write('<div id="description' + str(qclusternr) + '" style="text-align:left; position:absolute; top:' + str(int(63 + (51.7 * (hitclusternumber - (page - 1) * 50)))) + 'px; left:10px; font-size:10px; font-style:italic">' + cdescription + '</div>\n')
            if hitclusternumber == 1 + (page - 1) * 50:
                a = 0
                for j in queryclustergenes:
                    htmloutfile.write('<div id="all_' + str(qclusternr) + "_0_" + str(a) + '_div" class="hidden popup" style="position:absolute; top:' + str(100) + 'px; left:' + str(int(float(qrel_starts[a])*0.875)) + 'px; z-index:2;">\n')
                    htmloutfile.write(queryclustergenesdetails[j][3].replace("_"," ").replace("&","&amp;") + "\n")
                    link = "http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins&amp;PROGRAM=blastp&amp;BLAST_PROGRAMS=blastp&amp;QUERY=" + j + "&amp;LINK_LOC=protein&amp;PAGE_TYPE=BlastSearch"
                    if j != queryclustergenesdetails[j][4] and testaccession(j) == "y":
                        htmloutfile.write('<br/>Accession: <a href="http://www.ncbi.nlm.nih.gov/protein/' + j + '" target="_blank">' + j + "</a>\n")
                    htmloutfile.write("<br/>Location: " + str(queryclustergenesdetails[j][0]) + "-" + str(queryclustergenesdetails[j][1]) + "\n")
                    if testaccession(j) == "y":
                        htmloutfile.write("<br/><a href=\"" + link + "\" target=\"_blank\"> NCBI BlastP on this gene </a>\n")
                    if j in colorschemedict and colorschemedict[j] in musclegroups:
                        htmloutfile.write("<br/><a href=\"fasta" + os.sep + "orthogroup" + str(colorschemedict[j]) + "_muscle.fasta\" target=\"_blank\"> Muscle alignment of this gene with homologs </a>\n")
                    htmloutfile.write("</div>\n\n")
                    htmloutfile.write('<div id="all_' + str(qclusternr) + "_0_" + str(a) + '_divtext" class="hidden genenames" style="position:absolute; top:' + str(75) + 'px; left:' + str(int(float((float(qrel_starts[a])+float(qrel_ends[a]))/2)*0.9375)) + 'px;">\n')
                    if queryclustergenesdetails[j][4] != "" and queryclustergenesdetails[j][4] != "no_locus_tag":
                        htmloutfile.write(queryclustergenesdetails[j][4])
                    else:
                        htmloutfile.write(j)
                    htmloutfile.write("</div>\n\n")
                    a+= 1
            a = 0
            for j in hitclustergenes:
                if ((hitclustergenesdetails[j][0] in hstarts or hitclustergenesdetails[j][0] in hends) and (hitclustergenesdetails[j][1] in hends or hitclustergenesdetails[j][1] in hstarts)) or ((hitclustergenesdetails[j][1] in invertedhstarts or hitclustergenesdetails[j][1] in invertedhends) and (hitclustergenesdetails[j][0] in invertedhends or hitclustergenesdetails[j][0] in invertedhstarts)):
                    htmloutfile.write('<div id="all_' + str(qclusternr) + "_" + str(hitclusternumber) + "_" + str(a) + '_div" class="hidden popup" style="position:absolute; top:' + str(int(100 + 51.7 * (hitclusternumber - (page - 1) * 50))) + 'px; left:' + str(int(float(hrel_starts[a])*0.875)) + 'px; z-index:2;">\n')
                    htmloutfile.write(hitclustergenesdetails[j][3].replace("_"," ").replace("&","&amp;") + "\n")
                    link = "http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins&amp;PROGRAM=blastp&amp;BLAST_PROGRAMS=blastp&amp;QUERY=" + j + "&amp;LINK_LOC=protein&amp;PAGE_TYPE=BlastSearch"
                    if dbtype == "nucl":
                        htmloutfile.write('<br/>Accession: <a href="http://www.ncbi.nlm.nih.gov/nuccore/' + j.rpartition("_")[2] + '" target="_blank">' + j.rpartition("_")[2] + "</a>\n")
                    elif j != hitclustergenesdetails[j][4] and testaccession(j) == "y":
                        htmloutfile.write('<br/>Accession: <a href="http://www.ncbi.nlm.nih.gov/protein/' + j + '" target="_blank">' + j + "</a>\n")
                    htmloutfile.write("<br/>Location: " + str(hitclustergenesdetails[j][0]) + "-" + str(hitclustergenesdetails[j][1]) + "\n")
                    if j in cluster_blastdetails:
                        for blasthit in cluster_blastdetails[j]:
                            htmloutfile.write("<br/><br/><b>BlastP hit with " + blasthit[0] + "</b>\n<br/>Percentage identity: " + blasthit[1] + " %\n")
                            htmloutfile.write("<br/>BlastP bit score: " + blasthit[2] + "\n<br/>Sequence coverage: " + blasthit[3].partition(".")[0] + " %\n")
                            htmloutfile.write("<br/>E-value: " + blasthit[4] + "\n<br/>")
                    if testaccession(j) == "y" and dbtype != "nucl":
                        htmloutfile.write("<br/><a href=\"" + link + "\" target=\"_blank\"> NCBI BlastP on this gene </a>\n")
                    if j in colorschemedict and colorschemedict[j] in musclegroups:
                        htmloutfile.write("<br/><a href=\"fasta" + os.sep + "orthogroup" + str(colorschemedict[j]) + "_muscle.fasta\" target=\"_blank\"> Muscle alignment of this gene with homologs </a>\n")
                    htmloutfile.write("</div>\n\n")
                    htmloutfile.write('<div id="all_' + str(qclusternr) + "_" + str(hitclusternumber) + "_" + str(a) + '_divtext" class="hidden genenames" style="position:absolute; top:' + str(int(75 + 51.7 * (hitclusternumber - (page - 1) * 50))) + 'px; left:' + str(int(float((float(hrel_starts[a])+float(hrel_ends[a]))/2)*0.9375)) + 'px;">\n')
                    if hitclustergenesdetails[j][4] != "" and hitclustergenesdetails[j][4] != "no_locus_tag":
                        htmloutfile.write(hitclustergenesdetails[j][4])
                    else:
                        htmloutfile.write(j)
                    htmloutfile.write("</div>\n\n")
                    a += 1
        htmloutfile.write('</div>\n')
        htmloutfile.write('</div>\n\n')
    else:
        htmloutfile.write('<br/>No homologous gene clusters found.</div>\n')
    htmloutfile.write('</div>\n')
    htmloutfile.write('<div id="creditsbar' + str(i) + '" class="banner" style="position:absolute; width:' + str(int(0.98 * screenwidth)) +'px; align:\'left\'; height:75; top:2750px; left:0px; color:#000066; z-index:-1;">')
    htmloutfile.write('<div style="float:center; font-size:0.9em;">\n<div style="position:absolute; top:0px; left:30px;">\n<img src="images/ruglogo.gif" border="0"/>&nbsp;&nbsp;&nbsp;&nbsp;\n<img src="images/gbblogo.gif" border="0"/>&nbsp;&nbsp;&nbsp;&nbsp;\n</div>\n<div style="position:absolute; top:10px; left:340px;">\nDetecting sequence homology at the gene cluster level with MultiGeneBlast.\n<br/>Marnix H. Medema, Rainer Breitling &amp; Eriko Takano (2013)\n<br/><i>Molecular Biology and Evolution</i> , 30: 1218-1223.\n</div>\n</div>\n</div>')

def finalize_xhtml(htmloutfile, htmlparts):
    #Add final part of HTML file
    htmloutfile.write(htmlparts[-1])
    #Copy accessory files for HTML viewing
    #if sys.platform == ('win32'):
    #  copycommand1 = "copy/y vis\\* " + genomename + " > nul"
    #  copycommand2 = "copy/y vis\\html\\* " + genomename + "\\html > nul"
    #  copycommand3 = "copy/y vis\\images\\* " + genomename + "\\images > nul"
    #elif sys.platform == ('linux2'):
    #  copycommand1 = "cp vis/* " + genomename + " > /dev/null"
    #  copycommand2 = "cp -r vis/html " + genomename + "/html > /dev/null"
    #  copycommand3 = "cp -r vis/images " + genomename + "/images > /dev/null"
    #os.system(copycommand1)
    #os.system(copycommand2)
    #os.system(copycommand3)

    #Close open html file
    htmloutfile.close()

def create_xhtml_file(queryclusterdata, clusters, clusterblastpositiondata, nucname, page, pages, screenwidth, blastdetails, mgb_scores, musclegroups, colorschemedict, dbtype):
    htmloutfile, htmlparts = create_xhtml_template(queryclusterdata, page, pages)
    write_xhtml_output(htmloutfile, queryclusterdata, clusters, clusterblastpositiondata, nucname, page, screenwidth, blastdetails, mgb_scores, musclegroups, colorschemedict, dbtype)
    finalize_xhtml(htmloutfile, htmlparts)

def move_outputfiles(foldername, pages):
    global MGBPATH
    #Move output files to specified output folder. Overwriting when files/folders are already present there
    try:
        os.mkdir(foldername)
    except:
        pass
    try:
        shutil.rmtree(foldername + os.sep + "svg")
    except:
        pass
    try:
        shutil.rmtree(foldername + os.sep + "fasta")
    except:
        pass
    for page in range(pages):
        try:
            os.remove(foldername + os.sep + "displaypage" + str(page + 1) + ".xhtml")
        except:
            pass
        try:
            shutil.move("displaypage" + str(page + 1) + ".xhtml", foldername + os.sep + "displaypage" + str(page + 1) + ".xhtml")
        except:
            pass
    filestomove = ["clusterblast_output.txt", "svg", "fasta"]
    frame_update()
    for f in filestomove:
        try:
            os.remove(foldername + os.sep + f)
        except:
            try:
                shutil.rmtree(foldername + os.sep + f)
            except:
                pass
        try:
            shutil.move(f, foldername + os.sep + f)
        except:
            pass
    frame_update()
    filestocopy = ["style.css", "jquery.svg.js", "jquery-1.4.2.min.js", "jquery.svgdom.js"]
    for f in filestocopy:
        try:
            os.remove(foldername + os.sep + f)
        except:
            pass
        shutil.copy(MGBPATH + os.sep + f, foldername + os.sep + f)
    folderstocopy = ["images"]
    for f in folderstocopy:
        try:
            shutil.rmtree(foldername + os.sep + f)
        except:
            pass
        shutil.copytree(MGBPATH + os.sep + f, foldername + os.sep + f)


def main():
    #if the GUI is active or not.
    global GUI
    #path to temporary file directory. This is for windows a Temp directory on the c drive
    global TEMP
    #path to the multiGeneBlast source files
    global MGBPATH
    os.environ['BLASTDB'] = MGBPATH

    os.chdir(TEMP)
    GUI = "n"
    starttime = time.time()

    #Step 1: parse options into an Option Object
    user_options = get_arguments()

    #configure a logger to track what is happening over multiple files, make sure to do this after the
    #option parsing to save the log at the appropriate place
    setup_logger(user_options.outdir, starttime)
    logging.info("Step 1/11: Input has been parsed")

    #Step 2: Read GBK / EMBL file, select genes from requested region and output FASTA file
    query_proteins = read_query_file(user_options)
    logging.info("Step 2/11: Query input file has been read and parsed")

    #Step 3: Run internal BLAST
    query_cluster = internal_blast(user_options, query_proteins)
    logging.info("Step 3/11: Finished running internal blast")

    #Step 4: Run BLAST on genbank_mf database
    blast_output = db_blast(query_proteins, user_options)
    logging.info("Step 4/11: Finished running blast on the specified database.")

    #Step 5: Parse BLAST output
    blast_dict = parse_db_blast(user_options, query_proteins, blast_output)
    logging.info("Step 5/11: Finished parsing database blast")

    #Step 6: Load genomic databases into memory
    database = load_databases(query_proteins, blast_dict, user_options)
    logging.info("Step 6/11: Finished loading the relevant parts of the database.")

    #Step 7: Locate blast hits in genome and find clusters
    clusters = find_gene_clusters(blast_dict, user_options, database)
    logging.info("Step 7/11: Finished finding clusters.")
    #print([len(c.proteins) for c in clusters])

    #Step 8: Score Blast output on all loci
    score_clusters(clusters, query_proteins, user_options)
    logging.info("Step 8/11: All clusters have been scored.")

    #step 9: write the results to a text file
    write_txt_output(query_proteins, clusters, blast_output, user_options)
    logging.info("Step 9/11: Results have been written to a text file.")

    #Output. From here, iterate for every page
    # for page in [pagenr + 1 for pagenr in range(int(opts.pages))]:
    page = 1
    #Step 9: Write MultiGeneBlast SVGs
    write_svg_files(clusters, user_options, query_cluster, blast_dict)
    print(("Step 9/11, page " + str(page) + ": Time since start: " + str((time.time() - starttime))))
    return
    #Step 10: Create muscle alignments
    musclegroups = align_muscle(opts.muscle, colorschemedict, seqdict)
    print(("Step 10/11, page " + str(page) + ": Time since start: " + str((time.time() - starttime))))

    #Step 11: Create XHTML output file
    create_xhtml_file(queryclusterdata, clusters, clusterblastpositiondata, nucname, page, int(opts.pages), opts.screenwidth, blastdetails, mgb_scores, musclegroups, colorschemedict, opts.dbtype)
    print(("Step 11/11, page " + str(page) + ": Time since start: " + str((time.time() - starttime))))

    #Move all files to specified output folder
    move_outputfiles(opts.outputfolder, int(opts.pages))

    #Close log file
    print(("MultiGeneBlast successfully finished in " + str((time.time() - starttime)) + " seconds.\n"))

def setup_logger(outdir, starttime):
    """
    Function for setting up a logger that will write output to file as well as
    to sdout.

    :param outdir: output directory specified by the user. this is where the
    log file should end up.
    :param starttime: the time the program was started. The logger is created
    slightly later.
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
    logging.basicConfig(filename=log_file_loc, level=logging.DEBUG, format='%(levelname)s: %(asctime)s - %(message)s')

    #configure a handler that formats the logged events properly and prints the events to file as well as stdout
    handler = logging.StreamHandler(sys.stdout)
    formatter = MyFormatter('%(currentTime)s (%(passedTime)s sec) - %(message)s', starttime=starttime)
    formatter
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

if __name__ == '__main__':
    freeze_support()
    main()