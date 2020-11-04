#!/usr/bin/env python3

"""
The functions used to run multigeneblast. Also the file called to run multigeneblast from the command line

Original creator: Marnix Medena
Recent contributor: Bram van Wersch

Copyright (c) 2012 Marnix H. Medema
License: GNU General Public License v3 or later
A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.
"""


# imports
import tarfile
from multiprocessing import Process, freeze_support
import random
import argparse
from collections import OrderedDict
import pickle as pickle
import colorsys

# own imports
from databases import GenbankFile, embl_to_genbank, Protein
from visualisation import ClusterCollectionSvg, create_xhtml_file
from utilities import *
from constants import *

# constant specifically required for mbg
MGBPATH = get_mgb_path()
EXEC = get_exec_folder()


def write_fasta(protein_list, file):
    """
    Write a fasta file using a list of Protein objects to the current directory

    :param protein_list: a list of Protein objects
    :param file: the output file
    """
    try:
        with open(file, "w") as out_file:
            for prot in protein_list:
                out_file.write(prot.fasta_text())
    except Exception:
        logging.critical("No fasta file created for sequences with file name {}".format(file))
        raise MultiGeneBlastException("Cannot open file {}".format(file))
    logging.debug("Saved file {} at {}.".format(file, os.getcwd()))


#### STEP 1: PARSE OPTIONS ####
def get_arguments():
    """
    Parse the command line arguments using argparser.

    :return: an Option object that holds the options specified by the user
    """
    parser = argparse.ArgumentParser(description='Run multigeneblast on a specified database using a specified query.',
                                     epilog="-from, -to and -genes are mutually exclusive")
    parser.add_argument("-i", "-in", help="Query file name: GBK/EMBL file for homology search,FASTA file with multiple"
                                          " protein sequences for architecture search", required=True,
                        type=check_in_file, metavar="file path")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-f", "-from", help="Start position of query region", type=int, metavar="Int")
    group.add_argument("-g", "-genes", help="Accession codes of genes constituting query multigene module", nargs='+',
                       metavar="space-separted gene identifiers")

    parser.add_argument("-t", "-to", help="End position of query region", type=int,
                        metavar="Int")
    parser.add_argument("-o", "-out", help="Output folder path in which results will be stored", type=check_out_folder,
                        metavar="file path")
    parser.add_argument("-db", "-database", help="Diamond database to be queried", required=True, type=check_db_folder,
                        metavar="file path")
    parser.add_argument("-c", "-cores", help="Number of parallel CPUs to use for threading (default: all)",
                        type=determine_cpu_nr, default="all", metavar="Int")
    parser.add_argument("-hpg", "-hitspergene", help="Number of Blast hits per query gene to be taken into account "
                                                     "(default: 250), max = 10000", type=int, choices=range(50, 10001),
                        default=250, metavar="[50 - 10000]")
    parser.add_argument("-msc", "-minseqcov", help="Minimal %% coverage of a Blast hit on hit protein sequence to be "
                                                   "taken into account (default: 25)", type=int, choices=range(0, 101),
                        default=25, metavar="[0 - 100]")
    parser.add_argument("-mpi", "-minpercid", help="Minimal %% identity of a Blast hit on hit protein sequence to be "
                                                   "taken into account (default: 30)", type=int, choices=range(0, 101),
                        default=30, metavar="[0 - 100]")
    parser.add_argument("-dkb", "-distancekb", help="Maximum kb distance between two blast hits to be counted as"
                                                    " belonging to the same locus (default: 10)", default=20,
                        type=int, choices=range(0, 101), metavar="[0 - 100]")
    parser.add_argument("-sw", "-syntenywheight", help="Weight of synteny conservation in hit sorting score: "
                                                       "(default: 0.5)", type=restricted_float, default=0.5,
                        metavar="Float")
    parser.add_argument("-m", "-muscle", help="Generate a Muscle multiple sequence alignments of all hits of each input"
                                              " gene (default: n)", default="n", choices=["y", "n", "Y", "N"])
    parser.add_argument("-op", "-outpages", help="Maximum number of output pages (with 50 hits each) to be generated"
                                                 " (default: 5)", type=int, default=5, choices=range(1, 41),
                        metavar="[1 - 40]")
    parser.add_argument("-inf", "-information", help="What level of information should be printed while running the"
                                                     " program (default: basic).", choices=("none", "basic", "all"),
                        default="basic")

    name_space = parser.parse_args()

    # some final checks for certain arguments, raise error when appropriate
    # when defining a from without to or a to without a from
    if (name_space.f is not None and name_space.t is None) or \
            (name_space.f is None and name_space.t is not None):
        parser.error("When specifying -from also specify -to.")

    # if the -in argument is a fasta file make sure that no -from, -to or genes are defined
    # otherwise make sure that they are defined
    if (any(name_space.i.endswith(ext) for ext in FASTA_EXTENSIONS) and
            (name_space.t is not None or name_space.f is not None or name_space.g is not None)):
        parser.error("When providing a fasta file (running architecture mode) -f, -t and"
                     " -g arguments cannot be specified")
    elif (not any(name_space.i.endswith(ext) for ext in FASTA_EXTENSIONS) and
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

    # make sure relatively defined paths are also correct
    my_path = os.path.abspath(os.path.dirname(__file__))
    path = os.path.join(my_path, path)

    if not os.path.exists(path):
        raise argparse.ArgumentTypeError("File path to query is invalid")

    if not os.path.isfile(path):
        raise argparse.ArgumentTypeError("Provided query file is not a file.")
    root, ext = os.path.splitext(path)
    if ext.lower() not in GENBANK_EXTENSIONS + EMBL_EXTENSIONS + FASTA_EXTENSIONS:
        raise argparse.ArgumentTypeError("Please supply input file with valid GBK / EMBL extension (homology search)"
                                         " or FASTA extension (architecture search).")
    return path


def check_out_folder(path):
    """
    Check the provided output folder

    :param path: a string that is an absolute or relative path
    :raises ArgumentTypeError: when the path provided does not exist or the
    specified name for the folder contains illegal characters
    :return: the original path
    """
    # make sure relatively defined paths are also correct
    my_path = os.path.abspath(os.path.dirname(__file__))
    path = os.path.join(my_path, path)

    to_path, folder_name = os.path.split(path)
    if not os.path.exists(to_path):
        raise argparse.ArgumentTypeError("No valid output directory"
                                         " provided")
    # assert that the folder name does not contain illegal characters
    elif not folder_name.replace("_", "").isalnum():
        raise argparse.ArgumentTypeError("Output directory name is illegal.")
    return path


def check_db_folder(path):
    """
    Check the provided database folder to make sure that all the required files
    are present

    :param path: a string that is an absolute or relative path
    :raises ArgumentTypeError: when the path provided does not exist or the
    database folder does not contain all neccesairy files
    :return: the original path
    """

    # make sure relatively defined paths are also correct
    my_path = os.path.abspath(os.path.dirname(__file__))
    path = os.path.join(my_path, path)

    if not os.path.exists(path):
        raise argparse.ArgumentTypeError("Database path does not exist")
    to_path, db_file = os.path.split(path)

    # make sure the databae file is of the correct file type
    root, ext = os.path.splitext(path)
    dbname = root.split(os.sep)[-1]
    if ext == ".dmnd":
        expected_extensions = DATABASE_EXTENSIONS
    else:
        raise argparse.ArgumentTypeError("Incorrect extension {} for database file."
                                         "Should be .dmnd".format(ext))

    # make sure the db folder contains all files required
    db_folder = os.listdir(to_path)
    expected_files = [dbname + ext for ext in expected_extensions]
    for file in expected_files:
        if file not in db_folder:
            raise argparse.ArgumentTypeError("Expected a file named {} in "
                                             "the database folder".format(file))
    return path


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
        # the input file
        self.infile = arguments.i
        self.run_name = self.__get_run_name()
        # the output directory
        self.outdir = arguments.o
        # the database to use
        self.db_name, self.db, self.dbtype = self.__configure_db_and_type(arguments.db)

        self.architecture_mode = self.__check_architecture_mode()
        # values that define the query
        self.startpos = arguments.f
        self.endpos = arguments.t
        self.ingenes = arguments.g

        self.cores = arguments.c
        self.minseqcov = arguments.msc
        self.minpercid = arguments.mpi
        self.hitspergene = arguments.hpg
        self.distancekb = int(arguments.dkb * 1000)
        self.muscle = arguments.m.lower()
        self.pages = arguments.op
        self.syntenyweight = arguments.sw

        self.log_level = arguments.inf
        self.gui = False

    def __get_run_name(self):
        """
        Extract a name for the run based on the name of the input file

        :return: a string that is the name without extension or path
        """
        root, file = os.path.split(self.infile)
        file_name = file.split(".", 1)[0]
        return file_name

    def __check_architecture_mode(self):
        """
        Check if MultiGeneBlast is running in architecture or homolgy search mode

        :return: Boolean that is True if the input file ends in a fasta extension
        meaning MultiGeneBlast is running in architecture mode.
        """
        if any(self.infile.lower().endswith(ext) for ext in FASTA_EXTENSIONS):
            return True
        return False

    def __configure_db_and_type(self, database_path):
        to_path, db_file = os.path.split(database_path)

        # make sure the databae file is of the correct file type
        dbname, ext = os.path.splitext(db_file)
        if dbname.endswith("_nuc"):
            dbtype = "nucl"
        else:
            dbtype = "prot"
        return dbname, database_path, dbtype


#### STEP 2: READ INPUT FILE ####

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
        logging.debug("Started parsing genbank query file...")
        gb_file = GenbankFile(user_options.infile, protein_range=[user_options.startpos, user_options.endpos],
                              allowed_proteins=user_options.ingenes)
        query_proteins = gb_file.proteins
        logging.debug("Finished parsing.")
    else:
        logging.debug("Started parsing embl query file")
        genbank_file = embl_to_genbank(user_options.infile)
        gb_file = GenbankFile(user_options.infile, file_text=genbank_file,
                              protein_range=[user_options.startpos, user_options.endpos],
                              allowed_proteins=user_options.ingenes)
        query_proteins = gb_file.proteins
        logging.debug("Finished parsing.")

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

    # read the fasta file with headers that do not contain forbidden characters
    fasta_entries = fasta_to_dict(fasta_file)

    # make for each fasta sequence a protein object and save them in a dictionary
    total_lenght = 0
    query_proteins = {}
    for entry_name in fasta_entries:
        sequence = fasta_entries[entry_name]
        entry_protein = Protein(sequence, total_lenght, entry_name, "+")
        query_proteins[entry_protein.name] = entry_protein

        # +100 for drawing purposes to show the proteins at appropriate places
        total_lenght += entry_protein.nt_lenght + 100
    logging.debug("Finished parsing architecture input file")
    return query_proteins


#### STEP 3/4/5: RUNNING AND PARSING BLAST ####
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

    # Make Blast db for using the sequences saved in query.fasta in the previous step
    make_blast_db_command = "{}{}diamond makedb --in query.fasta --db {}{}query_db".format(EXEC, os.sep, TEMP, os.sep)
    logging.debug("Started making internal blast database...")
    try:
        run_commandline_command(make_blast_db_command, max_retries=5)
    except MultiGeneBlastException:
        raise MultiGeneBlastException("Data base can not be created. The input query is most likely invalid.")
    logging.debug("Finished making internal blast database.")

    # Run and parse BLAST search

    blast_search_command = "{}{}diamond blastp --db {}{}query_db --query {}{}query.fasta --outfmt 6" \
                           " --max-target-seqs 1000 --evalue 1e-05 --out {}{}internal_input.out" \
                           " --threads {}".format(EXEC, os.sep, TEMP, os.sep, TEMP, os.sep, TEMP, os.sep, user_options.cores)

    logging.debug("Running internal blastp...")
    try:
        run_commandline_command(blast_search_command, max_retries=5)
    except MultiGeneBlastException:
        raise MultiGeneBlastException("Running blast when finding internal homologs returns multiple errors.")
    logging.debug("Finished running internal blastp.")

    try:
        with open("{}{}internal_input.out".format(TEMP, os.sep), "r") as f:
            blastoutput = f.read()
    except Exception:
        logging.critical("Something went wrong reading the blast output file. Exiting...")
        raise MultiGeneBlastException("Something went wrong reading the blast output file.")
    # extract blast hits into a dictionary
    blast_lines = blastoutput.split("\n")[:-1]
    blast_dict = blast_parse(user_options, query_proteins, blast_lines)

    # because query proteins are sorted this is fine
    first_protein = list(query_proteins.values())[0]
    last_protein = list(query_proteins.values())[-1]
    query_cluster = Cluster(first_protein.start, last_protein.stop, first_protein.contig_id,
                            first_protein.contig_description)
    for results in blast_dict.values():
        for blast_result in results:
            query_cluster.add_protein(query_proteins[blast_result.subject], blast_result)
    return query_cluster


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

        # protein object that is the query
        self.query_protein = query_protein

        self.query = tabs[0]
        self.subject = tabs[1]
        self.percent_identity = float(tabs[2])
        self.allignment_lenght = int(tabs[3])
        # number of mismatches
        self.mismatch = int(tabs[4])
        # number of gap openings
        self.gapopen = int(tabs[5])

        # allignment in query
        self.query_start = int(tabs[6])
        self.query_stop = int(tabs[7])

        # allignment in subject
        self.subject_start = int(tabs[8])
        self.subject_stop = int(tabs[9])

        # expected value
        self.evalue = float(tabs[10])
        self.bit_score = float(tabs[11])

        # calculated values
        # percentage of the query that is covered by the subject
        self.percent_coverage = self.allignment_lenght / self.query_protein.aa_lenght * 100

    def summary(self):
        """
        The original blast line reported for the blast result

        :return: a string with tab separated values
        """
        return "\t".join(self.line)


def blast_parse(user_options, query_proteins, blast_lines):
    """
    Parse the output of a blast search using BlastResult objects.

    :param user_options: a Option object containing the options supplied by the user
    :param query_proteins: a dictionary that links gene_names against Protein objects
    :param blast_lines: the output produced by NCBI BLAST+ blastp algorithm version
    2.2.18 split on newlines
    :return: a dictionary that contains a key for every query that had a siginificant
    blast result with the values being BlastResult objects that have the key as
    query and are above user defined treshholds.
    """
    logging.debug("Started parsing blast output...")

    # Filter for unique blast comparissons
    query_subject_combinations = set()
    unique_blast_results = {}
    for line in blast_lines:
        tabs = line.split("\t")
        query = tabs[0]
        subject = tabs[1]
        query_subject_combination = query + subject
        if not (query_subject_combination in query_subject_combinations):
            query_subject_combinations.add(query_subject_combination)
            blast_result = BlastResult(tabs, query_proteins[query])
            unique_blast_results[blast_result.query + blast_result.subject] = blast_result

    # Filters blastlines to get rid of hits that do not meet criteria and save
    # hits that do in a dictionary that links a protein name to BlastHit objects.
    internal_homolog_dict = {}
    count = 0
    # TODO remove the or statement, now here for consistency sake. This statement is in the old version and removing it
    #  removes a number of clusters. It needs to go but will stay for now to keep comparissons accurate
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


def db_blast(user_options):
    """
    Run blast agianst the user defined database

    :param user_options: an Option object with user options
    :return: the blast output
    """
    logging.info("Running NCBI BLAST+ searches on the provided database {}..".format(user_options.db_name))
    command_start = "{}{}diamond blastp".format(EXEC, os.sep)

    complete_command = "{} --db {} --query {}{}query.fasta --outfmt 6 --max-target-seqs" \
                       " {} --evalue 1e-05 --out {}{}input.out --threads {}"\
        .format(command_start, user_options.db, os.getcwd(), os.sep, user_options.hitspergene,
                os.getcwd(), os.sep, user_options.cores)

    logging.debug("Started blasting against the provided database...")
    db_blast_process = Process(target=run_commandline_command, args=[complete_command, 0])
    db_blast_process.start()
    start_time = time.time()
    last_message_time = start_time
    while db_blast_process.is_alive():
        time_passed = time.time() - last_message_time
        if time_passed > 60:
            last_message_time = time.time()
            logging.info("{:.1f} minutes passed. Blasting still in progress."
                         .format((last_message_time - start_time) / 60))
    # when the process exits with an error
    if db_blast_process.exitcode != 0:
        raise MultiGeneBlastException("Blasting against the provided database returns multiple errors.")
    logging.debug("Finished blasting")

    try:
        with open("{}{}input.out".format(TEMP, os.sep), "r") as f:
            blastoutput = f.read()
    except Exception:
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
    blast_lines = blast_output.split("\n")[:-1]
    blast_dict = blast_parse(user_options, query_proteins, blast_lines)
    if len(blast_dict) == 0:
        logging.info("No blast hits encountered against the provided database."
                     " Maybe try and lower min_percentage_identity ({}) or "
                     "minimum_sequence_coverage ({})"
                     " Extiting...".format(user_options.minpercid,
                                           user_options.minseqcov))
        # make sure to exit with a non zero code for when using the gui
        sys.exit(1)
    return blast_dict


def filter_nuc_output(blast_output):
    """
    Function for removing overlapping blasthits that occur when running tblastn.
    Also renames subjects to be unique, the naming scheme for this is the accession
    of the contig underscore a unique number.

    :param blast_output: a large text containing the output generated by tblastn
    version 2.2.18
    :return: a list of lines, filtered based on the criteria

    Overlapping blast hits are regarded as hits that are within or largely within
    a larger hit. There is a small ALLOWED_OVERLAP to make sure that hits that
    are just overlapping are not left out of the results.

    This function is neccesairy because double hits cannot simply be filtered
    on query subject filtering because multiple allignments against a big contig
    in the raw nucleotide database have the same subject (the big contig).
    """
    logging.debug("Removing overlapping hits in big nucleotide contigs...")
    blast_lines = blast_output.split("\n")[:-1]
    new_blast_lines = []
    count = 0
    while len(blast_lines) > 0:
        compare_line = blast_lines.pop()
        c_tabs = compare_line.split("\t")
        c_start, c_stop = int(c_tabs[8]), int(c_tabs[9])
        if c_start > c_stop:
            c_stop, c_start = int(c_tabs[8]), int(c_tabs[9])
        for line_index in range(len(blast_lines) - 1, - 1, - 1):
            line = blast_lines[line_index]
            tabs = line.split("\t")
            start, stop = int(tabs[8]), int(tabs[9])
            if start > stop:
                stop, start = int(tabs[8]), int(tabs[9])
            # no overlap skip
            if stop - ALLOWED_OVERLAP < c_start or start + ALLOWED_OVERLAP > c_stop:
                continue
            # there is overlap and the the size of the line is bigger replace it
            elif stop - start > c_stop - c_start:
                c_tabs = tabs
                c_start, c_stop = start, stop
            # overlapping but not bigger remove
            else:
                blast_lines.pop(line_index)
        if "|" in c_tabs[1]:
            accession = c_tabs[1].split("|")[1]
        else:
            accession = c_tabs[1]
        subject = "{}_{}".format(accession, count)
        count += 1
        c_tabs[1] = subject
        new_blast_lines.append("\t".join(c_tabs))
    return new_blast_lines


#### STEP 6: LOAD RELEVANT PARTS OF THE DATABASE ####
def load_database(blast_dict, user_options):
    """
    General starting point for loading a database. The protocol for loading
    a protein database and a nucleotide database are quite different.

    :param blast_dict: a dictionary that contains all query genes with a list
    of BlastResults objects that where against these genes.
    :param user_options: Option object of user specified options
    :return: a GenbankFile object containing all the contigs with blast hits.
    """
    logging.info("Loading GenBank positional info into memory...")

    database = load_protein_database(blast_dict, user_options)
    return database


def load_protein_database(blast_dict, user_options):
    """
    Load each contig that has at least one blast hit.

    :param blast_dict: a dictionary that contains all query genes with a list
    of BlastResults objects that where against these genes.
    :param user_options: Option object of user specified options
    :return: a GenbankFile object containing all the contigs with blast hits.

    Finds all contigs that have a blast hit using the index_pickle that saves
    sets of proteins for each contig, then opens these contigs from the tar
    archive and puts them into a GenbankFile object.
    """
    db_path, _ = os.path.split(user_options.db)

    covered_contigs = set()
    with open("{}{}{}_database_index.pickle".format(db_path, os.sep, user_options.db_name), "rb") as index_file:
        list_of_contig_proteins = pickle.load(index_file)
    for results in blast_dict.values():
        for result in results:
            protein = result.subject
            for contig_nr, contig in enumerate(list_of_contig_proteins):
                if contig_nr in covered_contigs:
                    continue
                if protein in contig:
                    covered_contigs.add(contig_nr)
            # all contigs have been selected
            if len(covered_contigs) == len(list_of_contig_proteins):
                break
    contigs = []
    with tarfile.open("{}{}{}_contigs.tar.gz".format(db_path, os.sep, user_options.db_name), "r:gz") as tf:
        for contig_nr in covered_contigs:
            for entry in tf:
                if "{}/{}.pickle".format(user_options.db_name, contig_nr) == entry.name:
                    file_obj = tf.extractfile(entry)
                    contigs.append(pickle.load(file_obj))

    gb_collection = GenbankFile(contigs=contigs)
    return gb_collection


def load_nucleotide_database(blast_dict, user_options):
    """
    Load each contig that has at least one blast hit.

    :param blast_dict: a dictionary that contains all query genes with a list
    of BlastResults objects that where against these genes.
    :param user_options: Option object of user specified options
    :return: a GenbankFile object containing all the contigs with blast hits.

    First sort all blast hits in a dicitonary with contigs as keys, using the
    name of the blast hits to extract the contig. Then create a genbankfile
    text that can be read by a GenbankFile object to create a database with all
    relevant contigs and proteins.
    """
    db_path = os.environ["BLASTDB"]

    with open("{}{}{}_database_index.pickle".format(db_path, os.sep, user_options.db_name), "rb") as index_file:
        list_of_contig_accessions = pickle.load(index_file)

    hits_per_contig = {}
    for results in blast_dict.values():
        for result in results:
            accession = result.subject.split("_")[0]
            if accession in hits_per_contig:
                hits_per_contig[accession].append(result)
            else:
                hits_per_contig[accession] = [result]
    contigs = []
    with tarfile.open("{}{}{}_contigs.tar.gz".format(db_path, os.sep, user_options.db), "r:gz") as tf:
        for accession in hits_per_contig:
            accession_index = list_of_contig_accessions.index(accession)
            for entry in tf:
                if "{}/{}.pickle".format(user_options.db, accession_index) == entry.name:
                    file_obj = tf.extractfile(entry)
                    contigs.append(pickle.load(file_obj))
    file_text = make_genbank_from_nuc_contigs(contigs, hits_per_contig)
    gb_collection = GenbankFile(file_text=file_text)
    return gb_collection


def make_genbank_from_nuc_contigs(contigs, hits_per_contig):
    """
    Create a genbnak file from a dictionary blast hits and a list of
    NucleotideContig objects containing accessions and sequences

    :param contigs: a list of NucleotideContig objects
    :param hits_per_contig: a dictionary linking contig names to blast
    hits that are on the contig
    :return: a genbank file text, contining the minimal information for a
    GenbankFile object to read.
    """
    file_text = ""
    for contig in contigs:
        file_text += "DEFINITION  {}\nACCESSION   {}\n".format(contig.definition.replace("\n", ""), contig.accession)
        for index, hit in enumerate(hits_per_contig[contig.accession]):
            if hit.subject_start > hit.subject_stop:
                file_text += "     CDS             complement({}..{})\n".format(hit.subject_stop, hit.subject_start)
            else:
                file_text += "     CDS             {}..{}\n".format(hit.subject_start, hit.subject_stop)
            file_text += '                     /locus_tag="{}"\n'.format(hit.subject)
            file_text += '                     /codon_start=1\n'
            file_text += '                     /product=tBlastn hit on {}."\n'.format(contig.accession)
            file_text += '                     /protein_id=tBlastn_hit_{}"\n'.format(index)
        file_text += "\nORIGIN      \n"
        file_text += contig.dna_sequence
        file_text += "\n//\n"
    return file_text


#### STEP 7: FIND GENE CLUSTERS ####
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

    # add clusters to one big list and sort them
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
            # extract the subject proteins from the Blastresults using the database
            subject_protein = database.proteins[hit.subject]
            if hit.subject in query_per_blast_hit:
                query_per_blast_hit[hit.subject].append(hit)
            else:
                query_per_blast_hit[hit.subject] = [hit]
            if subject_protein.contig_id not in hits_per_contig:
                hits_per_contig[subject_protein.contig_id] = {subject_protein}
            else:
                hits_per_contig[subject_protein.contig_id].add(subject_protein)
    # sort each contig using start and stop locations
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
    # find all the clusters
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
                c = Cluster(start - extra_distance, protein.stop + extra_distance, protein.contig_id,
                            protein.contig_description)
                # add all the blast_hits to the cluster together with the query
                # gene they originate from
                for hit in blast_hits:
                    for result in query_per_blast_hit[hit.name]:
                        c.add_protein(hit, result)
                clusters.append(c)
                start = next_protein.start
                blast_hits = [next_protein]
            index += 1
        # make sure to add the final cluster
        c = Cluster(start - extra_distance, scaffold_proteins[index].stop + extra_distance,
                    scaffold_proteins[index].contig_id, scaffold_proteins[index].contig_description)
        blast_hits.append(scaffold_proteins[index])
        for hit in blast_hits:
            for result in query_per_blast_hit[hit.name]:
                c.add_protein(hit, result)
        clusters.append(c)

        # all cluster objects are from the same contig
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

        # because all clusters are sorted aswell as the proteins one loop suffices
        # to assign all proteins
        for protein in contig_proteins:
            cluster = clusters[cluster_index]
            if cluster.start <= protein.start <= cluster.stop:
                cluster.add_protein(protein)
                # if a protein is added make sure to add the same protein to a
                # different cluster that can be overlapping in the extra region part
                if cluster_index + 1 < len(clusters):
                    next_cluster = clusters[cluster_index + 1]
                    if next_cluster.start <= protein.stop <= next_cluster.stop:
                        next_cluster.add_protein(protein)

            elif protein.start > cluster.start:
                cluster_index += 1
                # when all clusters have been covered
                if cluster_index >= len(clusters):
                    break
                cluster = clusters[cluster_index]
                if protein.start >= cluster.start <= cluster.stop:
                    cluster.add_protein(protein)


class Cluster:
    """
    Tracks all proteins in a cluster. The proteins should start within the
    cluster bounds. The cluster can be sorted but only when specifically asked
    by the user. The cluster is then guaranteed sorted unitl new proteins are
    added.
    """
    COUNTER = 0

    def __init__(self, start, stop, contig, contig_description):
        """
        :param start: an integer that is the start of the cluster. Is allowed to
        be negative, this just means that the cluster starts from 0
        :param stop: an Integer that is the end of the cluster. Is allowed to be
        bigger then the contig. Simply means that the cluster stops at the end
        of the contig
        :param contig: the name of the contig this cluster is located on.
        """
        # track an increasing cluster number
        self.no = self.COUNTER
        Cluster.COUNTER += 1

        # start and stop are allowed to be negative or bigger then possible
        self.start = start
        self.stop = stop
        self.contig = contig
        self.contig_description = contig_description
        # all proteins in the cluster
        self.__proteins = OrderedDict()
        self.__protein_indexes = OrderedDict()
        self.__sorted = True

        # all protein names that are blast hits linked to the query protein that
        # that the blast hit was against
        self.__blast_proteins = OrderedDict()

        # a dictionary for tracking individual parts of the score
        self.__score = {"synteny": 0, "accumulated blast score": 0, "number unique hits": 0}

    def add_protein(self, protein, query_blast_hit=None):
        """
        Controlled way of adding proteins to the cluster.

        :param protein: a Protein object
        :param query_blast_hit: a potential Protein object that is part of the
        query that defines to what query gene the 'protein' has a blast hit.
        """
        if query_blast_hit is not None:
            if protein.name in self.__blast_proteins:
                self.__blast_proteins[protein.name].add(query_blast_hit)
            else:
                # create a set to auto filter duplicates
                self.__blast_proteins[protein.name] = {query_blast_hit}
        if protein.name not in self.__proteins:
            self.__proteins[protein.name] = protein
            self.__protein_indexes[protein.name] = len(self.__proteins) - 1
            self.__sorted = False

    def get_protein(self, name):
        """
        Get a protein object from the cluster by protein.name

        :param name: a string that is the name property of a protein object
        :return: a protein Object
        """
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
        Retrieve the start and stop coordinate of the cluster. This has to
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
        self.__score["accumulated blast score"] = accumulated_blast_score
        self.__score["number unique hits"] = number_unique_hits

    def segmented_score(self, names="all"):
        """
        Get a selection of the score values using string names of the score values. The default is all the
        values.

        :param names: optional 3 names 'synteny', 'accumulated_blast_score' and 'number unique hits'
        :return: a dictionary with the names as keys and the associated scores as numbers as values
        """
        if names == "all":
            names = ["synteny", "accumulated blast score", "number unique hits"]
        requested_scores = [val for key, val in self.__score.items() if key in names]
        return requested_scores

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
            self.__proteins = OrderedDict(sorted(self.__proteins.items(),
                                                 key=lambda item: (item[1].start, item[1].stop)))
            self.__protein_indexes = OrderedDict((name, index) for index, name in enumerate(self.__proteins))
            self.__blast_proteins = OrderedDict(sorted(self.__blast_proteins.items(),
                                                       key=lambda item: (self.__proteins[item[0]].start,
                                                                         self.__proteins[item[0]].stop)))
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
        # general information
        full_summary += "- Source contig: {}\n".format(self.contig)
        full_summary += "- Number of proteins with BLAST hits to this cluster: {}\n".format(len(self.__blast_proteins))
        full_summary += "- MultiGeneBlast 2.0 score: {:.7f}\n  - synteny score: {}\n  " \
                        "- accumulated blast bit score / 1.000.000: {:.7f}\n  " \
                        "- Number of unique blast hits: {}\n\n".\
            format(self.score, self.__score["synteny"], self.__score["accumulated blast score"],
                   self.__score["number unique hits"])

        # Table of all proteins in cluster
        full_summary += ">Table of genes in cluster:\n"
        full_summary += "Name\tstart\tstop\tstrand\tannotation\tlocus_tag\n"

        for protein in self.__proteins.values():
            full_summary += "{}\n".format(protein.summary())
        full_summary += "<\n"

        # Table of all the blast hits
        full_summary += "\n>Table of all blast hits:\n"
        full_summary += "query prot\tsubject prot\tperc ident\tallignment len" \
                        "\tmismatch\tgap open\tquery start\tquery stop\tsubject" \
                        " start\tsubject stop\te-value\tblast bit score\n"

        for result_set in self.__blast_proteins.values():
            for blast_result in result_set:
                full_summary += "{}\n".format(blast_result.summary())
        full_summary += "<\n"
        return full_summary


#### STEP 8: SCORE THE CLUSTERS ####
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

    # calculate the cumulative blast scores over all blast hits of a cluster
    # TODO figure out why these scores differ slightly from mgb. I have looked but not found :(
    for cluster in clusters:
        cum_blast_score = 0
        hitnr = 0
        hit_positions_dict = {}
        for result_set in cluster.blast_hit_proteins.values():
            for blast_result in result_set:
                cum_blast_score += blast_result.bit_score / 1_000_000.0
                # get the best scoring blast result for each query entry
                if blast_result.query not in hit_positions_dict:
                    hitnr += 1
                    hit_positions_dict[blast_result.query] = blast_result
                else:
                    prev_br = hit_positions_dict[blast_result.query]
                    if prev_br.evalue > blast_result.evalue:
                        hit_positions_dict[blast_result.query] = blast_result

        synteny_score = 0
        if not user_options.architecture_mode:
            # make a tuple of index in the query and index in the cluster
            hit_positions = []
            for query_protein in hit_positions_dict:
                blast_result = hit_positions_dict[query_protein]
                hit_positions.append((index_query_proteins[blast_result.query],
                                     cluster.index(blast_result.subject)))
            hit_positions.sort()
            synteny_score = score_synteny(hit_positions)

        # assign the score
        cluster.set_score(synteny_score * user_options.syntenyweight, cum_blast_score, hitnr)

    # sort the clusters on score
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


#### STEP 9: WRITE TEXT OUTPUT ####
def write_txt_output(query_proteins, clusters, user_options):
    """
    Write the output into a file containing information on all clusters.

    :param query_proteins: a list of Protein objects that where present in the
    provided query
    :param clusters: a list of Cluster objects
    :param user_options: an Option object with user define doptions

    Note the formatting of the file is defined as follows:
    - '>>' signify a header of a new piece of information, the text directly after
    it is the descriptor
    - '>' signifies the start of a table followed by the descriptor of what is in
    the table. The first line of the table should always be the header.
    - '<' indicates the end of a table
    - '-' indicates a piece of information, that can be regarded as general
    facts of a certain header
    """
    file_name = user_options.run_name + TEXT_OUT_NAME
    logging.info("Writing .mgb output file into {}...".format(os.path.join(user_options.outdir, file_name)))
    try:
        out_file = open("{}{}{}".format(user_options.outdir, os.sep, file_name), "w")
    except (OSError, IOError):
        logging.critical("Writing .mgb file has failed. Can not opern file {}"
                         .format(os.path.join(user_options.outdir, file_name)))
        raise MultiGeneBlastException("Can not open the output file {}."
                                      .format(os.path.join(user_options.outdir, file_name)))

    # write for what database
    out_file.write("ClusterBlast scores for: {}\n\n".format(user_options.db))

    # write all the query proteins
    out_file.write(">>Query proteins:\n")
    out_file.write(">Table of query proteins:\n")
    out_file.write("Name\tstart\tstop\tstrand\tannotation\tlocus_tag\n")
    for protein in query_proteins.values():
        out_file.write("{}\n".format(protein.summary()))
    out_file.write("<\n")

    # write general information about the clusters
    out_file.write("\n>>Significant clusters:\n")
    out_file.write(">Table of the summary of significant clusters:\n")
    out_file.write("cluster no\tcontig id\tscore\tnumber of unique blast hits\n")
    for cluster in clusters:
        out_file.write("{}\n".format(cluster.short_summary()))
    out_file.write("<\n\n")

    # write more comprehensive information for each cluster
    for cluster in clusters:
        out_file.write(">> Details Cluster {}\n".format(cluster.no))
        out_file.write("{}\n".format(cluster.summary()))
    out_file.close()


#### STEP 10: alligning muscle ####
def align_muscle(gene_color_dict, database, query_proteins, outdir):
    """
    Create multiple sequence allignments for all blast hits against all query proteins that have blast htis

    :param gene_color_dict: a dictionary linking protein names to rgb colors
    :param database: a GenbankFile object containing all relevant data
    :param query_proteins: a dictionary of Protein objects with all the query proteins
    :param outdir: the output directorty specified by the user
    """
    # create folder to store resutls
    logging.info("Started to allign with muscle. This will take a while...")
    try:
        os.mkdir("muscle_info")
    except(IOError, OSError):
        pass
    try:
        os.mkdir(outdir + os.sep + "muscle_info")
    except(IOError, OSError):
        pass
    color_per_genes = {}
    for gene, color in gene_color_dict.items():
        if color in color_per_genes:
            color_per_genes[color].append(gene)
        else:
            color_per_genes[color] = [gene]
    start_time = time.time()
    last_message_time = start_time
    count = 1
    total = len(color_per_genes)
    for color, color_group in color_per_genes.items():
        color = color.replace("#", "")
        file_name = "{}{}muscle_info" + os.sep + "group{}.fasta".format(TEMP, os.sep, color)
        with open(file_name, "w") as f:
            for prot in color_group:
                try:
                    sequence = database.proteins[prot].sequence
                except KeyError:
                    sequence = query_proteins[prot].sequence
                f.write(">{}\n".format(prot))
                f.write(sequence + "\n")

        complete_command = "{}{}{} -quiet -in {}{}{} -out {}{}muscle_info{}group{}_muscle.fasta".\
            format(EXEC, os.sep, MUSCLE_PROG_NAME, TEMP, os.sep, file_name, outdir, os.sep, os.sep, color)
        muscle_process = Process(target=run_commandline_command, args=[complete_command, 0])
        muscle_process.start()
        while muscle_process.is_alive():
            time_passed = time.time() - last_message_time
            if time_passed > 60:
                last_message_time = time.time()
                logging.info("{:.1f} minutes passed. Muscle alligning still in process."
                             .format((last_message_time - start_time)/60))
        # when the process exits with an error
        if muscle_process.exitcode != 0:
            logging.warning("Running Muscle returns multiple errors. Muscle output cannot be generated as a result.")
        logging.debug("Finished making muscle allignment {}/{}".format(count, total))
        count += 1


#### STEP 11: CREATE SVG VISUAL OUTPUT ####
def write_svg_files(clusters, gene_color_dict, user_options, query_cluster, page_nr):
    """
    Write clusters into svgs of comparissons between the query and each
    individual cluster as well as an svg as a complete overview of all clusters
    on a page

    :param clusters: a list of Cluster objects
    :param gene_color_dict: a dictionary linking protein names to rgb colors
    :param user_options: an option object with user defined options
    :param query_cluster: a Cluster ibject with the user defined query
    :param page_nr: integer of the current page
    :return: a dictionary of ClusterCollectionSvg objects that contain all the
    information to draw the clusters as well as location information to place
    tooltips in the right place
    """
    logging.debug("Writing visualization SVGs and XHTML...")

    svg_images = {}
    for index, cluster in enumerate(clusters):
        svg_clusters = [cluster]
        collection = ClusterCollectionSvg(query_cluster, svg_clusters, gene_color_dict, user_options)
        svg_images["clusterblast_{}".format(cluster.no)] = collection

    # Create svgs for multiple gene cluster alignment
    collection = ClusterCollectionSvg(query_cluster, clusters, gene_color_dict, user_options)
    svg_images["clusterblast_page{}_all".format(page_nr)] = collection

    return svg_images


def calculate_colorgroups(blast_dict, query_cluster):
    """
    Get distinct colors for all the groups of genes present in the query cluster

    :param blast_dict: a dictionary of Blast result objects
    :param query_cluster: a Cluster object containing the user defined query
    :return: a dictionary that links proteins with blast hits to
    colors so all proteins with the same blast subject have the same color
    """
    # first get groups internal homologs in the query cluster
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

    # get all blast proteins associated to the query proteins in the same groups
    color_groups = []
    for st in groups:
        color_group = set()
        for protein_name in st:
            if protein_name not in blast_dict:
                continue
            color_group.add(protein_name)
            for blast_result in blast_dict[protein_name]:
                color_group.add(blast_result.subject)
        if len(color_group) > 0:
            color_groups.append(list(color_group))
    # Generate RGB scheme
    rgb_color_scheme = generate_distinct_colours(len(color_groups))

    # generate a dictionary linking all proteins with blast hits to colors
    gene_color_dict = OrderedDict()
    for gr_indx, group in enumerate(color_groups):
        for protein_name in group:
            gene_color_dict[protein_name] = rgb_color_scheme[gr_indx]
    return gene_color_dict


def generate_distinct_colours(count):
    """
    Create random distinct non-white colors

    :param count: an integer for the amount of requested colors
    :return: a list of colors in #hex format

    Function from antismash:
    https://github.com/antismash/antismash/blob/98eabb4f0597e51608b49ba270af3a7b23bd5b4e/antismash/modules/clusterblast/svg_builder.py#L512
    """
    # make sure that random processes are consistent to allow reproducible colors
    random.seed(1)
    rgbs = [colorsys.hsv_to_rgb(i/count, .9, .85) for i in range(count)]
    colours = []
    for rgb in rgbs:
        red, green, blue = [str(hex(int(i*255)))[2:] for i in rgb]  # [2:] strips the 0x
        colours.append("#{}{}{}".format(red, green, blue))
    random.shuffle(colours)
    return colours


#### MOVE OUTPUT FILES ####
def move_outputfiles(outdir, pages):
    """
    Delete temporary folders and move files to the output folder

    :param outdir: The name output folder, defined  by the user
    :param pages: the amount of pages of visual output produced
    """
    logging.info("Cleaning up the last files and moving output...")

    # create folder for visuals
    try:
        os.mkdir(outdir + os.sep + "visual")
    except Exception:
        pass

    # move the page visuals from temp to outdir
    for page_indx in range(pages):
        page_nr = page_indx + 1
        try:
            os.remove(outdir + os.sep + "visual" + os.sep + "displaypage{}.xhtml".format(page_nr))
        except Exception:
            pass
        try:
            shutil.move("displaypage{}.xhtml".format(page_nr), outdir + os.sep + "visual" + os.sep +
                        "displaypage{}.xhtml".format(page_nr))
        except Exception:
            pass

    # copy visual files and folders for the xhtml page
    filestocopy = ["style.css", "jquery.svg.js", "jquery-1.4.2.min.js", "jquery.svgdom.js"]
    for f in filestocopy:
        # make sure to remove original files
        try:
            os.remove(outdir + os.sep + "visual" + os.sep + f)
        except Exception:
            pass
        try:
            shutil.copy(MGBPATH + os.sep + "visual_copys" + os.sep + f, outdir + os.sep + "visual" + os.sep + f)
        except IOError:
            raise MultiGeneBlastException("{} cannot be copied from {}. The file is missing."
                                          .format(f, MGBPATH + os.sep + "visual_copys"))

    folderstocopy = ["images"]
    for f in folderstocopy:
        # make sure to remove original files
        try:
            shutil.rmtree(outdir + os.sep + "visual" + os.sep + f)
        except Exception:
            pass
        shutil.copytree(MGBPATH + os.sep + "visual_copys" + os.sep + f, outdir + os.sep + "visual" + os.sep + f)


def get_page_sizes(total_clusters, max_pages):
    """
    Get the sizes of each page that is created by the user

    :param total_clusters: an Integer that is the size of the total amount of
    clusters produced
    :param max_pages: the maximum amount of pages allowed by the user
    :return: a list of sizes of each page. The page sizes are HITS_PER_PAGE
    escept for the last page that can be smaller
    """
    page_total, last_page_lenght = divmod(total_clusters, HITS_PER_PAGE)
    if page_total >= max_pages:
        page_total = max_pages - 1
        last_page_lenght = HITS_PER_PAGE
    page_sizes = [HITS_PER_PAGE for _ in range(page_total)] + [last_page_lenght]
    return page_sizes


def main():
    """
    Starting point of the program
    """

    # set the current directory to the temporary directory
    setup_temp_folder()
    starttime = time.time()

    # Step 1: parse options into an Option Object
    user_options = get_arguments()

    # configure a logger to track what is happening over multiple files, make sure to do this after the
    # option parsing to save the log at the appropriate place(outdir)
    setup_logger(user_options.outdir, starttime, user_options.log_level)
    logging.info("Step 1/12: Input has been parsed")

    # Step 2: Read GBK / EMBL file, select genes from requested region and output FASTA file
    query_proteins = read_query_file(user_options)
    logging.info("Step 2/12: Query input file has been read and parsed")

    # Step 3: Run internal BLAST
    query_cluster = internal_blast(user_options, query_proteins)
    logging.info("Step 3/12: Finished running internal blast")

    # Step 4: Run BLAST on genbank_mf database
    blast_output = db_blast(user_options)
    logging.info("Step 4/12: Finished running blast on the specified database.")

    # Step 5: Parse BLAST output
    blast_dict = parse_db_blast(user_options, query_proteins, blast_output)
    logging.info("Step 5/12: Finished parsing database blast")

    # Step 6: Load genomic databases into memory
    database = load_database(blast_dict, user_options)
    logging.info("Step 6/12: Finished loading the relevant parts of the database.")

    # Step 7: Locate blast hits in genome and find clusters
    clusters = find_gene_clusters(blast_dict, user_options, database)
    logging.info("Step 7/12: Finished finding clusters. {} clusters found in total".format(len(clusters)))

    # Step 8: Score Blast output on all loci
    score_clusters(clusters, query_proteins, user_options)
    logging.info("Step 8/12: All clusters have been scored.")

    # step 9: write the results to a text file
    write_txt_output(query_proteins, clusters, user_options)
    logging.info("Step 9/12: Results have been written to the mgb file.")

    logging.info("Creating visual output...")
    page_sizes = get_page_sizes(len(clusters), user_options.pages)

    # get the gene color dict neccesairy for all svgs and potential muscle allignments
    gene_color_dict = calculate_colorgroups(blast_dict, query_cluster)

    if user_options.muscle == "y":
        align_muscle(gene_color_dict, database, query_proteins, user_options.outdir)
    logging.info("Step 10/12: Muscle allignments created (if -muscle 'y')")

    for page_indx, page_size in enumerate(page_sizes):
        logging.debug("Creating visual information for cluster {} to {}".format(page_indx * HITS_PER_PAGE + 1,
                                                                                min((page_indx + 1) *
                                                                                    HITS_PER_PAGE, len(clusters))))
        page_clusters = clusters[page_indx * HITS_PER_PAGE: min((page_indx + 1) * HITS_PER_PAGE, len(clusters))]

        svg_images = write_svg_files(page_clusters, gene_color_dict, user_options, query_cluster, page_indx + 1)
        logging.info("Step 11/12: Visual images {}/{} created.".format(page_indx + 1, len(page_sizes)))

        # Step 11: Create XHTML output file
        create_xhtml_file(page_clusters, query_cluster, user_options, svg_images, page_sizes, page_indx,
                          gene_color_dict)
        logging.info("Step 12/12: XHTML page {}/{} created.".format(page_indx + 1, len(page_sizes)))

    # Move all files to specified output folder
    move_outputfiles(user_options.outdir, len(page_sizes))
    logging.info("MultiGeneBlast succesfully finished. Output can be found at {}.".format(user_options.outdir))


if __name__ == '__main__':
    freeze_support()
    main()
