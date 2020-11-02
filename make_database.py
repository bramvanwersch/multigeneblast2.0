#!/usr/bin/env python3

"""
Script that can be used to make databases that are compatible for searches with blast

Original creator: Marnix Medena
Recent contributor: Bram van Wersch
"""

import argparse
import logging
import time

from databases import ProteinDataBase, NucleotideDataBase
from constants import *
from utilities import remove_illegal_characters, setup_logger, run_commandline_command,\
    setup_temp_folder

MGBPATH = get_mgb_path()
EXEC = get_exec_folder()


def get_arguments():
    """
    Get arguments supplied by the user using argparse

    :return: a list of three values that are the user defined options
    """
    parser = argparse.ArgumentParser(description='Create a database for offline '
                                                 'search in MultiGeneBlast.')
    parser.add_argument("-i", "-in", help="GBK/EMBL files for creating the database", nargs='+',
                        metavar="space-separted file paths", required=True, type=check_in_file)
    parser.add_argument("-n", "-name", help="The name for the database files.", required=True,
                        type=remove_illegal_characters)
    parser.add_argument("-o", "-out", help="Optional output folder for the database files.", type=check_out_folder,
                        default=CURRENTDIR + os.sep + "database")
    parser.add_argument("-t", "-dbtype", help="The type of the database, protein or nucleotide",
                        choices=["prot", "nucl"], default="prot")
    parser.add_argument("-inf", "-information", help="What level of information should be printed while "
                                                     "running the program (default: basic).",
                        choices=("none", "basic", "all"), default="basic")
    name_space = parser.parse_args()
    return name_space.n, name_space.o, name_space.i, name_space.t, name_space.inf


def check_in_file(in_file):
    """
    Check if a provided input files has the right extensions

    :param in_file: a string path to a file
    :return: the original file
    :raises ArgumentTypeError: when the file does not exist or is not a file or
    has not the correct extension.
    """
    my_path = os.path.abspath(os.path.dirname(__file__))
    in_file = os.path.join(my_path, in_file)
    root, ext = os.path.splitext(in_file)
    if not os.path.exists(in_file):
        raise argparse.ArgumentTypeError("No valid path provided for {}".format(in_file))
    elif not os.path.isfile(in_file):
        raise argparse.ArgumentTypeError("{} is not a file".format(in_file))
    elif ext.lower() not in EMBL_EXTENSIONS + GENBANK_EXTENSIONS + FASTA_EXTENSIONS:
        raise argparse.ArgumentTypeError("No correct extension, {}, for input file {}.".format(ext, root))
    return in_file


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


def clean_outdir(dbname, outdir, dbtype="prot"):
    """
    Clean the output directory if files that have the same name as files that
    are going to be added

    :param dbname: the name of the database
    :param outdir: the directory where the database files need to go
    :param dbtype: the type of the database, either prot or nucl
    """
    files = os.listdir(outdir)

    if dbtype == "prot":
        database_alias_file = dbname + ".pal"
    else:
        database_alias_file = dbname + ".nal"
    remove_files = [dbname + ext for ext in PROT_DATABASE_EXTENSIONS]
    remove_files.append(database_alias_file)
    for file in files:
        if file in remove_files:
            os.remove(outdir + os.sep + file)
            logging.debug("Removed {} because it is going to be created.".format(file))


def write_nal_pal_file(dbname, outdir, dbtype):
    """
    Write a database alias file so Blast+ programs can acces it.

    :param dbname: a string as name of the database
    :param outdir: a string as path to the output directory
    :param dbtype: the type of the database, 'prot' or 'nucl'
    """
    logging.info("Writing nal/pal file...")
    ext = ".pal"
    if dbtype == "nucl":
        ext = ".nal"
    with open(outdir + os.sep + dbname + ext, "w") as f:
        f.write("TITLE " + dbname + "\nDBLIST " + dbname + "\n")


def main():
    """
    Main function called when running the program
    """

    starttime = time.time()

    # create temp multigeneblast temp folder and set the working directory to that folder
    setup_temp_folder()

    # parse options
    dbname, outdir, inputfiles, db_type, log_level = get_arguments()

    # setup a logger
    setup_logger(outdir, starttime, level=log_level)
    logging.info("Step 1/6: Parsed options.")

    # make sure to clear all files in the destination folder that have the same
    # name as files that are going to be added
    clean_outdir(dbname, outdir)
    logging.info("Step 2/6: Cleaned the output directory of potential duplicate files.")

    # create a database object
    base_path = os.path.abspath(os.path.dirname(__file__))
    if db_type == "prot":
        db = ProteinDataBase(base_path, inputfiles)
    else:
        db = NucleotideDataBase(base_path, inputfiles)
    db.create(outdir, dbname)
    logging.info("Step 3/6: Created the MultiGeneBlast database")

    # write the database to fasta for makeblastdb
    logging.info("Writing MultiGeneBlast database to fasta...")
    with open("{}{}{}_dbbuild.fasta".format(TEMP, os.sep, dbname), "w") as f:
        f.write(db.get_fasta())
    logging.info("Step 4/6: Written MultiGeneBlast database to fasta format.")

    # Create Blast database, set the outdir as the current directory to amke sure that the files end up in the
    # right place
    logging.info("Creating Blast+ database...")

    # make sure to change the working directory for the makeblastdb files to edn up in the right place
    os.chdir(outdir)
    command = "{}{}makeblastdb -dbtype {} -out {} -in {}{}{}_dbbuild.fasta".format(EXEC, os.sep, db_type, dbname, TEMP,
                                                                                   os.sep, dbname)

    run_commandline_command(command, max_retries=0)
    logging.info("Step 5/6: Blast+ database created.")

    # write nal/pal file as main entry for the database
    write_nal_pal_file(dbname, outdir, db_type)
    logging.info("Step 6/6: Created nal/pal file.")

    logging.info("Database was succesfully created at {}.".format(outdir))


if __name__ == "__main__":
    main()
