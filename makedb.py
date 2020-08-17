#!/usr/bin/env python
## Copyright (c) 2012 Marnix H. Medema
## University of Groningen
## Department of Microbial Physiology / Groningen Bioinformatics Centre
## License: GNU General Public License v3 or later
## A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

import os
import sys
import pickle as pickle
import argparse
import logging
import time


from dblib.parse_gbk import parse_gbk_embl, fix_wgs_master_record, fix_supercontig_record
from dblib.utils import fasta_names, make_blast_db, clean_up, get_accession
from dblib.generate_genecords_tar import generate_genecords_tar
from dblib.generate_proteininfo_tar import generate_proteininfo_tar

from genbank_parsing import DataBase
from constants import *
from utilities import MultiGeneBlastException, remove_illegal_characters, setup_logger, run_commandline_command


def get_arguments():
    parser = argparse.ArgumentParser(description='Create a database for offline '
                                                 'search in MultiGeneBlast.')
    parser.add_argument("-i", "-in", help="GBK/EMBL files for creating the "
                                          "database", nargs='+'
                        , metavar="space-separted file paths", required=True,
                        type=check_in_file)
    parser.add_argument("-n", "-name", help="The name for the database files.",
                        required=True, type=remove_illegal_characters)
    parser.add_argument("-o", "-out", help="Optional output folder for the "
                                           "database files.", type=check_out_folder,
                        default=CURRENTDIR + os.sep + "database")
    name_space = parser.parse_args()
    return name_space.n, name_space.o, name_space.i

def check_in_file(in_file):
    """
    Check if all provided input files have the right extensions

    :param files: a list of file locations
    :return: the original list of files
    """
    root, ext = os.path.splitext(in_file)
    if not os.path.exists(in_file):
        raise argparse.ArgumentTypeError("No valid path provided for {}".format(in_file))
    elif not os.path.isfile(in_file):
        raise argparse.ArgumentTypeError("{} is not a file".format(in_file))
    elif ext.lower() not in EMBL_EXTENSIONS + GENBANK_EXTENSIONS:
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
    #make sure relatively defined paths are also correct
    my_path = os.path.abspath(os.path.dirname(__file__))
    path = os.path.join(my_path, path)

    to_path, folder_name = os.path.split(path)
    if not os.path.exists(to_path):
        raise argparse.ArgumentTypeError("No valid output directory"
                                         " provided")
    #assert that the folder name does not contain illegal characters
    elif not folder_name.replace("_", "").isalnum():
        raise argparse.ArgumentTypeError("Output directory name is illegal.")
    return path

def write_nal_pal_file(dbname, outdir, dbtype = "prot"):
    logging.info("Writing nal/pal file...")
    ext = ".pal"
    if dbtype == "nucl":
        ext = ".nal"
    with open(outdir + os.sep + dbname + ext, "w") as f:
        f.write("TITLE " + dbname + "\nDBLIST " + dbname + "\n")


def main():

    starttime = time.time()
    #parse options
    dbname, outdir, inputfiles = get_arguments()

    #setup a logger
    setup_logger(outdir, starttime)
    logging.info("Step 1/6: Parsed options.")

    #create a database object
    base_path = os.path.abspath(os.path.dirname(__file__))
    db = DataBase(base_path, inputfiles)
    db.create(outdir, dbname)
    logging.info("Step 2/6: Created the MultiGeneBlast database")

    #write the database to fasta for makeblastdb
    logging.info("Writing MultiGeneBlast database to fasta...")
    with open(outdir + os.sep + dbname + "_dbbuild.fasta", "w") as f:
        f.write(db.get_fasta())
    logging.info("Step 3/6: Written MultiGeneBlast database to fasta format.")

    #Create Blast database, set the outdir as the current directory to amke sure that the files end up in the right place
    logging.info("Creating Blast+ database...")
    os.chdir(outdir)
    command = "{}\\exec_new\\makeblastdb.exe -dbtype prot -out {} -in {}{}{}_dbbuild.fasta".format(base_path, dbname, outdir, os.sep, dbname)
    run_commandline_command(command, max_retries=0)
    logging.info("Step 4/6: Blast+ database created.")

    #write nal/pal file as main entry for the database
    write_nal_pal_file(dbname, outdir)
    logging.info("Step 5/6: Created nal/pal file.")

    #cleaning up the fasta file
    os.remove(outdir + os.sep + dbname + "_dbbuild.fasta")
    logging.info("Step 6/6: Cleaned up left over files.")
    logging.info("Database was succesfully created.")


if __name__ == "__main__":
    main()
