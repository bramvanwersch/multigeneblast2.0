#!/usr/bin/env python3

"""
Functions for testing multigeneblast with a number of different scenarios that test all the functionality of
multigeneblast

Original creator: Marnix Medena
Recent contributor: Bram van Wersch

Copyright (c) 2012 Marnix H. Medema
License: GNU General Public License v3 or later
A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.
"""


import os
import shutil

from constants import get_mgb_path

MGBPATH = get_mgb_path()


def test_multigeneblast(command, outfolder, remove=True):
    assert os.system(command) == 0

    expected_files = ['run.log', 'jquery-1.4.2.min.js', 'jquery.svg.js', 'jquery.svgdom.js', 'style.css']
    my_path = os.path.abspath(os.path.dirname(__file__))
    db_folder = os.path.join(my_path, "test_data{}{}".format(os.sep, outfolder))

    # check if all files are present
    for path, subdirs, files in os.walk(db_folder):
        # ignore checking the image directory
        if os.path.basename(path) == "images":
            continue
        for file in files:
            if os.path.basename(path) == "visual":
                assert file in expected_files or (file.startswith("displaypage") and file.endswith(".xhtml"))
            elif file.endswith(".mgb"):
                assert file.endswith("_cluster_text.mgb")
            else:
                assert file in expected_files

    # remove the folder
    if remove:
        shutil.rmtree(db_folder)


if __name__ == "__main__":
    # command using gb query
    print("Starting run with genbank file:")
    command1 = "{}{}multigeneblast.py -in test_data{}test_querys{}V_maris_query.gb -from 0 -to 59320 -out test_data{}" \
               "test_run_gbk -db test_data{}test_data_bases{}v_maris.dmnd -c 4 -hpg 200 -msc 25 -mpi 30 -dkb 20 -sw " \
               "0.5 -m n -op 3 -inf all".format(MGBPATH, os.sep, os.sep, os.sep, os.sep, os.sep, os.sep)
    test_multigeneblast(command1, "test_run_gbk")
    print("Finished running gb query test...")
    print()

    # command using embl query
    print("Starting run with embl file:")
    command2 = "{}{}multigeneblast.py -in test_data{}test_querys{}V_maris_query.embl -from 0 -to 59320 -out test_data" \
               "{}test_run_embl -db test_data{}test_data_bases{}v_maris.dmnd -c 4 -hpg 200 -msc 25 -mpi 30 -dkb 20 " \
               "-sw 0.5 -m n -op 3 -inf all".format(MGBPATH, os.sep, os.sep, os.sep, os.sep, os.sep, os.sep)
    test_multigeneblast(command2, "test_run_embl")
    print("Finished running embl query test...")
    print()

    # command using architecture search fasta
    print("Starting run with architecture mode(fasta file):")
    command3 = "{}{}multigeneblast.py -in test_data{}test_querys{}V_maris_query_partial.fasta -out test_data{}test_" \
               "run_fasta -db test_data{}test_data_bases{}v_maris.dmnd -c 4 -hpg 200 -msc 25 -mpi 30 -dkb 20 -sw 0.5" \
               " -m n -op 3 -inf all".format(MGBPATH, os.sep, os.sep, os.sep, os.sep, os.sep, os.sep)
    test_multigeneblast(command3, "test_run_fasta")
    print("Finished running architecture mode test...")
    print()

    # command using gb query and raw nucleotide database
    print("Starting run with genbank file and raw nucleotide database:")
    command4 = "{}{}multigeneblast.py -in test_data{}test_querys{}V_maris_query.gb -from 0 -to 59320 -out test_" \
               "data{}test_run_gbk -db test_data{}test_data_bases{}v_maris_nuc.nal -c 4 -hpg 200 -msc 25 -mpi 30 " \
               "-dkb 20 -sw 0.5 -m n -op 3 -inf all".format(MGBPATH, os.sep, os.sep, os.sep, os.sep, os.sep, os.sep)
    test_multigeneblast(command4, "test_run_gbk")
    print("Finished running gb query test...")
    print()

    # command using embl query and raw nucleotide database
    print("Starting run with embl file and raw nucleotide database:")
    command5 = "{}{}multigeneblast.py -in test_data{}test_querys{}V_maris_query.embl -from 0 -to 59320 -out test_" \
               "data{}test_run_embl -db test_data{}test_data_bases{}v_maris_nuc.nal -c 4 -hpg 200 -msc 25 -mpi 30" \
               " -dkb 20 -sw 0.5 -m n -op 3 -inf all".format(MGBPATH, os.sep, os.sep, os.sep, os.sep, os.sep, os.sep)
    test_multigeneblast(command5, "test_run_embl")
    print("Finished running embl query test...")
    print()

    # command using architecture search fasta and raw nucleotide database
    print("Starting run with architecture mode(fasta file) and raw nucleotide database:")
    command6 = "{}{}multigeneblast.py -in test_data{}test_querys{}V_maris_query_partial.fasta -out test_data{}test_" \
               "run_fasta -db test_data{}test_data_bases{}v_maris_nuc.nal -c 4 -hpg 200 -msc 25 -mpi 30 -dkb 20" \
               " -sw 0.5 -m n -op 3 -inf all".format(MGBPATH, os.sep, os.sep, os.sep, os.sep, os.sep, os.sep)
    test_multigeneblast(command6, "test_run_fasta")
    print("Finished running architecture mode test...")
    print()

    # command testing on a raw nucleotide database made from fasta embl and genbank files
    print("Starting run with raw nucleotide database made from fasta embl and genbank files:")
    command7 = "{}{}multigeneblast.py -in test_data{}test_querys{}V_maris_query.gb -out test_data{}test_run_fasta" \
               " -from 0 -to 59320 -db test_data{}test_data_bases{}mix.nal -c 4 -hpg 200 -msc 25 -mpi 30 -dkb 20" \
               " -sw 0.5 -m n -op 3 -inf all".format(MGBPATH, os.sep, os.sep, os.sep, os.sep, os.sep, os.sep)
    test_multigeneblast(command7, "test_run_fasta")
    print("Finished mixed test...")
    print()

    print("7/7 test_data succesfully finished.")
