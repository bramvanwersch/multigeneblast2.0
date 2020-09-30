#!/usr/bin/python3

# import statements
import argparse
import gzip
import os
import sys
import urllib.request
import urllib.error

from Bio import SearchIO

# needed this or else utilities import did not work
sys.path.insert(0, "/mnt/d/Uni/Thesis_MultiGeneBlast/multigeneblast2.0")
from utilities import *


def check_pfam_db(path, file_names):
    """Check f Pfam-A db exists else download

    :param path: String, path where to check
    :param file_names: list of strings, names of file in a list
    """
    url_ls = ["ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam33.1/Pfam-A.hmm.gz",
              "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam33.1/Pfam-A.hmm.dat.gz"]
    if os.path.exists(path + "Pfam-A.hmm.gz"):
        print("Pfam database found")
    else:
        print("Fetching database from Pfam release: 33.1 ")
        counter = 0
        for url in url_ls:
            try:
                urllib.request.urlretrieve(url, path + file_names[counter])
            except FileNotFoundError:
                print("Error: Path or file does not exists")
            except urllib.error.URLError or urllib.error.HTTPError:
                print("Error: Internet connection problem")
            counter += 1


def get_full_accession_number(keys_file, db_path):
    """Get full accession number of Pfam profiles

    :param keys_file: String, Path to file with acc-nr
    :param db_path: String, Path to dat.gz file with the full acc-nr
    :return: key_lines: List, string of full acc-number
    """
    # Read the file with incomplete acc-numbers
    file = open(keys_file, 'r')
    key_lines = file.readlines()
    file.close()
    # Read dat.gz file with complete acc-numbers
    dat_gz_file = gzip.open(db_path + 'Pfam-A.hmm.dat.gz', 'r')
    content = str(dat_gz_file.read()).split("\\n")
    # Select the incomplete ones from the dat.gz info
    # Only appends to list when it is found in dat.gz file
    profile_ls = []
    for text in content:
        for key in key_lines:
            if key.strip() in text:
                profile_ls.append(text.split(" ")[-1])
    return key_lines


def fetch_profiles(keys, db_folder):
    """Fetch hmm profiles from db and save in a file

    """
    print("Fetching profiles from Pfam-A file")
    get_full_accession_number(keys, db_folder)
    #command_fetch_profile = "hmmfetch -o {}\key.hmm -f {} {}".format(db_path,
    #                                            db_path+file_names[0], keys)
    #subprocess.run(command_fetch_profile, shell=True)


def run_hmmsearch(hmm_prof, db):
    """Run hmmsearch of hmmer3

    :param hmm_prof: String, Path to rofile to search
    :param db: String, Path to db in fasta format
    """
    print("Preforming hmmsearch")
    # -o results.txt -> add to save in file
    command_run_hmmsearch = "hmmsearch {} {}".format(hmm_prof, db)
    subprocess.run(command_run_hmmsearch, shell=True)


def main():
    # this is for testing, change when implementing in main code #
    pfam_db_folder = "/mnt/d/Uni/Thesis_MultiGeneBlast/"
    names_file = ["Pfam-A.hmm.gz", "Pfam-A.hmm.dat.gz"]

    # step 1: check if pfam db is present and if not download it
    check_pfam_db(pfam_db_folder, names_file)

    # Step 2: Run hmmsearch
    # 2a: Run with pfam accession name(s): use hmmfetch to get profiles from
    # Pfam-A db.
    key_file = "/mnt/d/Uni/Thesis_MultiGeneBlast/key_file.txt"
    fetch_profiles(key_file, pfam_db_folder)

    # 2b: Run hmmsearch
    path_key = "/mnt/d/Uni/Thesis_MultiGeneBlast/key.hmm"
    path_db = "/mnt/d/Uni/Thesis_MultiGeneBlast/UP000008308_263358.fasta.gz"
    #run_hmmsearch(path_key, path_db)

    # step 3: parse results with biopython package


if __name__ == '__main__':
    main()