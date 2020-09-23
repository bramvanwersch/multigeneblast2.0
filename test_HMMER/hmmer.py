#!/usr/bin/python3

# import statements
import os
import sys
import urllib.request
import urllib.error
import argparse
from Bio import SearchIO

# needed this or else utilities import did not work
sys.path.insert(0, "/mnt/d/Uni/Thesis_MultiGeneBlast/multigeneblast2.0")
from utilities import *


def check_pfam_db(path):
    """Check f Pfam-A db exists else download

    :param path: String, path where to check
    """
    url = "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam33.1/Pfam-A.hmm.gz"
    if os.path.exists(path):
        print("Pfam database found")
    else:
        print("Fetching database from Pfam release: 33.1 ")
        try:
            urllib.request.urlretrieve(url, path)
            print("Saved:", path)
        except FileNotFoundError:
            print("Error: Path or file does not exists")
        except urllib.error.URLError or urllib.error.HTTPError:
            print("Error: Internet connection problem")


def fetch_profiles(keys, db):
    """Fetch hmm profiles from db and save in a file

    :param keys: String, Path to file with accession numbers
    :param db: String, path to Pfam db
    """
    print("Fetching profiles from Pfam-A file")
    command_fetch_profile = "hmmfetch -o key.hmm -f {} {}".format(db, keys)
    try:
        run_commandline_command(command_fetch_profile, max_retries=1)
    except MultiGeneBlastException:
        raise MultiGeneBlastException("Key not found in file")
    print("done fetching profiles")


def main():
    # this is for testing, change when implementing in main code #
    pfam_db = "/mnt/d/Uni/Thesis_MultiGeneBlast/Pfam-A.hmm.gz"

    # step 1: check if pfam db is present and if not download it
    check_pfam_db(pfam_db)

    # Step 2: Run hmmsearch
    # 2a: Run with pfam accession name(s): use hmmfetch to get information from
    # Pfam-A db.
    key_file = "/mnt/d/Uni/Thesis_MultiGeneBlast/key_file.txt"
    fetch_profiles(key_file, pfam_db)

    # 2b: Run with own created hmm profile

    # step 3: parse results with biopython


if __name__ == '__main__':
    main()