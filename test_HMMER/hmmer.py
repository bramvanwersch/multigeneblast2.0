#!/usr/bin/python3

# import statements
import os
import subprocess
import urllib.request
import urllib.error

from Bio import SearchIO

def check_pfam_db(path):
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

def main():
    # this variable now for testing, change when implementing in main code #
    file_path = "D:/Uni/Thesis_MultiGeneBlast/Pfam-A.hmm.gz"

    # step 1: check if pfam db is present and if not download it
    check_pfam_db(file_path)

    # Step 2: Run hmmsearch
    # 2a: Run with pfam accession name: use hmmfetch to get information
    # 2b: Run with own created hmm profile

    # step 3: parse results with biopython


if __name__ == '__main__':
    main()