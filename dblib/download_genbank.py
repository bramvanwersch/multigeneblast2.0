#!/usr/bin/env python
## Copyright (c) 2012 Marnix H. Medema
## University of Groningen
## Department of Microbial Physiology / Groningen Bioinformatics Centre
## License: GNU General Public License v3 or later
## A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

#Script that retrieves GenBank database from NCBI
from ftplib import FTP
import os
import sys
import gzip
from .utils import log

def handleDownload(block):
    genbank_download_file.write(block)

def read_in_chunks(file_object, chunk_size=1024):
    """Lazy function (generator) to read a file piece by piece.
    Default chunk size: 1k."""
    while True:
        data = file_object.read(chunk_size)
        if not data:
            break
        yield data

def extract_file(gzfile, GUI="n"):
  print("Extracting " + gzfile)
  log("Extracting " + gzfile)
  try:
    gzipfile = gzip.open(gzfile,"rb")
    outfile = open(gzfile.rpartition(".gz")[0],"wb")
    for chunk in read_in_chunks(gzipfile):
      outfile.write(chunk)
    outfile.close()
    gzipfile.close()
  except:
    if GUI == "n":
      print("Extraction of " + gzfile + "failed. Press enter to continue; press 'Q' to quit.")
      log("Extraction of " + gzfile + "failed. Press enter to continue; press 'Q' to quit.")
      user_response = input()
      if user_response.upper() == "Q":
        sys.exit(1)
    else:
      log("Extraction of " + gzfile + "failed. Exiting", exit=True)

def parse_options(dbname, args):
  #Parse options
  gbktypes = ["BCT", "PLN", "PHG", "SYN", "PAT", "ENV", "CON", "WGS"]
  chosengbktypes = []
  if len(args) < 1:
    print("""
    InputError: No arguments provided.

    Please supply a database name and one or more GBK types as arguments:
    BCT   bacterial sequences
    PLN   plant/fungal/algal sequences
    PHG   bacteriophage sequences
    SYN   synthetic sequences
    PAT   patent sequences
    ENV   environmental sequences
    WGS   whole-genome shotgun sequences
    CON   assembly data

    Example: python download_genbank.py <database name> BCT PAT WGS\n
    """)
    default = input("Use default GBK types BCT, CON, WGS and default database name 'genbank_mf'? (y/n)")
    while default != "y" and default != "n":
      default = input("Use default GBK types BCT, CON, WGS? (y/n)")
    if default == "y":
      chosengbktypes = ["BCT", "CON", "WGS"]
      dbname = "genbank_mf"
    elif default == "n":
      sys.exit(1)
  else:
    for arg in args:
      if arg.upper() not in gbktypes:
        print("Error: '" + arg + "' is not a GenBank sequence type.\n")
        print("""
        MAKEGBDB usage:

        Please supply a database name and one or more GBK types as arguments:
        BCT   bacterial sequences
        PLN   plant/fungal/algal sequences
        PHG   bacteriophage sequences
        SYN   synthetic sequences
        PAT   patent sequences
        ENV   environmental sequences
        WGS   whole-genome shotgun sequences
        CON   assembly data

        Example: python makegbdb.py <database name> BCT PAT WGS
        """)
        sys.exit(1)
      else:
        chosengbktypes.append(arg.upper())
  return gbktypes, chosengbktypes, dbname

def create_directories(chosengbktypes):
  #Creates local genbank / wgs directories if it doesn't exist
  try:
    os.mkdir("genbank")
  except:
    pass
  if "WGS" in chosengbktypes:
    try:
      os.mkdir("wgs")
    except:
      pass

def ftp_connect():
  #Connect to NCBI FTP site and change to genbank directory
  try:
    ftp = FTP('ftp.ncbi.nlm.nih.gov')   # connect to host, default port
    ftp.login()               # user anonymous, passwd anonymous@
    ftp.cwd("genbank")
    print("Connected succesfully to ftp.ncbi.nlm.nih.gov")
    log("Connected succesfully to ftp.ncbi.nlm.nih.gov")
    return ftp
  except:
    try:
      ftp = FTP('bio-mirror.net')   # connect to host, default port
      ftp.login()               # user anonymous, passwd anonymous@
      ftp.cwd("biomirror")
      ftp.cwd("genbank")
      print("Connection to ftp.ncbi.nlm.nih.gov failed; connected succesfully to bio-mirror.net")
      log("Connection to ftp.ncbi.nlm.nih.gov failed; connected succesfully to bio-mirror.net")
      return ftp
    except:
      print("Connection to ftp servers failed. Please check your internet connection. Exiting.")
      log("Connection to ftp servers failed. Please check your internet connection. Exiting.")
      sys.exit(1)

def download_general_genbank_files(ftp, chosengbktypes, GUI="n"):
  #Find all basic genbank files; download all basic GenBank entries selected
  dirlines = []
  downloadedfiles = []
  if chosengbktypes != ["WGS"]:
    print("Getting directory contents...")
    log("Getting directory contents...")
    try:
      ftp.retrlines('LIST',dirlines.append)     # list directory contents
    except:
      print("Getting directory contents failed. Exiting...")
      log("Getting directory contents failed. Exiting...")
      sys.exit(1)
    gzfiles = [item for sublist in [dirline.split(" ") for dirline in dirlines] for item in sublist if ".gz" in item]
    for gzfile in gzfiles:
      for gbktype in chosengbktypes:
        if "gb" + gbktype.lower() in gzfile:
          print("Downloading " + gzfile)
          log("Downloading " + gzfile)
          try:
            global genbank_download_file
            genbank_download_file = open("genbank/" + gzfile,"wb")
            ftp.retrbinary("RETR " + gzfile, handleDownload)
            downloadedfiles.append(gzfile)
            genbank_download_file.close()
          except:
            try:
              disconnect(ftp)
              ftp = ftp_connect()
              genbank_download_file = open("genbank/" + gzfile,"wb")
              ftp.retrbinary("RETR " + gzfile, handleDownload)
              downloadedfiles.append(gzfile)
              genbank_download_file.close()
            except:
              if GUI == "n":
                print("Downloading of " + gzfile + "failed. Press enter to continue; press 'Q' to quit.")
                log("Downloading of " + gzfile + "failed. Press enter to continue; press 'Q' to quit.")
                user_response = input()
                if user_response.upper() == "Q":
                  sys.exit(1)
              else:
                log("Downloading of " + gzfile + "failed.", exit=True)

def download_wgs_genbank_files(ftp, chosengbktypes, GUI="n", dbtype="prot"):
  dirlines = []
  downloadedwgsfiles = []
  #Find and download WGS entries
  if "WGS" in chosengbktypes:
    ftp.cwd("wgs")
    print("Getting directory contents...")
    log("Getting directory contents...")
    try:
      ftp.retrlines('LIST',dirlines.append)     # list directory contents
    except:
      print("Getting directory contents failed. Exiting...")
      log("Getting directory contents failed. Exiting...")
      sys.exit(1)
    if dbtype == "prot":
        gnpgzfiles = [item for sublist in [dirline.split(" ") for dirline in dirlines] for item in sublist if (".gnp.gz" in item or "mstr.gbff.gz" in item)]
    else:
        gnpgzfiles = [item for sublist in [dirline.split(" ") for dirline in dirlines] for item in sublist if (".fsa_nt.gz" in item)]
    for gnpgzfile in gnpgzfiles:
      if "wgs" in gnpgzfile:
        print("Downloading " + gnpgzfile)
        log("Downloading " + gnpgzfile)
        try:
          global genbank_download_file
          genbank_download_file = open("wgs/" + gnpgzfile,"wb")
          ftp.retrbinary("RETR " + gnpgzfile, handleDownload)
          downloadedwgsfiles.append(gnpgzfile)
          genbank_download_file.close()
        except:
          if GUI == "n":
            print("Downloading of " + gnpgzfile + "failed. Press enter to continue; press 'Q' to quit.")
            log("Downloading of " + gnpgzfile + "failed. Press enter to continue; press 'Q' to quit.")
            user_response = input()
            if user_response.upper() == "Q":
              sys.exit(1)
          else:
            log("Downloading of " + gnpgzfile + "failed.", exit=True)

def disconnect(ftp):
  print("Finished downloading and extracting GenBank files. Closing connection...")
  log("Finished downloading and extracting GenBank files. Closing connection...")
  #Disconnects from FTP server
  try:
    ftp.quit()
  except:
    ftp.close()

def find_downloaded_files(chosengbktypes):
  downloadedgbkfiles = [gbkfile for gbkfile in os.listdir("genbank") if ".seq" in gbkfile]
  if "WGS" in chosengbktypes:
    downloadedwgsfiles = [wgsfile for wgsfile in os.listdir("wgs") if (".gnp" in wgsfile or ".gbff" in wgsfile)]
  else:
    downloadedwgsfiles = []
  return downloadedgbkfiles, downloadedwgsfiles

def download_genbank_files(dbname, arguments, GUI="n", dbtype="prot"):
  gbktypes, chosengbktypes, dbname = parse_options(dbname, arguments)
  create_directories(chosengbktypes)
  #downloadedgbkfiles, downloadedwgsfiles = find_downloaded_files(chosengbktypes) #May be used to implement a --continue function after a broken-off attempt
  ftp = ftp_connect()
  download_general_genbank_files(ftp, chosengbktypes, GUI)
  download_wgs_genbank_files(ftp, chosengbktypes, GUI, dbtype)
  disconnect(ftp)

if __name__ == "__main__":
    download_genbank_files(sys.argv)