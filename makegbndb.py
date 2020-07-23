#!/usr/bin/env python
## Copyright (c) 2012 Marnix H. Medema
## University of Groningen
## Department of Microbial Physiology / Groningen Bioinformatics Centre
## License: GNU General Public License v3 or later
## A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

#Main script that allows downloading and reformatting of (parts of the) GenBank database for MultiGeneBlast

import sys
import os

#Find path to mgb files if run from another directory
pathfolders = os.environ['PATH'].split(os.pathsep)
pathfolders.append(os.getcwd())
pathfolders.append(os.getcwd().rpartition(os.sep)[0])
CURRENTDIR = os.getcwd()
MGBPATH = ""
for folder in pathfolders:
  try:
    if "read_input_gui.py" in os.listdir(folder) and "guilib.py" in os.listdir(folder) and "empty.xhtml" in os.listdir(folder) and "multigeneblast.py" in os.listdir(folder) and "mgb_gui.py" in os.listdir(folder):
      MGBPATH = folder
      break
  except:
    pass
if MGBPATH == "":
  print("Error: Please add the MultiGeneBlast installation directory to your $PATH environment variable before running the executable from another folder.")
  sys.exit(1)
#Set other environment variables
os.environ['EXEC'] = MGBPATH + os.sep + "exec"
os.environ['PATH'] = os.environ['EXEC'] + os.pathsep + os.environ['PATH']

import shutil
from multiprocessing import Process, freeze_support
from time import sleep
from dblib.utils import *
from dblib.download_genbank import *
from dblib.parse_seq import parse_nseq
from dblib.parse_gnp import parse_ngnp

def get_genbank_db(dbname, args, GUI="n", dbtype="prot"):
  #Download GenBank files based on arguments; otherwise download default GenBank files
  download_genbank_files(dbname, args, GUI, dbtype)
  
def main():
  freeze_support()
  logfile = open("makedb.log","w")
  args = sys.argv
  if len(args) > 2:
    dbname = args[1]
    args = args[2:]
    #Checks to prevent overwriting
    if dbname + ".nal" in os.listdir(".") or dbname + ".pal" in os.listdir("."):
      print("Error: a database with this name already exists")
      sys.exit()
    if "genbank" in os.listdir("."):
      print("A folder named 'genbank' (which is used by MultiGeneBlast to store downloaded GenBank files) is already present in this folder.\n Overwrite? <y/n>")
      user_response = eval(input())
      if user_response.upper() == "N":
        sys.exit(1)
      else:
        shutil.rmtree("genbank")
    if "WGS" in args and "wgs" in os.listdir("."):
      print("A folder named 'wgs' (which is used by MultiGeneBlast to store downloaded WGS GenBank files) is already present in this folder.\n Overwrite? <y/n>")
      user_response = eval(input())
      if user_response.upper() == "N":
        sys.exit(1)
      else:
        shutil.rmtree("wgs")
  else:
    dbname = ""
    args = []
  
  get_genbank_db(dbname, args, dbtype="nucl")  ##Can be commented out for testing of local parsing
  
  #Parse normal SEQ entries, descriptions
  parse_nseq(dbname, find_seqfiles("genbank"))

  if "WGS" in args:
    #Parse WGS entries, descriptions
    parse_ngnp(dbname, find_ngnpfiles("wgs"))
  
  #Combine descriptions files
  combine_txt_files(dbname)
  
  #Combine FASTA files
  combine_fasta_files(dbname)
  
  #Write txt file of names from combined FASTA file
  fasta_names(dbname + "_all.fasta", dbname + "_all.txt")

  #Create Blast database and TAR archives
  try:
    if "\n" not in open(dbname + "_all.fasta","r").read(10000):
      print("Error making BLAST database; no suitable sequences found in input.")
      log("Error making BLAST database; no suitable sequences found in input.", exit=True)
  except:
    pass
  make_blast_db(dbname, dbname + "_all.fasta", dbtype="nucl")
  
  #Clean up
  clean_up(dbname, dbtype="nucl")
  
def maingui(args, dbname, frame, outbox):
  freeze_support()
  logfile = open("makedb.log","w")
  startpath = os.getcwd()
  logfile.close()
  freeze_support()
  makegbdbprocess = Process(target=gui_runmgbndb, args=[args, dbname])
  makegbdbprocess.start()
  loglinesprocessed = []
  while True:
      processrunning = "n"
      if makegbdbprocess.is_alive():
          processrunning = "y"
      if processrunning == "y":
          sleep(1)
          loglines = [line for line in open(os.getcwd() + os.sep + "makedb.log","r").readlines() if line not in loglinesprocessed]
          for line in loglines:
            loglinesprocessed.append(line)
            outbox.text_insert(line)
          frame.update()
      elif processrunning == "y":
          pass
      else:
          break
  loglines = [line for line in open(os.getcwd() + os.sep + "makedb.log","r").readlines() if line not in loglinesprocessed]
  for line in loglines:
    loglinesprocessed.append(line)
    outbox.text_insert(line)
  frame.update()

def gui_runmgbndb(args, dbname):
  try:
    get_genbank_db(dbname, args, GUI="y", dbtype="nucl")  ##Can be commented out for testing of local parsing

    #Parse normal SEQ entries, descriptions
    parse_nseq(dbname, find_seqfiles("genbank"))

    if "WGS" in args:
      #Parse WGS entries, descriptions
      parse_ngnp(dbname, find_ngnpfiles("wgs"))

    #Combine descriptions files
    combine_txt_files(dbname)

    #Combine FASTA files
    combine_fasta_files(dbname)

    #Write txt file of names from combined FASTA file
    fasta_names(dbname + "_all.fasta", dbname + "_all.txt")

    #Create Blast database and TAR archives
    if "\n" not in open(dbname + "_all.fasta","r").read(10000):
      print("Error making BLAST database; no suitable sequences found in input.")
      log("Error making BLAST database; no suitable sequences found in input.", exit=True)
    make_blast_db(dbname, dbname + "_all.fasta", dbtype="nucl")

    #Clean up
    clean_up(dbname, dbtype="nucl")
  except:
    log("Error: MakeGBDB process terminated prematurely.\n", exit=True)

if __name__ == "__main__":
  main()
