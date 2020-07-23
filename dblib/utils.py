#!/usr/bin/env python
## Copyright (c) 2012 Marnix H. Medema
## University of Groningen
## Department of Microbial Physiology / Groningen Bioinformatics Centre
## License: GNU General Public License v3 or later
## A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

import os
import sys
import shutil
import fileinput
from string import ascii_letters
from multiprocessing import Process, freeze_support
import subprocess

def find_files(namepart):
  fnames = [fname for fname in os.listdir(".") if namepart in fname]
  return fnames
  
def find_seqfiles(path = "."):
  currentdir = os.getcwd()
  os.chdir(path)
  if path == ".":
    seqfiles = [fname for fname in find_files(".seq") if "con" not in fname]
  else:
    seqfiles = [path + "/" + fname for fname in find_files(".seq") if "con" not in fname]
  os.chdir(currentdir)
  return seqfiles
  
def find_confiles(path = "."):
  currentdir = os.getcwd()
  os.chdir(path)
  if path == ".":
    confiles = [fname for fname in find_files(".seq") if "con" in fname]
  else:
    confiles = [path + "/" + fname for fname in find_files(".seq") if "con" in fname]
  os.chdir(currentdir)
  return confiles

def find_ngnpfiles(path = "."):
  currentdir = os.getcwd()
  os.chdir(path)
  if path == ".":
    gnpfiles = [fname for fname in find_files(".fsa_nt")]
  else:
    gnpfiles = [path + "/" + fname for fname in find_files(".fsa_nt")]
  os.chdir(currentdir)
  return gnpfiles

def find_gnpfiles(path = "."):
  currentdir = os.getcwd()
  os.chdir(path)
  if path == ".":
    gnpfiles = [fname for fname in find_files(".gnp")]
  else:
    gnpfiles = [path + "/" + fname for fname in find_files(".gnp")]
  os.chdir(currentdir)
  return gnpfiles

def find_gbff_files(path = "."):
  currentdir = os.getcwd()
  os.chdir(path)
  if path == ".":
    gbff_files = [fname for fname in find_files(".gbff")]
  else:
    gbff_files = [path + "/" + fname for fname in find_files(".gbff")]
  os.chdir(currentdir)
  return gbff_files

def combine_txt_files(dbname):
  outfile = open(dbname + "_all.txt","w")
  txtfiles = [fname for fname in find_files("txt") if dbname in fname and fname != dbname + "_all.txt"]
  for txtfile in txtfiles:
    for line in fileinput.input(txtfile):
      outfile.write(line)
  outfile.close()
  outfile = open(dbname + "_all_descrs.txt","w")
  txtfiles = [fname for fname in find_files("txt") if "_descrs" in fname and dbname in fname and fname != dbname + "_all_descrs.txt"]
  for txtfile in txtfiles:
    for line in fileinput.input(txtfile):
      outfile.write(line)
  outfile.close()

def combine_fasta_files(dbname):
  outfile = open(dbname + "_all.fasta","w")
  fastafiles = [fname for fname in find_files("fasta") if dbname + "_" in fname]
  for fastafile in fastafiles:
    for line in fileinput.input(fastafile):
      outfile.write(line)
  outfile.close()

def fasta_names(fastafile, outfile):
  outfile = open(outfile,"w")
  infolines = []
  for line in fileinput.input(fastafile):
    if ">" in line:
      outfile.write(line)
  outfile.close()

def makeblastdb(dbname, infile, dbtype):
  new_env = os.environ.copy()
  makeblastdbcommand = "makeblastdb -dbtype " + dbtype + " -out " + dbname + " -in " + infile
  output = subprocess.Popen(makeblastdbcommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=new_env)
  output = output.stdout.read()
  if "Error" in output:
    if "it is empty" in output:
      print("Error making BLAST database; no suitable sequences found in input.")
      log("Error making BLAST database; no suitable sequences found in input.", exit=True)
    else:
      print("Error making BLAST database.")
      log("Error making BLAST database.", exit=True)

def make_blast_db(dbname, infile, frame=None, GUI="n", dbtype="prot"):
  #print "Constructing Blast+ database"
  makeblastdbprocess = Process(target=makeblastdb, args=[dbname, infile, dbtype])
  makeblastdbprocess.start()
  while True:
      processrunning = "n"
      if makeblastdbprocess.is_alive():
          processrunning = "y"
      if processrunning == "y" and GUI == "y":
          frame.update()
      elif processrunning == "y":
          pass
      else:
          break

def testaccession(accession):
  #Test if accession number is probably real GenBank/RefSeq acc nr
  numbers = list(range(0,10))
  letters = []
  for i in ascii_letters:
    letters.append(i)
  nrnumbers = 0
  nrletters = 0
  for i in accession:
    if i in letters:
      nrletters += 1
    try:
      j = int(i)
      if j in numbers:
        nrnumbers += 1
    except:
      pass
  test = "y"
  if nrnumbers < 3 or nrletters < 1:
    test = "n"
  return test

def clean_up(dbname, frame=None, outbox=None, GUI="n", exit=False, dbtype="prot"):
  tags = ["_con.fasta", "_seq.fasta", "_wgs.fasta", "_all.fasta", "_all.txt", "_descrs.txt", "_wgs_descrs.txt", "_con_descrs.txt", "_dbbuild.fasta"]
  if exit:
    tags += [".cords.tar", ".phr", ".pal", ".pin", ".psq", "_all.txt", "_all_descrs.txt"]
  for tag in tags:
    if dbname + tag in os.listdir("."):
      try:
        os.remove(dbname + tag)
      except:
        pass
    if GUI == "y":
      frame.update()
  if dbtype == "prot":
    ext = ".pal"
  else:
    ext = ".nal"
  if dbname + ext not in os.listdir(".") and not exit:
    try:
      palfile = open(dbname + ext,"w")
      palfile.write("TITLE " + dbname + "\nDBLIST " + dbname + "\n")
      palfile.close()
    except:
      if GUI == "y":
        outbox.text_insert("Failure in generating NAL/PAL file for database\n")
  if GUI == "y":
    frame.update()
  if "proteininfo" in os.listdir("."):
    try:
      shutil.rmtree("proteininfo")
    except:
      pass
  if "genbank" in os.listdir("."):
    try:
      shutil.rmtree("genbank")
    except:
      pass
  if GUI == "y":
    frame.update()
  if "wgs" in os.listdir("."):
    try:
      shutil.rmtree("wgs")
    except:
      pass
  if GUI == "y":
    frame.update()

def log(message, exit=False, retcode=1):
    logfile = open('makedb.log', 'a', 1)
    logfile.write(message + '\n')
    logfile.close()
    if exit:
      sys.exit(retcode)

def fix_accession(accession):
    illegal_chars  = '''!"#$%&()*+,:;=>?@[]^`'{|} '''
    for char in accession:
      if char in illegal_chars:
        accession = accession.replace(char, "_")
    accession = accession.replace("\n","")
    return accession

def get_accession(entry, entrynr, entrytype, accessions, seqfilename, dbname):
    seqfilename = seqfilename.rpartition(os.sep)[2]
    entry = entry.replace("\r","\n")
    if entrytype == "gbk":
      if "ACCESSION   " in entry:
        accession = fix_accession(entry.partition("ACCESSION   ")[2].partition("\n")[0].partition(" ")[0])
        if len(accession) < 2:
          if "LOCUS       " in entry:
            accession = fix_accession(entry.partition("LOCUS       ")[2].partition("\n")[0].partition(" ")[0])
            if len(accession) < 2:
              accession = fix_accession(seqfilename + "_" + str(entrynr))
        if accession in accessions:
          accession = fix_accession(seqfilename + "_" + str(entrynr))
      elif "LOCUS       " in entry:
        accession = fix_accession(entry.partition("LOCUS       ")[2].partition("\n")[0].partition(" ")[0])
        if len(accession) < 2 or accession in accessions:
          accession = fix_accession(seqfilename + "_" + str(entrynr))
      else:
        accession = fix_accession(seqfilename + "_" + str(entrynr))
      if accession in accessions:
        print("Error: multiple input files with same accession number:", accession)
        clean_up(dbname)
        sys.exit()
      else:
        accessions.append(accession)
      return accession, accessions
    elif entrytype == "embl":
      if "AC   " in entry:
        accession = fix_accession(entry.partition("AC   ")[2].partition("\n")[0].replace(";","").partition(" ")[0])
        if len(accession) < 2:
          if "ID   " in entry:
            accession = fix_accession(entry.partition("ID   ")[2].partition("\n")[0].partition(";")[0])
            if len(accession) < 2:
              accession = fix_accession(seqfilename + "_" + str(entrynr))
        if accession in accessions:
          accession = fix_accession(seqfilename + "_" + str(entrynr))
      elif "ID   " in entry:
        accession = fix_accession(entry.partition("ID   ")[2].partition("\n")[0].partition(";")[0])
        if len(accession) < 2 or accession in accessions:
          accession = fix_accession(seqfilename + "_" + str(entrynr))
      else:
        accession = fix_accession(seqfilename + "_" + str(entrynr))
      if accession in accessions:
        print("Error: multiple input files with same accession number:", accession)
        clean_up(dbname)
        sys.exit()
      else:
        accessions.append(accession)
      return accession, accessions
