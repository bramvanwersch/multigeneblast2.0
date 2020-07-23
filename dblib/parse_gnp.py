#!/usr/bin/env python
## Copyright (c) 2012 Marnix H. Medema
## University of Groningen
## Department of Microbial Physiology / Groningen Bioinformatics Centre
## License: GNU General Public License v3 or later
## A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

from .utils import *
from .gbkparser import *
from .download_genbank import extract_file
import os

def parse_gnp(dbname, gnpfiles):
  #For all files
  outputfile = open(dbname + "_wgs.fasta","w")
  for gnpfile in gnpfiles:
    extract_file(gnpfile)
    gnpfile = gnpfile.replace(".gz","")
    print("Parsing " + gnpfile)
    log("Parsing " + gnpfile)
    allnames = []
    allseqs = []
    #Open gnp file  
    infile = open(gnpfile,"r")
    #Divide in GBK files
    gbks = ["LOCUS       " + gbk for gbk in infile.read().split("LOCUS       ")[1:] if gbk.count("     CDS           ") > 0]
    infile.close()
    #For all resultant GBKs, parse sequence
    for infile in gbks:
      proteins = gnp2proteins(infile)
      if proteins != "nocds":
        genomic_accnr = proteins[1]
        protseqlength = proteins[2]
        proteins = proteins[0]
        allnames = allnames + proteins[0]
        allseqs = allseqs + proteins[1]
      del infile
    writefasta(allnames,allseqs,outputfile)
    os.remove(gnpfile)
    del gbks
  outputfile.close()

def parse_gnp_descr(dbname, gbff_files):
  #For all files
  outputfile = open(dbname + "_wgs_descrs.txt","w")
  descriptiondict = {}
  for gbff_file in gbff_files:
    extract_file(gbff_file)
    gbff_file = gbff_file.replace(".gz","")
    print("Parsing " + gbff_file + " description")
    log("Parsing " + gbff_file + " description")
    allnames = []
    allseqs = []
    #Open seq file  
    infile = open(gbff_file,"r")
    text = infile.read().split("\n//")[0]
    infile.close()
    accession = ""
    description = ""
    if "ACCESSION   " in text:
      j = text.split("ACCESSION   ")[1].split("\n")[0]
      accession = j.split(" ")[0]
      if len(accession) < 6:
        accession = ""
      else:
        accession = accession[:6]
    if testaccession(accession) == "n":
      if "LOCUS       " in text:
        j = text.split("LOCUS       ")[1].split("\n")[0]
        accession = j.split(" ")[0]
        if len(accession) < 6:
          accession = ""
        else:
          accession = accession[:6]
    if "DEFINITION  " in text:
      description = text.split("DEFINITION  ")[1].split("\nACCESSION   ")[0].split("\nVERSION     ")[0].split("\nKEYWORDS    ")[0].split("\nSOURCE      ")[0].replace("\n","")
      while "  " in description:
        description = description.replace("  "," ")
    if accession != "":
      descriptiondict[accession] = description
    writestring = ""
    os.remove(gbff_file)
  for j in list(descriptiondict.keys()):
    writestring += j + "\t" + descriptiondict[j] + "\n"
  outputfile.write(writestring)
  outputfile.close()

def parse_ngnp(dbname, gnpfiles):
  #For all files
  outputfile1 = open(dbname + "_wgs.fasta","w")
  outputfile2 = open(dbname + "_wgs_descrs.txt","w")
  for gnpfile in gnpfiles:
    extract_file(gnpfile)
    gnpfile = gnpfile.replace(".gz","")
    print("Parsing " + gnpfile)
    log("Parsing " + gnpfile)
    #Open gnp file  
    infile = open(gnpfile,"r")
    #Divide in FASTA entries
    tempfastas = infile.read().split(">")[1:]
    infile.close()
    #For all resultant FASTA entry, parse sequence & description
    for entry in tempfastas:
      sequence = entry.partition("\n")[2].replace("\n","")
      accession = entry.partition("\n")[0].rpartition("|")[0].rpartition("|")[2].rpartition(".")[0]
      description = entry.partition("\n")[0].rpartition("|")[2]
      while len(description) > 0 and description[0] == " ":
        description = description[1:]
      outputfile2.write(accession + "\t" + description + "\n")
      writefasta([accession], [sequence], outputfile1)
    os.remove(gnpfile)
  outputfile1.close()
  outputfile2.close()