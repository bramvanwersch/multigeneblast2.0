#!/usr/bin/env python
## Copyright (c) 2012 Marnix H. Medema
## University of Groningen
## Department of Microbial Physiology / Groningen Bioinformatics Centre
## License: GNU General Public License v3 or later
## A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

import os
from .utils import *
from .gbkparser import *
from .download_genbank import extract_file

def parse_confiles(dbname, confiles):
  for confile in confiles:
    parse_con(dbname, confile)

def parse_con(dbname, confiles):
  #For all files
  confastafiles = []
  outputfile2 = open(dbname + "_con_descrs.txt","w")
  for confile in confiles:
    #Parsing CON sequences
    extract_file(confile)
    confile = confile.replace(".gz","")
    outputfile1 = open(confile.rpartition(".")[0] + ".fasta","w")
    confastafiles.append(confile.rpartition(".")[0] + ".fasta")
    print("Parsing " + confile)
    allnames = []
    allseqs = []
    #Open seq file  
    infile = open(confile,"r")
    #Divide in GBK files
    tempgbks = []
    first = "y"
    gbk = ""
    #Process all entries but last
    for line in infile:
      if "LOCUS       " in line and first == "n":
        if gbk.count("     CDS           ") > 2:
          entry = "LOCUS       " + gbk
          proteins = gbk2proteins(entry)
          if proteins != "nocds":
            genomic_accnr = proteins[1]
            proteins = proteins[0]
            writefasta(proteins[0],proteins[1],outputfile1)
        gbk = line
      elif "LOCUS       " in line and first == "y":
        first = "n"
        gbk = line
      else:
        gbk += line
    #Process last entry
    entry = "LOCUS       " + gbk
    proteins = gbk2proteins(entry)
    if proteins != "nocds":
      genomic_accnr = proteins[1]
      proteins = proteins[0]
      writefasta(proteins[0],proteins[1],outputfile1)
    infile.close()
    outputfile1.close()
    #Parsing CON descriptions
    print("Parsing " + confile + " descriptions")
    #Open seq file  
    infile = open(confile,"r")
    #Divide in GBK files
    tempgbks = []
    first = "y"
    gbk = ""
    for line in infile:
      if "LOCUS       " in line and first == "n":
        text = gbk.split("\n//")[0]
        accession = ""
        description = ""
        if "ACCESSION   " in text:
          j = text.split("ACCESSION   ")[1].split("\n")[0]
          accession = j.split(" ")[0]
          if len(accession) < 4:
            accession = ""
        if testaccession(accession) == "n":
          if "LOCUS       " in text:
            j = text.split("LOCUS       ")[1].split("\n")[0]
            accession = j.split(" ")[0]
            if len(accession) < 4:
              accession = ""
        if "DEFINITION  " in text:
          description = text.split("DEFINITION  ")[1].split("\nACCESSION   ")[0].split("\nVERSION     ")[0].split("\nKEYWORDS    ")[0].split("\nSOURCE      ")[0].replace("\n","")
          while "  " in description:
            description = description.replace("  "," ")
        if accession != "":
          outputfile2.write(accession + "\t" + description + "\n")
        gbk = line
      elif "LOCUS       " in line and first == "y":
        first = "n"
        gbk = line
      else:
        gbk += line
    text = gbk.split("\n//")[0]
    accession = ""
    description = ""
    if "ACCESSION   " in text:
      j = text.split("ACCESSION   ")[1].split("\n")[0]
      accession = j.split(" ")[0]
      if len(accession) < 4:
        accession = ""
    if testaccession(accession) == "n":
      if "LOCUS       " in text:
        j = text.split("LOCUS       ")[1].split("\n")[0]
        accession = j.split(" ")[0]
        if len(accession) < 4:
          accession = ""
    if "DEFINITION  " in text:
      description = text.split("DEFINITION  ")[1].split("\nACCESSION   ")[0].split("\nVERSION     ")[0].split("\nKEYWORDS    ")[0].split("\nSOURCE      ")[0].replace("\n","")
      while "  " in description:
        description = description.replace("  "," ")
    if accession != "":
      outputfile2.write(accession + "\t" + description + "\n")
    infile.close()
    os.remove(confile)
  outputfile2.close()
  merge_confastafiles(dbname, confastafiles)

def merge_confastafiles(dbname, confastafiles):
  print("Merging CON files")
  outfile = open(dbname + "_con.fasta","w")
  for fname in confastafiles:
    content = open(fname,"r")
    outfile.write(content.read())
    content.close()
    os.remove(fname)
  outfile.close()
