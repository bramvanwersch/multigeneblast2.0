#!/usr/bin/env python
## Copyright (c) 2012 Marnix H. Medema
## University of Groningen
## Department of Microbial Physiology / Groningen Bioinformatics Centre
## License: GNU General Public License v3 or later
## A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

from .utils import *
from .gbkparser import *
from .download_genbank import extract_file

def parse_seq(dbname, seqfiles):
  #For all files
  outputfile1 = open(dbname + "_seq.fasta","w")
  outputfile1.close()
  outputfile2 = open(dbname + "_descrs.txt","w")
  outputfile2.close()
  for seqfile in seqfiles:
    #Parse seq file sequences
    extract_file(seqfile)
    seqfile = seqfile.replace(".gz","")
    outputfile1 = open(dbname + "_seq.fasta","a")
    print("Parsing " + seqfile)
    log("Parsing " + seqfile)
    allnames = []
    allseqs = []
    #Open seq file  
    infile = open(seqfile,"r")
    #Divide in GBK files
    tempgbks = []
    first = "y"
    gbk = ""
    #Process all entries but last
    for line in infile:
      if "LOCUS       " in line[:12] and "LOCUS" in line[:5] and first == "n":
        if gbk.count("     CDS           ") > 2:
          entry = "LOCUS       " + gbk
          proteins = gbk2proteins(entry)
          if proteins != "nocds":
            genomic_accnr = proteins[1]
            dnaseqlength = proteins[2]
            proteins = proteins[0]
            writefasta(proteins[0],proteins[1],outputfile1)
        gbk = line
      elif "LOCUS       " in line[:12] and "LOCUS" in line[:5] and first == "y":
        first = "n"
        gbk = line
      else:
        gbk += line
    #Process last entry
    entry = "LOCUS       " + gbk
    proteins = gbk2proteins(entry)
    if proteins != "nocds":
      genomic_accnr = proteins[1]
      dnaseqlength = proteins[2]
      proteins = proteins[0]
      writefasta(proteins[0],proteins[1],outputfile1)
    infile.close()
    outputfile1.close()
    #Parse seq file descriptions
    outputfile2 = open(dbname + "_descrs.txt","a")
    print("Parsing " + seqfile + " descriptions")
    log("Parsing " + seqfile + " descriptions")
    #Open seq file  
    infile = open(seqfile,"r")
    #Divide in GBK files
    tempgbks = []
    first = "y"
    gbk = ""
    for line in infile:
      if "LOCUS       " in line[:12] and "LOCUS" in line[:5] and first == "n":
        if gbk.count("     CDS           ") > 2:
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
      elif "LOCUS       " in line[:12] and "LOCUS" in line[:5] and first == "y":
        first = "n"
        gbk = line
      else:
        gbk += line
    if gbk.count("     CDS           ") > 2:
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
    outputfile2.close()
    os.remove(seqfile)

def parse_nseq(dbname, seqfiles):
  #For all files
  outputfile1 = open(dbname + "_seq.fasta","w")
  outputfile1.close()
  outputfile2 = open(dbname + "_descrs.txt","w")
  outputfile2.close()
  for seqfile in seqfiles:
    #Parse seq file sequences
    extract_file(seqfile)
    seqfile = seqfile.replace(".gz","")
    outputfile1 = open(dbname + "_seq.fasta","a")
    outputfile2 = open(dbname + "_descrs.txt","a")
    print("Parsing " + seqfile)
    log("Parsing " + seqfile)
    allnames = []
    allseqs = []
    #Open seq file  
    infile = open(seqfile,"r")
    #Divide in GBK files
    tempgbks = []
    first = "y"
    gbk = ""
    #Process all entries but last
    for line in infile:
      if "LOCUS       " in line and first == "n":
        if gbk.count("     CDS           ") > 2:
          entry = "LOCUS       " + gbk
          description = entry.partition("DEFINITION  ")[2].partition("ACCESSION   ")[0].replace("\n"," ").replace("            "," ").replace("  "," ").replace("  "," ")
          if "ACCESSION   " in entry:
            j = entry.split("ACCESSION   ")[1].split("\n")[0]
            accession = j.split(" ")[0]
            if len(accession) < 4:
              accession = ""
          if testaccession(accession) == "n":
            if "LOCUS       " in entry:
              j = entry.split("LOCUS       ")[1].split("\n")[0]
              accession = j.split(" ")[0]
              if len(accession) < 4:
                accession = ""
          sequence = cleandnaseq(entry.partition("\nORIGIN")[2].partition("\n")[2])
          writefasta([accession],[sequence],outputfile1)
          outputfile2.write(accession + "\t" + description + "\n")
        gbk = line
      elif "LOCUS       " in line and first == "y":
        first = "n"
        gbk = line
      else:
        gbk += line
    #Process last entry
    entry = "LOCUS       " + gbk
    description = entry.partition("DEFINITION  ")[2].partition("ACCESSION   ")[0].replace("\n"," ").replace("            "," ").replace("  "," ").replace("  "," ")
    accession = entry.partition("ACCESSION   ")[2].partition("\n")[0].partition(" ")[0]
    sequence = cleandnaseq(entry.partition("\nORIGIN")[2].partition("\n")[2])
    writefasta([accession],[sequence],outputfile1)
    outputfile2.write(accession + "\t" + description + "\n")
    infile.close()
    outputfile1.close()
    outputfile2.close()
    os.remove(seqfile)