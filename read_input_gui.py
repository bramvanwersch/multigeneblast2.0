#!/usr/bin/env python
## Copyright (c) 2010,2011 Marnix H. Medema
## Department of Microbial Physiology / Groningen Bioinformatics Centre
## University of Groningen
## License: GNU General Public License v3 or later
## A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

from tkinter.messagebox import showerror
from multigeneblast import cleandnaseq, parse_dna_from_embl

def parsegenes_gui(genes):
  genelist = []
  genenr = 0
  for i in genes:
    i = i.split("     gene            ")[0]
    genenr += 1
    #Find gene name for each gene, preferably locus_tag, than gene, than protein_ID
    a = 0
    b = 0
    genename = ""
    nrlines = len(i.split("\n"))
    while b == 0:
      line = i.split("\n")[a]
      if "protein_id=" in line:
        genename = (line.split("protein_id=")[1][1:-1]).replace(" ","_")
        genename = genename.replace("\\","_")
        genename = genename.replace("/","_")
        b += 1
      elif "protein_id=" in line.lower():
        genename = (line.lower().split("protein_id=")[1][1:-1]).replace(" ","_")
        genename = genename.replace("\\","_")
        genename = genename.replace("/","_")
        b += 1
      elif a == (nrlines - 1):
        genename = ""
        b += 1
      else:
        a += 1
    if len(genename) > 1:
      accnr = genename
    else:
      accnr = "no_accession_number_found"
    a = 0
    b = 0
    nrlines = len(i.split("\n"))
    while b == 0:
      line = i.split("\n")[a]
      genername = ""
      if "gene=" in line:
        genename = (line.split("gene=")[1][1:-1]).replace(" ","_")
        genename = genename.replace("\\","_")
        genename = genename.replace("/","_")
        genername = genename
        b += 1
      elif "gene=" in line.lower():
        genename = (line.lower().split("gene=")[1][1:-1]).replace(" ","_")
        genename = genename.replace("\\","_")
        genename = genename.replace("/","_")
        genername = genename
        b += 1
      elif a == (nrlines - 1):
        if genername == "":
          genername = "prot_ID_" + str(genenr)
        b += 1
      else:
        a += 1
    a = 0
    b = 0
    nrlines = len(i.split("\n"))
    while b == 0:
      line = i.split("\n")[a]
      if "locus_tag=" in line:
        genename = (line.split("locus_tag=")[1][1:-1]).replace(" ","_")
        genename = genename.replace("\\","_")
        genename = genename.replace("/","_")
        b += 1
      elif "locus_tag=" in line.lower():
        genename = (line.lower().split("locus_tag=")[1][1:-1]).replace(" ","_")
        genename = genename.replace("\\","_")
        genename = genename.replace("/","_")
        b += 1
      elif a == (nrlines - 1):
        if genename == "":
          genename = "prot_ID_" + str(genenr)
        b += 1
      else:
        a += 1
    genelist.append(genename)
  return genelist

def gbk2proteins_gui(gbkfile):
  try:
    file = open(gbkfile,"r")
  except:
    showerror("File Error","Specified GBK file not found.")
    return 0,0
  filetext = file.read()
  filetext = filetext.replace("\r","\n")
  if "     CDS             " not in filetext or "\nORIGIN" not in filetext:
    showerror("File Error","GBK file not properly formatted, no sequence found or no CDS annotation found.")
    return 0,0
  cdspart = filetext.split("\nORIGIN")[0]
  #Extract DNA sequence and calculate reverse complement of it
  dnaseq = filetext.split("\nORIGIN")[1]
  dnaseq = cleandnaseq(dnaseq)
  dnaseqlength = len(dnaseq)
  if dnaseqlength < 1:
      showerror("File Error","No sequence found in GBK/EMBL file. Please provide an annotated nucleotide GBK/EMBL file with a DNA sequence.")
      return 0,0
  #Extract genes
  genes = cdspart.split("     CDS             ")
  genes = genes[1:]
  genes = parsegenes_gui(genes)
  return genes, dnaseqlength

def embl2proteins_gui(emblfile):
  try:
    file = open(emblfile,"r")
  except:
    showerror("File Error","Specified EMBL file not found.")
    return 0,0
  filetext = file.read()
  filetext = filetext.replace("\r","\n")
  if "FT   CDS " not in filetext or ("\nSQ" not in filetext):
      showerror("File Error","EMBL file not properly formatted, no sequence found or no CDS annotation found.")
      return 0,0
  cdspart = filetext.split("\nSQ  ")[0]
  #Extract DNA sequence and calculate reverse complement of it
  dnaseq = parse_dna_from_embl(filetext.split("\nSQ  ")[1])
  dnaseq = cleandnaseq(dnaseq)
  sequence = dnaseq
  if (sequence.count('A') + sequence.count('a') + sequence.count('C') + sequence.count('c') + sequence.count('G') + sequence.count('g') + sequence.count('T') + sequence.count('t')) < (0.5 * len(sequence)):
      showerror("File Error","Protein GBK/EMBL file provided. Please provide nucleotide GBK/EMBL file.")
      return 0,0
  dnaseqlength = len(dnaseq)
  if dnaseqlength < 1:
      showerror("File Error","No sequence found in GBK/EMBL file. Please provide an annotated nucleotide GBK/EMBL file with a DNA sequence.")
      return 0,0
  #Extract genes
  genes = cdspart.split("FT   CDS             ")
  genes = genes[1:]
  genes = parsegenes_gui(genes)
  return genes, dnaseqlength

def read_input_file_gui(infile):
  if len(infile)== 0:
    return 0,0
  ext = infile.rpartition(".")[2]
  if ext.lower() in ["gbk","gb","genbank"]:
    genes, dnaseqlength = gbk2proteins_gui(infile)
  elif ext.lower() in ["embl","emb"]:
    genes, dnaseqlength = embl2proteins_gui(infile)
  return genes, dnaseqlength

def checkfasta(infile):
  content = open(infile,"r").read()
  entries = [">" + entry for entry in content.split(">")][1:]
  if len(entries) < 2:
    showerror("File Error","Please provide a FASTA file with multiple entries, representing the amino acid sequences of your query genes.")
    return 0
  for entry in entries:
    entryseq = entry.partition("\n")[2]
    if (entryseq.count('A') + entryseq.count('a') + entryseq.count('C') + entryseq.count('c') + entryseq.count('G') + entryseq.count('g') + entryseq.count('T') + entryseq.count('t')) > (0.5 * len(entryseq)):
      showerror("File Error","Nucleotide FASTA sequences provided. Please provide a multi-FASTA file with amino acid sequences.")
      return 0
  return 1