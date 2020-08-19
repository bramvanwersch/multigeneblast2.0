#!/usr/bin/env python
## Copyright (c) 2012 Marnix H. Medema
## University of Groningen
## Department of Microbial Physiology / Groningen Bioinformatics Centre
## License: GNU General Public License v3 or later
## A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

from string import ascii_letters
import string

def reverse_complement(sequence, frame=None, outbox=None, GUI="n"):
    complement = string.maketrans('atcgn', 'tagcn')
    return sequence.lower().translate(complement)[::-1]

def parsegenes(genes):
  genedict = {}
  genelist = []
  joinlist = []
  joindict = {}
  accessiondict = {}
  locustagdict = {}
  genenr = 0
  for i in genes:
    i = i.split("     gene            ")[0]
    join = "no"
    genenr += 1
    #Find gene location info for each gene
    if "complement" in i.split("\n")[0].lower() and i.split("\n")[0][-1] == ")":
      location = i.split("\n")[0]
    elif "complement" in i.split("\n")[0].lower() and i.split("\n")[0][-1] != ")":
      location = i.split("   /")[0]
      while ")" not in location.replace(" ","")[-3:]:
        locationlist = location.split("\n")
        locationlist = locationlist[:-1]
        location = ""
        for i in locationlist:
          location = location + "i"
      location = location.replace("\n","")
      location = location.replace(" ","")
    elif "join" in i.split("\n")[0].lower() and i.split("\n")[0][-1] == ")":
      location = i.split("\n")[0]
    elif "join" in i.split("\n")[0].lower() and i.split("\n")[0][-1] != ")":
      location = i.split("/")[0]
      while ")" not in location.replace(" ","")[-3:]:
        locationlist = location.split("\n")
        locationlist = locationlist[:-1]
        location = ""
        for i in locationlist:
          location = location + "i"
      location = location.replace("\n","")
      location = location.replace(" ","")
    else:
      location = i.split("\n")[0]
    #location info found in embl file, now extract start and end positions
    if "complement" in location.lower():
      location = location.lower()
      location = location.split("complement(")[1][:-1]
      if "join(" in location.lower():
        join = "yes"
        location = location.lower()
        location2 = location.split("join(")[1][:-1]
        start = location2.split(",")[0]
        start = start.split("..")[0]
        start = start.replace("<","").replace(")","").replace("(","")
        end = location2.split(",")[-1]
        if ".." in end:
          end = end.split("..")[1]
        end = end.replace(">","").replace(")","").replace("(","")
        joinedparts = location2.split(",")
        joinedparts2 = []
        for j in joinedparts:
          newjoinedpart = j.replace("<","")
          newjoinedpart = newjoinedpart.replace(">","").replace("(","")
          joinedparts2.append(newjoinedpart)
      else:
        start = location.split("..")[0]
        start = start.replace("<","").replace(")","").replace("(","")
        end = location.split("..")[1]
        end = end.replace(">","").replace(")","").replace("(","")
      strand = "-"
    else:
      if "join(" in location.lower():
        join = "yes"
        location = location.lower()
        location2 = location.split("join(")[1][:-1]
        start = location2.split(",")[0]
        start = start.split("..")[0]
        start = start.replace("<","").replace(")","").replace("(","")
        end = location2.split(",")[-1]
        if ".." in end:
          end = end.split("..")[1]
        end = end.replace(">","").replace(")","").replace("(","")
        joinedparts = location2.split(",")
        joinedparts2 = []
        for j in joinedparts:
          newjoinedpart = j.replace("<","")
          newjoinedpart = newjoinedpart.replace(">","")
          joinedparts2.append(newjoinedpart)
      else:
        start = location.split("..")[0]
        start = start.replace("<","").replace(")","").replace("(","")
        end = location.split("..")[1]
        end = end.replace(">","").replace(")","").replace("(","")
      strand = "+"
    if ":" in end:
      end = end.split(":")[1]
    if ":" in start:
      start = start.split(":")[1]
    if "," in end:
      end = end.split(",")[1]
    if "." in end:
      end = end.split(".")[1]
    if "," in start:
      start = start.split(",")[0]
    if "." in start:
      start = start.split(".")[0]
    if not start.isdigit() or not end.isdigit():
      continue
    if int(start) > int(end):
      start2 = end
      end2 = start
      start = start2
      end = end2
    #Correct for alternative codon start positions
    if "codon_start=" in i.lower():
      codonstart = i.lower().split("codon_start=")[1][0]
      if strand == "+":
        start = str(int(start) +  (int(codonstart) - 1))
      elif strand == "-":
        end = str(int(end) - (int(codonstart) - 1))
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
    #Find gene name or locus tag
    a = 0
    b = 0
    while b == 0:
      line = i.split("\n")[a]
      locustag = ""
      if "locus_tag=" in line:
        locustag = (line.split("locus_tag=")[1][1:-1]).replace(" ","_")
        locustag = locustag.replace("\\","_")
        locustag = locustag.replace("/","_")
        b += 1
      elif "locus_tag=" in line.lower():
        locustag = (line.lower().split("locus_tag=")[1][1:-1]).replace(" ","_")
        locustag = locustag.replace("\\","_")
        locustag = locustag.replace("/","_")
        b += 1
      elif a == (nrlines - 1):
        if locustag == "":
          locustag = "none"
        b += 1
      else:
        a += 1
    a = 0
    b = 0
    while b == 0:
      line = i.split("\n")[a]
      if "gene=" in line:
        locustag = (line.split("gene=")[1][1:-1]).replace(" ","_")
        locustag = locustag.replace("\\","_")
        locustag = locustag.replace("/","_")
        b += 1
      elif "gene=" in line.lower():
        locustag = (line.lower().split("gene=")[1][1:-1]).replace(" ","_")
        locustag = locustag.replace("\\","_")
        locustag = locustag.replace("/","_")
        b += 1
      elif a == (nrlines - 1):
        if locustag == "":
          locustag = "none"
        b += 1
      else:
        a += 1
    if locustag != "none":
      locustagdict[accnr.rpartition(".")[0]] = locustag
    if accnr == "no_accession_number_found" and locustag != "none":
      accnr = locustag
    if genename == "" and locustag != "none":
      genename = locustag
    #print genename
    #Find sequence for each gene
    a = 0                                             ###Not all gbks contain protein sequences as translations, therefore sequences from gene clusters are now extracted from the __database_file_label at a later stage if sequence is not in gbk
    b = 0
    sequence = ""
    while b < 2:
      line = i.split("\n")[a]
      if "translation=" in line:
        sequence = line.split("translation=")[1][1:]
        b += 1
        a += 1
        if line.count('"') > 1:
          sequence = line.split("translation=")[1][1:-1]
          b = 2
      elif "translation=" in line.lower():
        sequence = line.lower().split("translation=")[1][1:]
        b += 1
        a += 1
        if line.count('"') > 1:
          sequence = line.lower().split("translation=")[1][1:-1]
          b = 2
      elif a == (nrlines - 1):
        sequence = ""
        b = 2
      elif b == 1:
        if '"' in line:
          seqline = line.replace(" ","")
          seqline = seqline.split('"')[0]
          sequence = sequence + seqline
          b += 1
        else:
          seqline = line.replace(" ","")
          sequence = sequence + seqline
        a += 1
      else:
        a += 1
    sequence = sequence.upper()
    #Quality-check sequence
    forbiddencharacters = ["'",'"','=',';',':','[',']','>','<','|','\\',"/",'*','-','_','.',',','?',')','(','^','#','!','`','~','+','{','}','@','$','%','&']
    for z in forbiddencharacters:
      if z in sequence:
        sequence = ""
    #Find annotation for each gene
    a = 0
    b = 0
    while b == 0:
      line = i.split("\n")[a]
      if "product=" in line:
        annotation = line.split("product=")[1][1:]
        annotation = annotation.replace(" ","_")
        if len(annotation) > 0 and annotation[-1] == '"':
          annotation = annotation[:-1]
        else:
          nextline = i.split("\n")[a + 1]
          if "                     " in nextline and "/" not in nextline and '"' in nextline:
            annotation = annotation + " " + nextline
            while "  " in annotation:
              annotation = annotation.replace("  ", " ")
        if len(annotation) > 0 and annotation[-1] == '"':
          annotation = annotation[:-1]
        b += 1
      elif "product=" in line.lower():
        annotation = line.lower().split("product=")[1][1:]
        annotation = annotation.replace(" ","_")
        if len(annotation) > 0 and annotation[-1] == '"':
          annotation = annotation[:-1]
        else:
          nextline = i.split("\n")[a + 1]
          if "                     " in nextline and "/" not in nextline and '"' in nextline:
            annotation = annotation + " " + nextline
            while "  " in annotation:
              annotation = annotation.replace("  ", " ")
        if len(annotation) > 0 and annotation[-1] == '"':
          annotation = annotation[:-1]
        b += 1
      elif a == (nrlines - 1):
        annotation = "not_annotated"
        b += 1
      else:
        a += 1
    accessiondict[genename] = accnr
    if join == "yes":
      joinlist.append(genename)
      joindict[genename] = joinedparts2
    #Save data to dictionary
    if len(genename) > 1:
      genedict[genename] = [start,end,strand,annotation,sequence]
      genelist.append(genename)
  return [genelist, genedict, joinlist, joindict, accessiondict, locustagdict]

def parsegenesgnp(genes,genpept,protaccession,sequence,annotation):
  genedict = {}
  genelist = []
  joinlist = []
  joindict = {}
  accessiondict = {}
  genenr = 0
  for i in genes:
    i = i.split("     gene            ")[0]
    join = "no"
    genenr += 1
    #Find gene location info for each gene
    if "complement" in i.split("\n")[0].lower() and i.split("\n")[0][-1] == ")":
      location = i.split("\n")[0]
    elif "complement" in i.split("\n")[0].lower() and i.split("\n")[0][-1] != ")":
      location = i.split("   /")[0]
      while ")" not in location.replace(" ","")[-3:]:
        locationlist = location.split("\n")
        locationlist = locationlist[:-1]
        location = ""
        for i in locationlist:
          location = location + "i"
      location = location.replace("\n","")
      location = location.replace(" ","")
    elif "join" in i.split("\n")[0].lower() and i.split("\n")[0][-1] == ")":
      location = i.split("\n")[0]
    elif "join" in i.split("\n")[0].lower() and i.split("\n")[0][-1] != ")":
      location = i.split("/")[0]
      while ")" not in location.replace(" ","")[-3:]:
        locationlist = location.split("\n")
        locationlist = locationlist[:-1]
        location = ""
        for i in locationlist:
          location = location + "i"
      location = location.replace("\n","")
      location = location.replace(" ","")
    else:
      location = i.split("\n")[0]
    #location info found in embl file, now extract start and end positions
    if "complement" in location.lower():
      location = location.lower()
      location = location.split("complement(")[1][:-1]
      if "join(" in location.lower():
        join = "yes"
        location = location.lower()
        location2 = location.split("join(")[1][:-1]
        start = location2.split(",")[0]
        start = start.split("..")[0]
        start = start.replace("<","")
        end = location2.split(",")[-1]
        if ".." in end:
          end = end.split("..")[1]
        end = end.replace(">","")
        joinedparts = location2.split(",")
        joinedparts2 = []
        for j in joinedparts:
          newjoinedpart = j.replace("<","")
          newjoinedpart = newjoinedpart.replace(">","")
          joinedparts2.append(newjoinedpart)
      else:
        start = location.split("..")[0]
        start = start.replace("<","")
        end = location.split("..")[1]
        end = end.replace(">","")
      strand = "-"
    else:
      if "join(" in location.lower():
        join = "yes"
        location = location.lower()
        location2 = location.split("join(")[1][:-1]
        start = location2.split(",")[0]
        start = start.split("..")[0]
        start = start.replace("<","")
        end = location2.split(",")[-1]
        if ".." in end:
          end = end.split("..")[1]
        end = end.replace(">","")
        joinedparts = location2.split(",")
        joinedparts2 = []
        for j in joinedparts:
          newjoinedpart = j.replace("<","")
          newjoinedpart = newjoinedpart.replace(">","")
          joinedparts2.append(newjoinedpart)
      else:
        start = location.split("..")[0]
        start = start.replace("<","")
        if ".." in location:
          end = location.split("..")[1]
        else:
          end = location
        end = end.replace(">","")
      strand = "+"
    if int(start) > int(end):
      start2 = end
      end2 = start
      start = start2
      end = end2
    #Correct for alternative codon start positions
    if "codon_start=" in i.lower():
      codonstart = i.lower().split("codon_start=")[1][0]
      if strand == "+":
        start = str(int(start) +  (int(codonstart) - 1))
      elif strand == "-":
        end = str(int(end) - (int(codonstart) - 1))
    #Quality-check sequence
    forbiddencharacters = ["'",'"','=',';',':','[',']','>','<','|','\\',"/",'*','-','_','.',',','?',')','(','^','#','!','`','~','+','{','}','@','$','%','&']
    for z in forbiddencharacters:
      if z in sequence:
        sequence = ""
    #Find correct coding info
    codinginfo = genpept.split('coded_by="')[1].split('"')[0]
    if "complement(" in codinginfo:
      if "join(" in codinginfo:
        codinginfo = codinginfo.split("complement(")[1][:-1].replace("\n","").replace(" ","").replace("\t","")
        codinginfo = codinginfo.split("join(")[1][:-1].strip()
        strand = "-"
        nucacc = codinginfo.split(":")[0].split(".")[0]
        start = codinginfo.split(":")[1].split("..")[0]
        end = codinginfo.split(":")[-1].split("..")[-1]
      else:
        codinginfo = codinginfo.split("complement(")[1][:-1].replace("\n","").replace(" ","").replace("\t","")
        strand = "-"
        nucacc = codinginfo.split(":")[0].split(".")[0]
        start = codinginfo.split(":")[1].split("..")[0]
        end = codinginfo.split(":")[-1].split("..")[-1]
    else:
      if "join(" in codinginfo:
        codinginfo = codinginfo.split("join(")[1][:-1].replace("\n","").replace(" ","").replace("\t","")
        strand = "+"
        nucacc = codinginfo.split(":")[0].split(".")[0]
        start = codinginfo.split(":")[1].split("..")[0]
        end = codinginfo.split(":")[-1].split("..")[-1]
      else:
        codinginfo = codinginfo.replace("\n","").replace(" ","").replace("\t","")
        strand = "+"
        nucacc = codinginfo.split(":")[0].split(".")[0]
        start = codinginfo.split(":")[1].split("..")[0]
        end = codinginfo.split(":")[-1].split("..")[-1]
    start = start.replace(">","").replace("<","")
    end = end.replace(">","").replace("<","")
    #Save data to dictionary
    if len(protaccession) > 1:
      genedict[protaccession] = [start,end,strand,annotation,sequence.upper()]
      genelist.append(protaccession)
  return [genelist, genedict, joinlist, joindict, accessiondict, nucacc]

def extractprotfasta(genelist,genedict,dnaseq,rc_dnaseq,joinlist,joindict,accessiondict,nucleotidename, locustagdict):
  names = []
  seqs = []
  for i in genelist:
    genename = i
    if genename in locustagdict:
      locustag = locustagdict[genename]
    if genename.partition(".")[0] in locustagdict:
      locustag = locustagdict[genename.partition(".")[0]]
    else:
      locustag = "no_locus_tag"
    #If suitable translation found in gbk, use that
    if len(genedict[i][4]) > 5:
      protseq = genedict[i][4]
      i = genedict[i]
    #If no suitable translation found in gbk, extract from DNA sequence
    else:
      i = genedict[i]
      y = int(i[0])
      z = int(i[1])
      if i[2] == "+":
        if genename in joinlist:
          geneseq = ""
          for j in joindict[genename]:
            if ":" in j:
              j = j.split(":")[1]
            j = j.replace(")","").replace("(","")
            partstart = int(j.split("..")[0])
            if ".." in j:
              partend = int(j.split("..")[1])
            else:
              partend = int(j)
            geneseqpart = dnaseq[(partstart - 1):partend]
            geneseq = geneseq + geneseqpart
        else:
          geneseq = dnaseq[(y - 1):z]
        protseq = translate(geneseq)
      elif i[2] == "-":
        if genename in joinlist:
          geneseq = ""
          joinlistrev = joindict[genename]
          joinlistrev.reverse()
          for j in joinlistrev:
            if ":" in j:
              j = j.split(":")[1]
            j = j.replace(")","").replace("(","")
            partstart = int(j.split("..")[0])
            if ".." in j:
              partend = int(j.split("..")[1])
            else:
              partend = int(j)
            geneseqpart = rc_dnaseq[(len(rc_dnaseq) - partend):(len(rc_dnaseq) - partstart + 1)]
            geneseq = geneseq + geneseqpart
        else:
          geneseq = rc_dnaseq[(len(rc_dnaseq) - z):(len(rc_dnaseq) - y + 1)]
        protseq = translate(geneseq)
    name = nucleotidename + "|" + i[0] + "-" + i[1] + "|" + i[2] + "|" + genename.split(".")[0] + "|" + i[3] + "|" + locustag
    seqs.append(protseq)
    names.append(name)
  proteins = [names,seqs,genelist,genedict,accessiondict]
  return proteins

def extractprotfastacontig(genelist,genedict,dnaseq,rc_dnaseq,joinlist,joindict,accessiondict,nucleotidename, locustagdict):
  names = []
  seqs = []
  for i in genelist:
    genename = i
    if genename in locustagdict:
      locustag = locustagdict[genename]
    else:
      locustag = "no_locus_tag"
    #If suitable translation found in gbk, use that
    if len(genedict[i][4]) > 5:
      protseq = genedict[i][4]
      i = genedict[i]
      name = nucleotidename + "|" + i[0] + "-" + i[1] + "|" + i[2] + "|" + genename.split(".")[0] + "|" + i[3] + "|" + locustag
      seqs.append(protseq)
      names.append(name)
  proteins = [names,seqs,genelist,genedict,accessiondict]
  return proteins

def extractprotfastagnp(genelist,genedict,protseq,accessiondict,nucleotidename,protaccession, locustagdict):
  names = []
  seqs = []
  for i in genelist:
    genename = i
    #If suitable translation found in gbk, use that
    if len(genedict[i][4]) > 5:
      protseq = genedict[i][4]
    i = genedict[i]
    if genename in locustagdict:
      locustag = locustagdict[genename]
    else:
      locustag = "no_locus_tag"
    name = nucleotidename + "|" + i[0] + "-" + i[1] + "|" + i[2] + "|" + genename + "|" + i[3] + "|" + locustag
    seqs.append(protseq)
    names.append(name)
  proteins = [names,seqs,genelist,genedict,accessiondict]
  return proteins

def translate(sequence):
  #Translation table standard genetic code; according to http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
  transldict = { 'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',
                 'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C', 
                 'TTA': 'L', 'TCA': 'S', 'TAA': 'X', 'TGA': 'X', 
                 'TTG': 'L', 'TCG': 'S', 'TAG': 'X', 'TGG': 'W', 
                 'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R', 
                 'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R', 
                 'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R', 
                 'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R', 
                 'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S', 
                 'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S', 
                 'ATA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R', 
                 'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R', 
                 'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G', 
                 'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G', 
                 'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G', 
                 'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G',
                 'ttt': 'F', 'tct': 'S', 'tat': 'Y', 'tgt': 'C',
                 'ttc': 'F', 'tcc': 'S', 'tac': 'Y', 'tgc': 'C',
                 'tta': 'L', 'tca': 'S', 'taa': 'X', 'tga': 'X',
                 'ttg': 'L', 'tcg': 'S', 'tag': 'X', 'tgg': 'W',
                 'ctt': 'L', 'cct': 'P', 'cat': 'H', 'cgt': 'R',
                 'ctc': 'L', 'ccc': 'P', 'cac': 'H', 'cgc': 'R',
                 'cta': 'L', 'cca': 'P', 'caa': 'Q', 'cga': 'R',
                 'ctg': 'L', 'ccg': 'P', 'cag': 'Q', 'cgg': 'R',
                 'att': 'I', 'act': 'T', 'aat': 'N', 'agt': 'S',
                 'atc': 'I', 'acc': 'T', 'aac': 'N', 'agc': 'S',
                 'ata': 'I', 'aca': 'T', 'aaa': 'K', 'aga': 'R',
                 'atg': 'M', 'acg': 'T', 'aag': 'K', 'agg': 'R',
                 'gtt': 'V', 'gct': 'A', 'gat': 'D', 'ggt': 'G',
                 'gtc': 'V', 'gcc': 'A', 'gac': 'D', 'ggc': 'G',
                 'gta': 'V', 'gca': 'A', 'gaa': 'E', 'gga': 'G',
                 'gtg': 'V', 'gcg': 'A', 'gag': 'E', 'ggg': 'G'}
  triplets = []
  triplet = ""
  a = 0
  for i in sequence:
    if a < 2:
      a += 1
      triplet = triplet + i
    elif a == 2:
      triplet = triplet + i
      triplets.append(triplet)
      triplet = ""
      a = 0
  protseq = ""
  aanr = 0
  for i in triplets:
    aanr += 1
    if aanr == 1:
      protseq = protseq + "M"
    else:
      if "n" in i or "N" in i or i not in list(transldict.keys()):
        protseq = protseq + "X"
      else:
        protseq = protseq + transldict[i]
  if  len(protseq) > 0 and protseq[-1] == "*":
    protseq = protseq[:-1]
  return protseq

def writefasta(names,seqs,out_file):
  e = 0
  f = len(names) - 1
  try:
    while e <= f:
      out_file.write(">")
      out_file.write(names[e])
      out_file.write("\n")
      out_file.write(seqs[e])
      out_file.write("\n")
      e += 1
  except(IOError,OSError,NotImplementedError):
    print("FASTA file not created.", file=sys.stderr)

def cleandnaseq(dnaseq):
  dnaseq = dnaseq.replace(" ","")
  dnaseq = dnaseq.replace("\t","")
  dnaseq = dnaseq.replace("\n","")
  dnaseq = dnaseq.replace("0","")
  dnaseq = dnaseq.replace("1","")
  dnaseq = dnaseq.replace("2","")
  dnaseq = dnaseq.replace("3","")
  dnaseq = dnaseq.replace("4","")
  dnaseq = dnaseq.replace("5","")
  dnaseq = dnaseq.replace("6","")
  dnaseq = dnaseq.replace("7","")
  dnaseq = dnaseq.replace("8","")
  dnaseq = dnaseq.replace("9","")
  dnaseq = dnaseq.replace("/","")
  return dnaseq

def gbk2proteins(gbkfile):
  filetext = gbkfile
  filetext = filetext.replace("\r","\n")
  if "\n     CDS             " not in filetext or ("\nORIGIN" not in filetext and "\nCONTIG" not in filetext):
    return "nocds"
  if "\nORIGIN" in filetext:
    cdspart = filetext.split("\nORIGIN")[0]
    #Extract DNA sequence and calculate reverse complement of it
    dnaseq = filetext.split("\nORIGIN")[1]
    dnaseq = cleandnaseq(dnaseq)
    dnaseqlength = len(dnaseq)
    rc_dnaseq = reverse_complement(dnaseq)
  elif "\nCONTIG" in filetext:
    cdspart = filetext.split("\nCONTIG")[0]
    dnaseq = ""
    rc_dnaseq = ""
    dnaseqlength = int(cdspart.partition("LOCUS       ")[2].partition(" bp")[0].rpartition(" ")[2])
  #Extract genes
  genes = cdspart.split("\n     CDS             ")
  genes = genes[1:]
  genesdetails = parsegenes(genes)
  genelist = genesdetails[0]
  genedict = genesdetails[1]
  joinlist = genesdetails[2]
  joindict = genesdetails[3]
  accessiondict = genesdetails[4]
  locustagdict = genesdetails[5]
  #Locate all genes on DNA sequence and translate to protein sequence
  textlines = filetext.split("\n//")[0]
  textlines = textlines.split("\n")
  accession = ""
  for i in textlines:
    if accession == "":
      if "ACCESSION   " in i:
        j = i.split("ACCESSION   ")[1]
        accession = j.split(" ")[0]
        if len(accession) < 4:
          accession = ""
  #Test if accession number is probably real GenBank/RefSeq acc nr
  numbers = [nr for nr in range(0,10)]
  letters = [i for i in ascii_letters]
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
  if nrnumbers < 3 or nrletters < 1:
    accession = ""
  #If not, try finding accession number under 'LOCUS' tag
  if accession == "":
    for i in textlines:
      if accession == "":
        if "LOCUS       " in i:
          j = i.split("LOCUS       ")[1]
          accession = j.split(" ")[0]
          if len(accession) < 4:
            accession = ""
  #Test again
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
  if nrnumbers < 3 or nrletters < 1:
    accession = ""
  #Get lists of all predicted proteins from the GBK entry
  if "\nORIGIN" in filetext:
    proteins = extractprotfasta(genelist,genedict,dnaseq,rc_dnaseq,joinlist,joindict,accessiondict,accession, locustagdict)
  elif "\nCONTIG" in filetext:
    proteins = extractprotfastacontig(genelist,genedict,dnaseq,rc_dnaseq,joinlist,joindict,accessiondict,accession, locustagdict)
  return [proteins,accession,dnaseqlength]

def cleanprotseq(protseq):
  protseq = protseq.replace(" ","")
  protseq = protseq.replace("\t","")
  protseq = protseq.replace("\n","")
  protseq = protseq.replace("0","")
  protseq = protseq.replace("1","")
  protseq = protseq.replace("2","")
  protseq = protseq.replace("3","")
  protseq = protseq.replace("4","")
  protseq = protseq.replace("5","")
  protseq = protseq.replace("6","")
  protseq = protseq.replace("7","")
  protseq = protseq.replace("8","")
  protseq = protseq.replace("9","")
  protseq = protseq.replace("/","")
  return protseq

def gnp2proteins(gbkfile):
  filetext = gbkfile
  filetext = filetext.replace("\r","\n")
  if "\n     CDS             " not in filetext or "\nORIGIN" not in filetext:
    return "nocds"
  cdspart = filetext.split("\nORIGIN")[0]
  #Extract DNA sequence and calculate reverse complement of it
  protseq = filetext.split("\nORIGIN")[1]
  protseq = cleanprotseq(protseq)
  sequence = protseq.upper()
  protseqlength = len(protseq)
  #Extract genes
  genes = cdspart.split("\n     CDS             ")
  genes = genes[1:]
  textlines = filetext.split("\n//")[0]
  textlines = textlines.split("\n")
  locustagdict = {}
  protaccession = ""
  for i in textlines:
    if protaccession == "":
      if "LOCUS       " in i:
        j = i.split("LOCUS       ")[1]
        protaccession = j.split(" ")[0]
        if len(protaccession) < 4:
          protaccession = ""
  if '/product="' in filetext:
    annotation = filetext.split('/product="')[1].split('"')[0].replace(" ","_").replace("\n","")
    while "__" in annotation:
      annotation = annotation.replace("__","_")
  else:
    annotation = "not annotated"
  if '/locus_tag="' in filetext:
    locustag = filetext.split('/locus_tag="')[1].split('"')[0].replace(" ","_").replace("\n","")
  else:
    locustag = "none"
  if '/gene="' in filetext:
    locustag = filetext.split('/gene="')[1].split('"')[0].replace(" ","_").replace("\n","")
  if locustag != "none":
    locustagdict[protaccession.rpartition(".")[0]] = locustag
  genesdetails = parsegenesgnp(genes,filetext,protaccession,protseq,annotation)
  genelist = genesdetails[0]
  genedict = genesdetails[1]
  joinlist = genesdetails[2]
  joindict = genesdetails[3]
  accessiondict = genesdetails[4]
  accession = genesdetails[5]
  #Locate all genes on DNA sequence and translate to protein sequence
  textlines = filetext.split("\n//")[0]
  textlines = textlines.split("\n")
  proteins = extractprotfastagnp(genelist,genedict,protseq,accessiondict,accession,protaccession, locustagdict)
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
  if nrnumbers < 3 or nrletters < 1:
    accession = ""
  return [proteins,accession,protseqlength]