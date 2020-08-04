import os
import sys
from string import ascii_letters
import string
from .utils import clean_up, get_accession
import urllib.request, urllib.error, urllib.parse
from urllib.request import Request, urlopen
from urllib.error import URLError, HTTPError
import http.client
from http.client import BadStatusLine,HTTPException
import time

# def reverse_complement(sequence, frame=None, outbox=None, GUI="n"):
#     complement = string.maketrans('atcgn', 'tagcn')
#     return sequence.lower().translate(complement)[::-1]


def reverse_complement(sequence, frame=None, outbox=None, GUI="n"):
  """
  Create the complement strand of a sequence. Use replacing for speed

  :param seq: a DNA sequence as a string
  :return: the complement of that string
  """
  complement = sequence
  complement = complement.replace("a", "x").replace("A", "X")
  complement = complement.replace("t", "a").replace("T", "A")
  complement = complement.replace("x", "t").replace("X", "T")

  complement = complement.replace("c", "x").replace("C", "X")
  complement = complement.replace("g", "c").replace("G", "C")
  complement = complement.replace("x", "g").replace("X", "G")
  return sequence.lower().translate(complement)[::-1]

def parsegenes(genes, filetype="genbank_derived", frame=None, outbox=None, GUI="n", nr=0):
  genedict = {}
  genelist = []
  joinlist = []
  joindict = {}
  accessiondict = {}
  locustagdict = {}
  genenr = 0
  for i in genes:
    if GUI == "y":
      frame.update()
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
        if ".." in location:
          end = location.split("..")[1]
        end = end.replace(">","").replace(")","").replace("(","")
      strand = "+"
    if ":" in end:
      end = end.split(":")[1]
    if ":" in start:
      start = start.split(":")[1]
    if "," in end:
      end = end.split(",")[1]
    if "," in start:
      start = start.split(",")[0]
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
    #Find gene name for each gene, preferably protein_ID, than locus_tag, than gene 
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
        genename = genename.replace('"','')
        b += 1
      elif "protein_id=" in line.lower():
        genename = (line.lower().split("protein_id=")[1][1:-1]).replace(" ","_")
        genename = genename.replace("\\","_")
        genename = genename.replace("/","_")
        genename = genename.replace('"','')
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
        locustag = locustag.replace('"','')
        b += 1
      elif "locus_tag=" in line.lower():
        locustag = (line.lower().split("locus_tag=")[1][1:-1]).replace(" ","_")
        locustag = locustag.replace("\\","_")
        locustag = locustag.replace("/","_")
        locustag = locustag.replace('"','')
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
        locustag = locustag.replace('"','')
        b += 1
      elif "gene=" in line.lower():
        locustag = (line.lower().split("gene=")[1][1:-1]).replace(" ","_")
        locustag = locustag.replace("\\","_")
        locustag = locustag.replace("/","_")
        locustag = locustag.replace('"','')
        b += 1
      elif a == (nrlines - 1):
        if locustag == "":
          locustag = "none"
        b += 1
      else:
        a += 1
    if locustag != "none":
      locustagdict[accnr.rpartition(".")[0]] = locustag
    if accnr == "no_accession_number_found" and locustag != "none" and filetype == "user_inputfile":
      accnr = locustag
    if genename == "" and locustag != "none":
      genename = locustag
    #Find sequence for each gene
    a = 0                                             ###Not all gbks contain protein sequences as translations, therefore sequences from gene clusters are now extracted from the database at a later stage if sequence is not in gbk
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
        if len(annotation) != 0 and annotation[-1] == '"':
          annotation = annotation[:-1]
        else:
          nextline = i.split("\n")[a + 1]
          if "                     " in nextline and "/" not in nextline and '"' in nextline:
            annotation = annotation + " " + nextline
            while "  " in annotation:
              annotation = annotation.replace("  ", " ")
        if len(annotation) == 0:
          annotation = "not_annotated"
        elif annotation[-1] == '"':
          annotation = annotation[:-1]
        b += 1
      elif "product=" in line.lower():
        annotation = line.lower().split("product=")[1][1:]
        annotation = annotation.replace(" ","_")
        if len(annotation) != 0 and annotation[-1] == '"':
          annotation = annotation[:-1]
        else:
          nextline = i.split("\n")[a + 1]
          if "                     " in nextline and "/" not in nextline and '"' in nextline:
            annotation = annotation + " " + nextline
            while "  " in annotation:
              annotation = annotation.replace("  ", " ")
        if len(annotation) == 0:
          annotation = "not_annotated"
        elif annotation[-1] == '"':
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
    #Remove illegal chars
    illegal_chars  = '''!"#$%&()*+,:;=>?@[]^`'{|} '''
    genename = "".join([char for char in genename if char not in illegal_chars])
    if len(genename) < 2:
      genename = "orf" + str(nr) + "_" + str(genenr)
    #Save data to dictionary
    if len(genename) > 1:
      genedict[genename] = [start,end,strand,annotation,sequence]
      genelist.append(genename)
  return [genelist, genedict, joinlist, joindict, accessiondict, locustagdict]

def cleandnaseq(dnaseq, frame=None, outbox=None, GUI="n"):
  dnaseq = dnaseq.replace(" ","")
  dnaseq = dnaseq.replace("\t","")
  dnaseq = dnaseq.replace("\n","")
  dnaseq = dnaseq.replace("0","")
  if GUI == "y":
    frame.update()
  dnaseq = dnaseq.replace("1","")
  dnaseq = dnaseq.replace("2","")
  dnaseq = dnaseq.replace("3","")
  dnaseq = dnaseq.replace("4","")
  if GUI == "y":
    frame.update()
  dnaseq = dnaseq.replace("5","")
  dnaseq = dnaseq.replace("6","")
  dnaseq = dnaseq.replace("7","")
  if GUI == "y":
    frame.update()
  dnaseq = dnaseq.replace("8","")
  dnaseq = dnaseq.replace("9","")
  dnaseq = dnaseq.replace("/","")
  if GUI == "y":
    frame.update()
  return dnaseq

def extractprotfasta(genelist,genedict,dnaseq,rc_dnaseq,joinlist,joindict,accessiondict,nucleotidename, locustagdict, frame=None, outbox=None, GUI="n"):
  names = []
  seqs = []
  for i in genelist:
    if GUI == "y":
      frame.update()
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
    name = nucleotidename + "|" + i[0] + "-" + i[1] + "|" + i[2] + "|" + genename.split(".")[0] + "|" + i[3].replace("\t","") + "|" + locustag
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

def parse_dna_from_embl(embl_string, frame, outbox, GUI):
    "Parse DNA sequence from EMBL input"
    seq_array = []
    lines = embl_string.split('\n')
    for line in lines:
        if GUI == "y":
          frame.update()
        if line.lower().find('sequence') > -1:
            continue
        line = line.strip()
        line = line.rstrip('0123456789')
        line = line.rstrip('/')
        line = line.strip()
        seq_array.append(line)
    return "".join(seq_array)

def embl2proteins(emblfile, nr=0, filetype="genbank_derived", frame=None, outbox=None, GUI="n", accession=None):
  filetext = emblfile
  filetext = filetext.replace("\r","\n")
  if "FT   CDS " not in filetext or ("\nSQ" not in filetext):
    if accession == None:
      accession = emblfile.partition("AC   ")[2].partition("\n")[0].replace(";","").partition(" ")[0].replace(" ","")
    if GUI == "n":
      print("No CDS or DNA sequence found in " + accession + ". Skipping file.")
    else:
      outbox.text_insert("No CDS or DNA sequence found in " + accession + ". Skipping file.\n")
      frame.update()
    return "nocds"
  cdspart = filetext.split("\nSQ  ")[0]
  #Extract DNA sequence and calculate reverse complement of it
  dnaseq = parse_dna_from_embl(filetext.split("\nSQ  ")[1], frame, outbox, GUI)
  dnaseq = cleandnaseq(dnaseq)
  sequence = dnaseq
  if (sequence.count('N') + sequence.count('n') + sequence.count('A') + sequence.count('a') + sequence.count('C') + sequence.count('c') + sequence.count('G') + sequence.count('g') + sequence.count('T') + sequence.count('t')) < (0.5 * len(sequence)):
      if GUI == "n":
        print("Protein GBK/EMBL file provided. Please provide nucleotide GBK/EMBL file.")
        sys.exit(1)
      else:
        outbox.text_insert("Protein GBK/EMBL file provided. Please provide nucleotide GBK/EMBL file.\n")
        frame.update()
        sys.exit(1)
  dnaseqlength = len(dnaseq)
  if GUI == "y":
    frame.update()
  rc_dnaseq = reverse_complement(dnaseq)
  if dnaseqlength < 1:
      if GUI == "n":
        print("No sequence found in GBK/EMBL file. Please provide an annotated nucleotide GBK/EMBL file with a DNA sequence.")
        sys.exit(1)
      else:
        outbox.text_insert("No sequence found in GBK/EMBL file. Please provide an annotated nucleotide GBK/EMBL file with a DNA sequence.\n")
        frame.update()
        sys.exit(1)
  #Extract genes
  genes = cdspart.split("FT   CDS             ")
  genes = genes[1:]
  genesdetails = parsegenes(genes, filetype, frame, outbox, GUI, nr)
  genelist = genesdetails[0]
  genedict = genesdetails[1]
  joinlist = genesdetails[2]
  joindict = genesdetails[3]
  accessiondict = genesdetails[4]
  locustagdict = genesdetails[5]
  #Locate all genes on DNA sequence and translate to protein sequence
  textlines = filetext.split("SQ   ")[0]
  textlines = textlines.split("\n")
  if accession == None:
      accession = ""
      for i in textlines:
        if accession == "":
          if "AC   " in i:
            j = i.split("AC   ")[1]
            j = j.replace(" ","")
            accession = j.split(";")[0]
            if len(accession) < 4:
              accession = "UNKNWN_" + str(nr)
  proteins = extractprotfasta(genelist,genedict,dnaseq,rc_dnaseq,joinlist,joindict,accessiondict,accession,locustagdict, frame, outbox, GUI)
  return [proteins, accession, dnaseqlength]

def gbk2proteins(gbkfile, nr=0, filetype="genbank_derived", frame=None, outbox=None, GUI="n", accession=None):
  filetext = gbkfile
  filetext = filetext.replace("\r","\n")
  if accession == None:
    accession = gbkfile.partition("ACCESSION   ")[2].partition("\n")[0].partition(" ")[0].replace(" ","")
  if "\n     CDS             " not in filetext or "\nORIGIN" not in filetext:
    if GUI == "n":
      print("No CDS or DNA sequence found in " + accession + ". Skipping file.")
      return "nocds"
    else:
      outbox.text_insert("No CDS or DNA sequence found in " + accession + ". Skipping file.\n")
      frame.update()
      return "nocds"
  cdspart = filetext.split("\nORIGIN")[0]
  #Extract DNA sequence and calculate reverse complement of it
  dnaseq = filetext.split("\nORIGIN")[1]
  dnaseq = cleandnaseq(dnaseq, frame, outbox, GUI)
  sequence = dnaseq
  if (sequence.count('N') + sequence.count('n') + sequence.count('A') + sequence.count('a') + sequence.count('C') + sequence.count('c') + sequence.count('G') + sequence.count('g') + sequence.count('T') + sequence.count('t')) < (0.5 * len(sequence)):
      if GUI == "n":
        print("Protein GBK/EMBL file provided. Please provide nucleotide GBK/EMBL file.")
        sys.exit(1)
      else:
        outbox.text_insert("Protein GBK/EMBL file provided. Please provide nucleotide GBK/EMBL file.\n")
        frame.update()
        sys.exit(1)
  dnaseqlength = len(dnaseq)
  rc_dnaseq = reverse_complement(dnaseq, frame, outbox, GUI)
  if dnaseqlength < 1:
      if GUI == "n":
        print("No sequence found in GBK/EMBL file. Please provide an annotated nucleotide GBK/EMBL file with a DNA sequence.")
        sys.exit(1)
      else:
        outbox.text_insert("No sequence found in GBK/EMBL file. Please provide an annotated nucleotide GBK/EMBL file with a DNA sequence.\n")
        frame.update()
        sys.exit(1)
  #Extract genes
  genes = cdspart.split("\n     CDS             ")
  genes = genes[1:]
  genesdetails = parsegenes(genes, filetype, frame, outbox, GUI, nr)
  genelist = genesdetails[0]
  genedict = genesdetails[1]
  joinlist = genesdetails[2]
  joindict = genesdetails[3]
  accessiondict = genesdetails[4]
  locustagdict = genesdetails[5]
  #Locate all genes on DNA sequence and translate to protein sequence
  textlines = filetext.split("\n//")[0]
  textlines = textlines.split("\n")
  if accession == None:
      accession = ""
      for i in textlines:
        if accession == "":
          if "ACCESSION   " in i:
            j = i.partition("ACCESSION   ")[2]
            accession = j.partition("\n")[0].partition(" ")[0]
            if len(accession) < 4:
              accession = "UNKWN" + str(nr)
  proteins = extractprotfasta(genelist,genedict,dnaseq,rc_dnaseq,joinlist,joindict,accessiondict, accession, locustagdict, frame, outbox, GUI)
  return [proteins,accession,dnaseqlength]

def check_files(seqfiles, dbname, outbox=None, GUI="n"):
  if GUI == "y":
    outbox.text_insert("Checking files for suitability: " + ", ".join([seqfile.rpartition(os.sep)[2].rpartition("/")[2] for seqfile in seqfiles]) + "\n")
  else:
    print("Checking files for suitability: " + ", ".join([seqfile.rpartition(os.sep)[2].rpartition("/")[2] for seqfile in seqfiles]))
  ##Check if files are present and check if they do not have the same accession numbers
  accessions = []
  for inputfilename in seqfiles:
    try:
      test_file = open(inputfilename,"rU")
    except:
      if GUI == "y":
        outbox.text_insert("File " + inputfilename + " not found.\n")
      else:
        print("File " + inputfilename + " not found.")
      clean_up(dbname)
      sys.exit(1)
    inputfile = open(inputfilename,"rU")
    ext = inputfilename.rpartition(".")[2]
    if ext.lower() in ["gbk","gb","genbank"]:
      infiletext = inputfile.read()
      if "LOCUS       " not in infiletext or infiletext.count("LOCUS       ") != infiletext.count("\nORIGIN") or "ORIGIN\n//" in infiletext.replace(" ",""):
        if not ("CONTIG      " in infiletext or "WGS         " in infiletext or "WGS_SCAFLD  " in infiletext):
          if GUI == "y":
            outbox.text_insert("Error: all GenBank entries should start with 'LOCUS' identifier and contain a sequence after 'ORIGIN'. Exiting.\n")
          else:
            print("Error: all GenBank entries should start with 'LOCUS' identifier and contain a sequence after 'ORIGIN'. Exiting.")
          clean_up(dbname)
          sys.exit(1)
      entries = infiletext.split("LOCUS       ")[1:]
      entrynr = 0
      for i in entries:
        if "DEFINITION" not in i:
            print("Warning: entry in " + inputfilename + " lacks GBK definition line; using file name instead.")
            description = inputfilename.rpartition(os.sep)[2].replace(".","_") + "_" + str(entrynr)
        else:
            definitionlines = i.partition("DEFINITION  ")[2]
            description = definitionlines.partition("\n")[0]
            definitionlines = definitionlines.partition("\n")[2]
            while " " == definitionlines.partition("\n")[0][0]:
                description += definitionlines.partition("\n")[0]
                definitionlines = definitionlines.partition("\n")[2]
            description = description.replace("\n"," ").replace("            "," ").replace("  "," ").replace("  "," ")
        accession = i.partition("ACCESSION   ")[2].partition("\n")[0].partition(" ")[0].replace(" ","")
        if len(accession) == 0:
          if GUI == "y":
            outbox.text_insert("Warning: entry in " + inputfilename + " lacks accession number; using file name instead.\n")
          else:
            print("Warning: entry in " + inputfilename + " lacks accession number; using file name instead.")
          accession = inputfilename.rpartition(os.sep)[2].replace(".","_") + "_" + str(entrynr)
          accessions.append(accession)
        else:
          accessions.append(accession)
        entrynr += 1
    elif ext.lower() in ["embl","emb"]:
      infiletext = inputfile.read()
      if "ID   " not in infiletext or infiletext.count("ID   ") != infiletext.count("SQ   ") or "other;\n//" in infiletext.replace(" ",""):
        if GUI == "y":
          outbox.text_insert("Error: all EMBL entries should start with 'ID' identifier and contain a sequence after 'SQ'. Exiting.\n")
        else:
          print("Error: all EMBL entries should start with 'ID' identifier and contain a sequence after 'SQ'. Exiting.")
        clean_up(dbname)
        sys.exit(1)
      entries = infiletext.split("ID   ")[1:]
      entrynr = 0
      for i in entries:
        description = i.partition("DE   ")[2].partition("\n")[0].replace(";","").replace("\n"," ")
        accession = i.partition("AC   ")[2].partition("\n")[0].replace(";","").partition(" ")[0].replace(" ","")
        if len(accession) == 0:
          if GUI == "y":
            outbox.text_insert("Warning: entry in " + inputfilename + " lacks accession number; using file name and entry nr instead.\n")
          else:
            print("Warning: entry in " + inputfilename + " lacks accession number; using file name and entry nr  instead.")
          accession = inputfilename.rpartition(os.sep)[2].replace(".","_") + "_" + str(entrynr)
          accessions.append(accession)
        else:
          accessions.append(accession)
      entrynr += 1
    inputfile.close()

def parse_gbk_embl(args, dbname, frame=None, outbox=None, GUI="n"):
  seqfiles = args
  check_files(seqfiles, dbname, outbox, GUI)
  ##Parse input files
  outputfile = open(dbname + "_dbbuild.fasta","w")
  #For all files
  nr = 1
  allnames = []
  all_accessions = []
  allseqs = []
  descriptions = {}
  for seqfile in seqfiles:
    if GUI == "y":
      outbox.text_insert("  Incorporating " + seqfile + "\n")
      frame.update()
    if GUI == "n":
      if seqfile not in os.listdir(".") and seqfile.rpartition(os.sep)[2] not in os.listdir(seqfile.rpartition(os.sep)[0]):
        outbox.text_insert("File " + seqfile + " not found. Skipping..\n")
        frame.update()
        continue
      else:
        print("  Incorporating " + seqfile + "\n")
    #Open seq file  
    infile = open(seqfile,"rU").read()
    if GUI == "y":
      frame.update()
    ext = seqfile.rpartition(".")[2]
    entries = []
    if ext.lower() in ["gbk","gb","genbank"]:
      #Divide in GBK entries
      tempgbks = infile.split("LOCUS       ")[1:]
      #For each entry, check if >2 CDS features.
      for i in tempgbks:
        entry = "LOCUS       " + i
        if "ORIGIN      \n        1 " not in i and ("\nWGS_SCAFLD  " in entry or "\nWGS         " in entry):
          gbks = fix_wgs_master_record(entry, dbname, frame, outbox, GUI)
          for gbk in gbks:
            if gbk.count("     CDS           ") >= 1:
              entries.append(gbk)
        elif "ORIGIN      \n        1 " not in i and "\nCONTIG      " in entry:
          entry = fix_supercontig_record(entry, dbname, frame, outbox, GUI)
          if entry.count("     CDS           ") >= 1:
            entries.append(entry)
        elif entry.count("     CDS           ") >= 1:
          entries.append(entry)
    elif ext.lower() in ["embl","emb"]:
      #Divide in EMBL entries
      tempembls = infile.split("ID   ")[1:]
      #For each entry, check if >2 CDS features.
      for i in tempembls:
        entry = "ID   " + i
        if entry.count("FT   CDS             ") >= 1:
          entries.append(entry)
    #For all resultant GBKs, parse sequence
    entrynr = 1
    for infile in entries:
      if GUI == "y":
        frame.update()
      if ext.lower() in ["gbk","gb","genbank"]:
        accession, all_accessions = get_accession(infile, entrynr, "gbk", all_accessions, seqfile, dbname)
        accession, sequence, description = get_sequence(infile, allnames, seqfile, entrynr, "gbk", dbname)
        proteins = gbk2proteins(infile, nr, "user_inputfile", frame, outbox, GUI, accession)
      elif ext.lower() in ["embl","emb"]:
        accession, all_accessions = get_accession(infile, entrynr, "embl", all_accessions, seqfile, dbname)
        accession, sequence, description = get_sequence(infile, allnames, seqfile, entrynr, "embl", dbname)
        proteins = embl2proteins(infile, nr, "user_inputfile", frame, outbox, GUI, accession)
      if proteins != "nocds":
        genomic_accnr = proteins[1]
        dnaseqlength = proteins[2]
        proteins = proteins[0]
        allnames = allnames + proteins[0]
        allseqs = allseqs + proteins[1]
      descriptions[accession] = description
      entrynr += 1
  if GUI == "y":
    frame.update()

    #TODO this needs a change for sure
  allnames = [name.split("|")[3] for name in allnames]
  writefasta(allnames,allseqs,outputfile)
  if len(allnames) == 0:
    if GUI == "y":
      outbox.text_insert("ERROR: No CDS features found in these entries. Exiting.\n")
      frame.update()
      sys.exit(1)
    else:
      print("ERROR: No CDS features found in these entries. Exiting.")
      clean_up(dbname)
      sys.exit(1)
  outputfile.close()
  return descriptions

def is_nucl_seq(sequence):
    if len(str(sequence).lower().replace("a","").replace("c","").replace("g","").replace("t","").replace("n","")) < 0.2 * len(sequence):
        return True
    else:
        return False

def checkfastafile(sequence, frame=None, outbox=None, GUI="n"):
    #Check input file formatting
    if not is_nucl_seq(sequence):
        if GUI == "y":
            outbox.text_insert("Warning: protein FASTA file provided. Entry skipped.\n")
            frame.update()
        else:
            print("Warning: protein FASTA file provided. Entry skipped.\n")
        return False
    return True

def nucparse_gbk_embl_fasta(args, dbname, frame=None, outbox=None, GUI="n"):
  seqfiles = args
  outputfile = open(dbname + "_dbbuild.fasta","w")
  #For all files
  allnames = []
  allseqs = []
  descriptions = {}
  for seqfile in seqfiles:
    if GUI == "y":
      outbox.text_insert("  Incorporating " + seqfile + "\n")
      frame.update()
    if GUI == "n":
      if seqfile not in os.listdir(".") and seqfile.rpartition(os.sep)[2] not in os.listdir(seqfile.rpartition(os.sep)[0]):
        outbox.text_insert("File " + seqfile + " not found. Skipping..\n")
        frame.update()
        continue
      else:
        print("  Incorporating " + seqfile + "\n")
    #Open seq file  
    infile = open(seqfile,"r").read()
    if GUI == "y":
      frame.update()
    ext = seqfile.rpartition(".")[2]
    entries = []
    if ext.lower() in ["gbk","gb","genbank"]:
      #Divide in GBK entries
      tempgbks = infile.split("LOCUS       ")[1:]
      for entry in tempgbks:
        if "\nWGS_SCAFLD  " in entry or "\nWGS         " in entry:
          gbks = fix_wgs_master_record(entry, dbname, frame, outbox, GUI)
          for gbk in gbks:
            entries.append(gbk)
        elif "\nCONTIG      " in entry:
          entries.append(fix_supercontig_record(entry, dbname, frame, outbox, GUI))
        else:
          entries.append("LOCUS       " + entry)
    elif ext.lower() in ["embl","emb"]:
      #Divide in EMBL entries
      tempembls = infile.split("ID   ")[1:]
      #For each entry, check if >2 CDS features.
      for entry in tempembls:
        entries.append("ID   " + entry)
    elif ext.lower() in ["fasta","fas","fa"]:
      #Divide in FASTA entries
      tempfastas = [">" + fasta for fasta in infile.split(">")[1:]]
      for entry in tempfastas:
        if "\n" in entry and len(entry.partition("\n")[2]) > 3:
          sequence = entry.partition(">")[2].partition("\n")[2]
          if checkfastafile(sequence, frame, outbox, GUI):
            entries.append(">" + entry)
    if len(entries) == 0:
      if GUI == "n":
        print("Error: protein FASTA files given instead of nucleotide FASTA files. Exiting.")
        clean_up(dbname)
        sys.exit(1)
      else:
        outbox.text_insert("Error: empty FASTA files or protein FASTA files given instead of nucleotide FASTA files. Exiting.\n")
        frame.update()
        sys.exit(1)
    #For all resultant GBKs, parse sequence
    entrynr = 1
    for entry in entries:
      if GUI == "y":
        frame.update()
      if ext.lower() in ["gbk","gb","genbank"]:
        accession, sequence, description = get_sequence(entry, allnames, seqfile, entrynr, "gbk", dbname)
      elif ext.lower() in ["embl","emb"]:
        accession, sequence, description = get_sequence(entry, allnames, seqfile, entrynr, "embl", dbname)
      elif ext.lower() in ["fasta","fas","fa"]:
        accession, sequence, description = get_sequence(entry, allnames, seqfile, entrynr, "fasta", dbname)
      if accession in allnames:
        if GUI == "n":
          print("Error: multiple input files with same name / accession number:", accession)
          clean_up(dbname)
          sys.exit()
        else:
          outbox.text_insert("Error: multiple input files with same name / accession number: " +  accession + "\n")
          frame.update()
          sys.exit(1)
      else:
          descriptions[accession] = description
          allnames.append(accession)
          allseqs.append(sequence)
      entrynr += 1
  if GUI == "y":
    frame.update()
  writefasta(allnames,allseqs,outputfile)
  outputfile.close()
  with open(dbname + "_dbbuild.fasta","r") as f:
    print(f.read())
  return descriptions

def get_sequence(entry, allnames, seqfile, entrynr, entrytype, dbname):
  seqfilename = seqfile.rpartition(os.sep)[2].replace(".","_")
  if entrytype == "gbk":
      if "DEFINITION" not in entry:
          description = seqfile.rpartition(os.sep)[2].replace(".","_") + "_" + str(entrynr)
      else:
          definitionlines = entry.partition("DEFINITION  ")[2]
          description = definitionlines.partition("\n")[0]
          definitionlines = definitionlines.partition("\n")[2]
          while " " == definitionlines.partition("\n")[0][0]:
              description += definitionlines.partition("\n")[0]
              definitionlines = definitionlines.partition("\n")[2]
          description = description.replace("\n"," ").replace("            "," ").replace("  "," ").replace("  "," ")
      accession, accessions = get_accession(entry, entrynr, entrytype, [], seqfilename, dbname)
      sequence = cleandnaseq(entry.partition("\nORIGIN")[2].partition("\n")[2])
  elif entrytype == "embl":
      description = entry.partition("DE   ")[2].partition("\n")[0].replace(";","")
      accession, accessions = get_accession(entry, entrynr, entrytype, [], seqfilename, dbname)
      sequence = cleandnaseq(entry.partition("SQ   ")[2].partition("\n")[2])
  elif entrytype == "fasta":
      accession = entry.partition(">")[2].partition("\n")[0].replace("|","_")
      description = accession
      sequence = entry.partition(">")[2].partition("\n")[2]
  return accession, sequence, description

#def cleandnaseq(dnaseq, frame=None, outbox=None, GUI="n"):
#    "Remove all whitespace from sequence and replace invalid chars with 'n'"
#    dnaseq = ''.join(c for c in dnaseq if c not in string.whitespace)
#
#    # Only allow characters in 'nucleotides', replace everything else
#    nucleotides = ["a", "c", "g", "t", "n"]
#    newseq = ''.join(c if c.lower() in nucleotides else 'n' for c in dnaseq)
#    return newseq


def fetch_entries_from_ncbi(efetch_url):
    urltry = "n"
    nrtries = 0
    output = ""
    while urltry == "n" and nrtries < 4:
        try:
            nrtries += 1
            time.sleep(3)
            req = urllib.request.Request(efetch_url)
            response = urllib.request.urlopen(req)
            output = response.read()
            if len(output) > 5:
                urltry = "y"
        except (IOError,http.client.BadStatusLine,URLError,http.client.HTTPException):
            print("Entry fetching from NCBI failed. Waiting for connection...")
            time.sleep(5)
    return output


def fix_wgs_master_record(gbk, dbname, frame=None, outbox=None, GUI="n"):
    #If seq_record is a WGS master record, parse out contig accession numbers and download these as separate seq_records
    if "WGS_SCAFLD  " in gbk:
        contigranges = gbk.rpartition("WGS_SCAFLD  ")[2].partition("\n")[0]
    elif "WGS         " in gbk:
        contigranges = gbk.rpartition("WGS         ")[2].partition("\n")[0]
    if ";" in contigranges:
        contigranges2 = []
        for contigrange in contigranges.split(";"):
            if "-" in contigrange:
                contigranges2.append(contigrange.split("-"))
            else:
                contigranges2.append([contigrange])
        contigranges = contigranges2
    elif "-" in contigranges:
        contigranges = [contigranges.split("-")]
    else:
        contigranges = [[contigranges]]
    allcontigs = []
    for contigrange in contigranges:
        if len(contigrange) == 1:
            allcontigs.extend(contigrange)
            continue
        startnumber, endnumber = '', ''
        alpha_tag = ''
        for char in contigrange[0].partition(".")[0]:
            if char.isdigit():
                startnumber += char
            else:
                alpha_tag += char
        for char in contigrange[1].partition(".")[0]:
            if char.isdigit():
                endnumber += char
        nrzeros = 0
        for char in startnumber:
            if char == "0":
                nrzeros += 1
            else:
                break
        contigrange = [alpha_tag + nrzeros * "0" + str(number) for number in range(int(startnumber), int(endnumber))]
        allcontigs.extend(contigrange)
    #Create contig groups of 50 (reasonable download size per download)
    nr_groups = len(allcontigs) / 50 + 1
    contig_groups = [allcontigs[i:i+50] for i in range(0, len(allcontigs), 50)]
    #Download contigs and parse into seq_record objects
    allgbks = []
    for contig_group in contig_groups:
        efetch_url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id='
        efetch_url = efetch_url + ",".join(contig_group) + '&rettype=gbwithparts&retmode=text'
        output = fetch_entries_from_ncbi(efetch_url)
        if not len(output) > 5:
            break
        if "Resource temporarily unavailable" in output[:200] or "<h1>Server Error</h1>" in output[:500] or "NCBI - WWW Error" in output[:500] or "LOCUS       " not in output[:500]:
            if GUI == "n":
                print('ERROR: NCBI server temporarily unavailable: downloading contigs failed.')
                clean_up(dbname)
                sys.exit(1)
            else:
                outbox.text_insert("ERROR: NCBI server temporarily unavailable: downloading contigs failed.\n")
                frame.update()
                sys.exit(1)
        gbks = ["LOCUS       " + locus for locus in output.split("LOCUS       ")[1:]]
        allgbks.extend(gbks)
    return allgbks

def fix_supercontig_record(gbk, dbname, frame=None, outbox=None, GUI="n"):
    #If seq_record is a supercontig record, reconstruct sequence and replace CONTIG feature by ORIGIN feature
    if "CONTIG      " in gbk:
        contig_info = gbk.rpartition("CONTIG      ")[2].partition("\n//")[0].replace("join(","").replace("\n","")
    if contig_info[-1] == ")":
        contig_info = contig_info[:-1]
    allcontigparts = [part.replace(" ","") for part in contig_info.split(",")]
    accessions = [part.partition(":")[0].partition(".")[0].replace(" ","") for part in contig_info.split(",") if "gap(" not in part]
    #Create contig groups of 50 (reasonable download size per download)
    contig_groups = [accessions[i:i+50] for i in range(0, len(accessions), 50)] 
    #Download contig sequences based on their accessions
    contigseqdict = {}
    for contig_group in contig_groups:
        efetch_url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id='
        efetch_url = efetch_url + ",".join(contig_group) + '&rettype=fasta&retmode=text'
        output = fetch_entries_from_ncbi(efetch_url)
        if not len(output) > 5:
            break
        if "Resource temporarily unavailable" in output[:200] or "<h1>Server Error</h1>" in output[:500]:
            if GUI == "n":
                print('ERROR: NCBI server temporarily unavailable: downloading contigs failed.')
                clean_up(dbname)
                sys.exit(1)
            else:
                outbox.text_insert("ERROR: NCBI server temporarily unavailable: downloading contigs failed.\n")
                frame.update()
                sys.exit(1)
        sequences = [seq for seq in output.split(">") if len(seq) > 5]
        for sequence in sequences:
            for contig_acc in contig_group:
                if contig_acc in sequence.partition("\n")[0]:
                    contigseqdict[contig_acc] = sequence.partition("\n")[2].replace("\n","")
    #Reconstruct supercontig sequence based on contig sequences and gaps
    fullsequence = ''
    for part in allcontigparts:
        if "gap(" in part:
            candidate_gap_int = part.partition('gap(')[2][:-1]
            if "unk" in candidate_gap_int:
                candidate_gap_int = candidate_gap_int.partition("unk")[2]
            if candidate_gap_int.isdigit():
                fullsequence += int(candidate_gap_int) * 'N'
            else:
                print('Parsing supercontig file failed: faulty gap identifier' + part)
                sys.exit(1)
        else:
            accession = part.partition(":")[0].partition(".")[0]
            sequence = contigseqdict[accession]
            if ":" in part and ".." in part:
                seqstart, seqend = part.partition(":")[2].split("..")
                if int(seqstart) > 0:
                    seqstart = int(seqstart) - 1
                sequence = sequence[int(seqstart) : int(seqend)]
            fullsequence += sequence
    #Add compiled sequence to seq_record
    newgbk = ""
    newgbk += gbk.partition("CONTIG      ")[0] + "ORIGIN      \n"
    x = 1
    while len(fullsequence) > 0:
        seqpart = fullsequence[:60].lower()
        fullsequence = fullsequence[60:]
        nrspaces = 9 - len(str(x))
        if len(seqpart) > 50:
            newgbk += nrspaces * " " + str(x) + " " + " ".join([seqpart[:10], seqpart[10:20], seqpart[20:30], seqpart[30:40], seqpart[40:50], seqpart[50:60]]) + "\n"
        elif len(seqpart) > 40:
            newgbk += nrspaces * " " + str(x) + " " + " ".join([seqpart[:10], seqpart[10:20], seqpart[20:30], seqpart[30:40], seqpart[40:50]]) + "\n"
        elif len(seqpart) > 30:
            newgbk += nrspaces * " " + str(x) + " " + " ".join([seqpart[:10], seqpart[10:20], seqpart[20:30], seqpart[30:40]]) + "\n"
        elif len(seqpart) > 20:
            newgbk += nrspaces * " " + str(x) + " " + " ".join([seqpart[:10], seqpart[10:20], seqpart[20:30]]) + "\n"
        elif len(seqpart) > 10:
            newgbk += nrspaces * " " + str(x) + " " + " ".join([seqpart[:10], seqpart[10:20]]) + "\n"
        else:
            newgbk += nrspaces * " " + str(x) + " " + seqpart + "\n"
        x += 60
    newgbk += "//"
    return newgbk