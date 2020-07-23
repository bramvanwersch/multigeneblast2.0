#!/usr/bin/env python
## Copyright (c) 2010 Marnix H. Medema
## University of Groningen
## Department of Microbial Physiology / Groningen Bioinformatics Centre

###INSTRUCTIONS:
###COMPATIBLE WITH BOTH LINUX AND WINDOWS
###TO RUN BATCH ANTISMASH RUNS, PLACE THIS SCRIPT IN THE FOLDER TOGETHER WITH A FILE 'annotationtable.txt' CONTAINING YOUR ANNOTATION TABLE.
###USAGE: python format_embl.py
#The table in annotationtable.txt should contain the following columns, which should be ordered in the order in which they appear on the DNA.
#1. FASTA file name of contig (which should be placed in same folder)
#2. Gene locus tag
#3. 5' gene or exon start
#4. 3' gene or exon stop
#5. Gene function / annotation

##Script for formatting an annotation table as EMBL format.
def integerfromstring(string):
  numbers = ['0','1','2','3','4','5','6','7','8','9']
  newstring = ""
  for i in string:
    if i in numbers:
      newstring = newstring + i
  return int(newstring)

import sys
import operator
import itertools

#Read in annotationtable.txt
try:
  annotationtable = [i.split("\t") for i in open("annotationtable.txt","r").read().split("\n") if len(i.split("\t")) > 3]
except:
  print("""Please create a TXT file 'annotationtable.txt', with CDS info for your contigs.
The table in 'annotationtable.txt' should contain the following columns, which should be ordered in the order in which they appear on the DNA.
1. FASTA file name of contig (which should be placed in same folder)
2. Gene locus tag
3. 5' gene or exon start
4. 3' gene or exon stop
5. Gene function / annotation""")
  sys.exit()
print("Reading annotations")

#Quality check for unicity of locus tags in annotation table
lastitem = ""
unique_items = []
for i in [j[0] for j in annotationtable]:
  if i == lastitem:
    continue
  elif i in unique_items:
    print("Not all locus tags in your table are unique for a gene. Please correct this.")
    sys.exit()
  else:
    unique_items.append(i)
  lastitem = i
#Quality check for real values in 3nd and 4rd columns in annotation table
annotationtable2 = []
for i in annotationtable:
  if not (i[2].isdigit() and i[3].isdigit()):
    continue
  else:
    annotationtable2.append(i)
annotationtable = annotationtable2

#Reading sequence files
sequencefiles = [i[0] for i in annotationtable]
sequencedict = {}
for sequencefile in sequencefiles:
  try:
    seq = open(sequencefile,"r").read().partition(">")[2].partition("\n")[2].replace(" ","").replace("\n","").replace("\t","")
  except:
    print("Sequence file ", sequencefile, "not found. Please modify annotationtable.txt")
    sys.exit()
  if ">" in seq:
    print(sequencefile,"MultiFASTA sequence file provided. Please provide a single DNA sequence FASTA file.")
  dnaseqlength = len(seq)
  sequencedict[sequencefile] = [seq.lower(), dnaseqlength]

#Split annotation table in multiple pieces, each representing a contig / chromosome
x = 1
annotationtablenames = []
annotationtables = {}
contigtable = []
for entry in annotationtable:
  #For first entry
  if x == 1:
    contigtable.append(entry[1:])
    lastcontig = entry[0]
    if x == len(annotationtable):
      annotationtables[entry[0]] = contigtable
      annotationtablenames.append(entry[0])
  #For last entry
  elif x == len(annotationtable):
    if lastcontig != entry[0]:
      annotationtables[lastcontig] = contigtable
      annotationtablenames.append(lastcontig)
      contigtable = [entry[1:]]
      annotationtables[entry[0]] = contigtable
      annotationtablenames.append(entry[0])
    else:
      contigtable.append(entry[1:])
      annotationtables[entry[0]] = contigtable
      annotationtablenames.append(entry[0])
  #For middle entries
  else:
    if lastcontig == entry[0]:
      contigtable.append(entry[1:])
    else:
      annotationtables[lastcontig] = contigtable
      annotationtablenames.append(lastcontig)
      contigtable = [entry[1:]]
      lastcontig = entry[0]
  x += 1

#For each contig, create an EMBL file
for annotationtablename in annotationtablenames:
  seq, dnaseqlength = sequencedict[annotationtablename]
  annotationtable = annotationtables[annotationtablename]
  orfnames = []
  positions = []
  firstline = "y"
  x = 0
  starts = []
  ends = []
  annotations = []
  for i in annotationtable:
    columns = i
    if len(i) > 1:
      if x == 0:
        if int(columns[1]) < int(columns[2]):
          strand = "+"
        else:
          strand = "-"
        starts.append(str(min([int(columns[1]), int(columns[2])])))
        ends.append(str(max([int(columns[1]), int(columns[2])])))
        #Save orf if this is the last line describing the first orf
        if len(annotationtable) == 1 or annotationtable[x + 1][0] != columns[0]:
          orfnames.append(columns[0])
          annotations.append(columns[3])
          if strand == "+":
            if len(starts) == 1:
              if starts[0] == '0':
                starts[0] = '1'
              if ends[0] == '0':
                ends[0] = '1'
              pos = starts[0] + ".." + ends[0]
              positions.append(pos)
            else:
              pos = "join("
              y = 0
              for i in starts:
                if i == '0':
                  i = '1'
                if ends[y] == '0':
                  ends[y] = '1'
                pos = pos + i + ".." + ends[y]
                if i != starts[-1]:
                  pos = pos + ","
                y += 1
              pos = pos + ")"
              positions.append(pos)
          elif strand == "-":
            if len(starts) == 1:
              if starts[0] == '0':
                starts[0] = '1'
              if ends[0] == '0':
                ends[0] = '1'
              pos = "complement(" + starts[0] + ".." + ends[0] + ")"
              positions.append(pos)
            else:
              pos = "complement(join("
              y = 0
              for i in starts:
                if i == '0':
                  i = '1'
                if ends[y] == '0':
                  ends[y] = '1'
                pos = pos + i + ".." + ends[y]
                if i != starts[-1]:
                  pos = pos + ","
                y += 1
              pos = pos + "))"
              positions.append(pos)
          starts = []
          ends = []
      #Save orf if this is the last line describing the orf
      elif x == (len(annotationtable) - 1) or (columns[0] != annotationtable[x + 1][0] and annotationtable[x + 1][0] not in orfnames):
        if int(columns[1]) < int(columns[2]):
          strand = "+"
        else:
          strand = "-"
        starts.append(str(min([int(columns[1]), int(columns[2])])))
        ends.append(str(max([int(columns[1]), int(columns[2])])))
        orfnames.append(columns[0])
        annotations.append(columns[3])
        if strand == "+":
          if len(starts) == 1:
            if starts[0] == '0':
              starts[0] = '1'
            if ends[0] == '0':
              ends[0] = '1'
            pos = starts[0] + ".." + ends[0]
            positions.append(pos)
          else:
            pos = "join("
            y = 0
            starts = [int(startpos) for startpos in starts]
            starts.sort()
            starts = [str(startpos) for startpos in starts]
            ends = [int(endpos) for endpos in ends]
            ends.sort()
            ends = [str(endpos) for endpos in ends]
            for i in starts:
              if i == '0':
                i = '1'
              if ends[y] == '0':
                ends[y] = '1'
              pos = pos + i + ".." + ends[y]
              if i != starts[-1]:
                pos = pos + ","
              y += 1
            pos = pos + ")"
            positions.append(pos)
        elif strand == "-":
          if len(starts) == 1:
            if starts[0] == '0':
              starts[0] = '1'
            if ends[0] == '0':
              ends[0] = '1'
            pos = "complement(" + starts[0] + ".." + ends[0] + ")"
            positions.append(pos)
          else:
            pos = "complement(join("
            y = 0
            starts = [int(startpos) for startpos in starts]
            starts.sort()
            starts = [str(startpos) for startpos in starts]
            ends = [int(endpos) for endpos in ends]
            ends.sort()
            ends = [str(endpos) for endpos in ends]
            for i in starts:
              if i == '0':
                i = '1'
              if ends[y] == '0':
                ends[y] = '1'
              pos = pos + i + ".." + ends[y]
              if i != starts[-1]:
                pos = pos + ","
              y += 1
            pos = pos + "))"
            positions.append(pos)
        starts = []
        ends = []
      else:
        starts.append(str(min([int(columns[1]), int(columns[2])])))
        ends.append(str(max([int(columns[1]), int(columns[2])])))
    x += 1
  out_file = open(annotationtablename.rpartition(".")[0] + '.embl',"w")
  a = 0
  print("Writing EMBL file based on",annotationtablename,"...")
  out_file.write("ID   A01; SV 1; linear; DNA; STD; FUN; " + str(dnaseqlength) + " BP.\nXX\n")
  out_file.write("AC   A01;\nXX\n")
  out_file.write("DE   unknown organism;\nXX\n")
  out_file.write("KW   none;\nXX\n")
  out_file.write("OS   unknown;\n")
  out_file.write("OC   unknown;\nXX\n")
  out_file.write("RN   [1]\n")
  out_file.write("RT   ;\n")
  out_file.write("RL   Unknown.\nXX\n")
  out_file.write("FH   Key             Location/Qualifiers\nFH\n")
  out_file.write("FT   source          1.." + str(dnaseqlength) + "\n")
  for i in orfnames:
    if abs(integerfromstring(positions[a].partition(",")[0].partition("..")[0]) - integerfromstring(positions[a].rpartition(",")[2].rpartition("..")[2])) < 300000:
      out_file.write("FT   gene            ")
      out_file.write(positions[a])
      out_file.write("\n")
      out_file.write('FT                   /gene="' + i + '"\n')
      out_file.write("FT   CDS             ")
      out_file.write(positions[a])
      out_file.write("\n")
      out_file.write('FT                   /gene="' + i + '"\n')
    a += 1
  out_file.write("XX\nSQ   Sequence " + str(dnaseqlength) + " BP; " + str(seq.count("a") + seq.count("A")) + " A; " + str(seq.count("c") + seq.count("C")) + " C; " + str(seq.count("g") + seq.count("G")) + " G; " + str(seq.count("t") + seq.count("T")) + " T; " + str(dnaseqlength - (seq.count("a") + seq.count("A") + seq.count("c") + seq.count("C") + seq.count("g") + seq.count("G") + seq.count("t") + seq.count("T"))) + " other;\n")
  seq2 = seq
  out_file.write("     ")
  grouplen=10
  textlen = len(seq)
  end = textlen - (textlen % grouplen)
  repeated_iterator = [iter(itertools.islice(seq, 0, end))] * grouplen
  parts = list(map(lambda *chars: ''.join(chars),*repeated_iterator))
  if dnaseqlength%grouplen != 0:
    parts.append(seq[-1 * (dnaseqlength%grouplen):])
  w = 1
  for l in parts:
    out_file.write(l + " ")
    if w == len(parts):
      if w%6 == 0 and dnaseqlength%60 != 0:
        out_file.write((" " * (10 - dnaseqlength%grouplen) + " " * (10 - len(str(dnaseqlength)))) + str(dnaseqlength) + "\n//")
      elif dnaseqlength%60 == 0:
        out_file.write((" " * (10 - len(str(dnaseqlength)))) + str(dnaseqlength) + "\n//")
      elif w%6 == 5 and dnaseqlength%grouplen == 0:
        out_file.write(("           " + " " * (10 - len(str(dnaseqlength)))) + str(dnaseqlength) + "\n//")
      elif dnaseqlength%grouplen != 0:
        out_file.write(" " * (10 - dnaseqlength%grouplen) + "          " * (6 - len(parts)%6) + " " * (6 - len(parts)%6) + (" " * (10 - len(str(dnaseqlength)))) + str(dnaseqlength) + "\n//")
      else:
        out_file.write("          " * (6 - len(parts)%6) + " " * (5 - len(parts)%6) + (" " * (10 - len(str(dnaseqlength)))) + str(dnaseqlength) + "\n//")
    elif w%6 == 0:
      out_file.write((" " * (10 - len(str(w * 10)))) + str(w * 10) + "\n     ")
    w += 1
  out_file.close()