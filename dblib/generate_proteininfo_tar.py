import pickle as pickle
import sys
import os
import tarfile
import shutil
import fileinput
from .utils import log, clean_up

def sortdictkeysbyvalues(dict):
    items = [(len(list(dict[key].keys())), key) for key, value in list(dict.items())]
    items.sort()
    return [key for value, key in items]

def generate_proteininfo_tar(dbname, frame=None, outbox=None, GUI="n"):
  if GUI == "n":
    print("Generating proteininfo TAR archive")
    log("Generating proteininfo TAR archive")
  #Create proteininfo directory if nonexistent
  try:
    os.mkdir("proteininfo")
    print(os.getcwd())
    print("trigger?")
  except:
    pass

  #Read protein info input from genbank_mf_all/txt;
  #In the mean time, start saving to pickle files in order to reduce memory load
  proteininfodict= {}
  z = 0
  for i in fileinput.input(dbname + "_all.txt"):
    if GUI == "y":
      frame.update()
    if len(i) > 0:
      if z < 10000:
        i = i.replace(">","").replace("\n","")
        tabs = i.split("^")
        protein = tabs[3]
        genome = tabs[0]
        illegal_chars  = '''!"#$%&()*+,:;=>?@[]^`'{|} '''
        tag = "".join([char for char in protein[:4] if char not in illegal_chars]).upper()
        if tag in proteininfodict:
            if protein in proteininfodict[tag]:
              #fileinput.close()
              #clean_up(dbname, exit=True)
              print("Warning: non-unique protein accession:" + protein + "\n")
              log("Warning: non-unique protein accession:" + protein + "\n", exit=False)
            proteininfodict[tag][protein] = i
        else:
          proteininfodict[tag] = {}
          proteininfodict[tag][protein] = i
        z += 1
      else:
        first = "no"
        i = i.replace(">","").replace("\n","")
        tabs = i.split("|")
        protein = tabs[3]
        genome = tabs[0]
        illegal_chars  = '''!"#$%&()*+,:;=>?@[]^`'{|} '''
        tag = "".join([char for char in protein[:4] if char not in illegal_chars]).upper()
        if tag in proteininfodict:
            if protein in proteininfodict[tag]:
              #fileinput.close()
              #clean_up(dbname, exit=True)
              print("Warning: non-unique protein accession:" + protein + "\n")
              log("Warning: non-unique protein accession:" + protein + "\n", exit=False)
            proteininfodict[tag][protein] = i
        else:
          proteininfodict[tag] = {}
          proteininfodict[tag][protein] = i
        for j in list(proteininfodict.keys()):
          if GUI == "y":
            frame.update()
          try:
            picklefile = open("proteininfo" + os.sep + j + ".pickle","rb")
            tagdict = pickle.load(picklefile)
            picklefile.close()
            tempdict = dict(list(proteininfodict[j].items()) + list(tagdict.items()))
            outfile = open("proteininfo" + os.sep + j + ".pickle","wb")
            pickle.dump(tempdict,outfile)
            outfile.close()
          except(IOError):
            outfile = open("proteininfo" + os.sep + j + ".pickle","wb")
            pickle.dump(proteininfodict[j],outfile)
            outfile.close()
        proteininfodict = {}
        z = 0

  #Sort dictionary by size
  sortedproteininfokeylist = sortdictkeysbyvalues(proteininfodict)

  #Pickle the last pieces into the proteininfo folder
  for j in sortedproteininfokeylist:
    if GUI == "y":
      frame.update()
    try:
      picklefile = open("proteininfo" + os.sep + j + ".pickle","rb")
      tagdict = pickle.load(picklefile)
      picklefile.close()
      tempdict = dict(list(proteininfodict[j].items()) + list(tagdict.items()))
      outfile = open("proteininfo" + os.sep + j + ".pickle","wb")
      pickle.dump(tempdict,outfile)
      outfile.close()
      del proteininfodict[j]
    except(IOError):
      outfile = open("proteininfo" + os.sep + j + ".pickle","wb")
      pickle.dump(proteininfodict[j],outfile)

      outfile.close()
      del proteininfodict[j]

  #Archive directory as TAR file and remove original directory
  # with open("proteininfo" + os.sep + j + ".pickle", "r") as f:
  #   print(f.read())
  if GUI == "y":
    frame.update()
  try:
    tar = tarfile.open(dbname + ".pinfo.tar", "w")
    tar.add("proteininfo")
    tar.close()
  except:
    print("Could not create TAR file from genecords folder. Please create archive manually.")
    log("Could not create TAR file from genecords folder. Please create archive manually.", exit=True)
  if GUI == "y":
    frame.update()
  shutil.rmtree("proteininfo")
