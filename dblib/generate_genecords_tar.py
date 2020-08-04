import pickle as pickle
import sys
import os
import tarfile
import shutil
import fileinput
from .utils import log

def sortdictkeysbyvalues(dict):
    items = [(sum([len(dict[key][subkey]) for subkey in list(dict[key].keys())]), key) for key, value in list(dict.items())]
    items.sort()
    return [key for value, key in items]

def save_to_pickle(tag, genome_info_dict):
  #Dump info to pickle files in genecords directory
    #print sum([len(genecordslistdict[j][subkey]) for subkey in genecordslistdict[j].keys()]) ##Can be used to check size of dict contents
    outfile = open("genecords" + os.sep + tag + ".pickle","wb")
    try:
      picklefile = pickle.Pickler(outfile)
      picklefile.fast = True
      picklefile.dump(genome_info_dict)
    except:
      print("Error")
      sys.exit()
    outfile.close()

def generate_genecords_tar(dbname, frame=None, outbox=None, GUI="n"):
  if GUI == "n":    
    print("Generating genecords TAR archive")
    log("Generating genecords TAR archive")
  #Create genecords directory if nonexistent
  try:
    os.mkdir("genecords")
  except:
    pass
  #genecordslistdict = {}

  #Read gene coordinates input from genbank_mf_all/txt
  passedtags = []
  genome_info = {}
  lasttag = ""
  for i in fileinput.input(dbname + "_all.txt"):
    if GUI == "y":
      frame.update()
    if len(i) > 0:
      i = i.replace(">","").replace("\n","")
      tabs = i.split("^")
      protein = tabs[3]
      genome = tabs[0]
      tag = genome[:5].upper()
      #If new genome is reached, save data in pickle file
      if tag != lasttag and not fileinput.isfirstline():
        if GUI == "y":
          frame.update()
        #Load previous data if available
        if lasttag in passedtags:
          pickle_file = open("genecords" + os.sep + lasttag + ".pickle", "rb")
          previous_data = pickle.load(pickle_file)
          for key in previous_data:
            if key in genome_info:
              genome_info[key].extend(previous_data[key])
            else:
              genome_info[key] = previous_data[key]
        else:
          passedtags.append(lasttag)
        save_to_pickle(lasttag, genome_info)
        genome_info = {}
      if genome in genome_info:
        genome_info[genome].append(i)
      else:
        genome_info[genome] = [i]
      lasttag = tag
        
  #Repeat data saving
  if lasttag in passedtags:
    pickle_file = open("genecords" + os.sep + lasttag + ".pickle", "rb")
    previous_data = pickle.load(pickle_file)
    for key in previous_data:
      if key in genome_info:
        genome_info[key].extend(previous_data[key])
      else:
        genome_info[key] = previous_data[key]
  else:
    passedtags.append(lasttag)
  save_to_pickle(lasttag, genome_info)

  
  fileinput.close()

  #Sort dictionary by size
  #sortedgenecordskeylist = sortdictkeysbyvalues(genecordslistdict)

  #Archive directory as TAR file and remove original directory
  # try:
  #   if GUI == "y":
  #     frame.update()
  #   tar = tarfile.open(dbname + ".cords.tar", "w")
  #   tar.add("genecords")
  #   tar.close()
  # except:
  #   print("Could not create TAR file from genecords folder. Please create archive manually.")
  #   log("Could not create TAR file from genecords folder. Please create archive manually.", exit=True)
  # if GUI == "y":
  #   frame.update()
  # try:
  #   shutil.rmtree("genecords")
  # except:
  #   pass
