#!/usr/bin/env python
## Copyright (c) 2012 Marnix H. Medema
## Department of Microbial Physiology / Groningen Bioinformatics Centre
## University of Groningen
## License: GNU General Public License v3 or later
## A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

import sys
import os
from multiprocessing import freeze_support

global GUI
global CURRENTDIR
global MGBPATH
global APPDATA
global TEMP
global DBPATH
global USERDIR


#Find path to mgb files if run from another directory
pathfolders = os.environ['PATH'].split(os.pathsep)
pathfolders.append(os.getcwd())
CURRENTDIR = os.getcwd()
MGBPATH = ""
for folder in pathfolders:
  try:
    if "read_input_gui.py" in os.listdir(folder) and "guilib.py" in os.listdir(folder) and "empty.xhtml" in os.listdir(folder) and "multigeneblast.py" in os.listdir(folder) and "mgb_gui.py" in os.listdir(folder):
      MGBPATH = folder
      break
  except:
    pass
try:
  if  MGBPATH == "" and os.sep in sys.argv[0] and "read_input_gui.py" in os.listdir(sys.argv[0].rpartition(os.sep)[0]) and "guilib.py" in os.listdir(sys.argv[0].rpartition(os.sep)[0]):
    MGBPATH = sys.argv[0].rpartition(os.sep)[0]
    os.chdir(MGBPATH)
except:
  pass
if MGBPATH == "":
  print("Error: Please add the MultiGeneBlast installation directory to your $PATH environment variable before running the executable from another folder.")
  enter = eval(input())
  sys.exit(1)
 
#Find path to Application Data
if sys.platform == ('win32'):
  APPDATA = os.environ['ALLUSERSPROFILE'] + os.sep + 'Application Data'
elif sys.platform == ('darwin'):
  APPDATA = os.path.expanduser("~") + "/Library/Application Support"
else:
  try:
    if os.path.exists(os.getcwd() + os.sep + "multigeneblast_data"):
      APPDATA = os.getcwd() + os.sep + "multigeneblast_data"
    else:
      os.mkdir(os.getcwd() + os.sep + "multigeneblast_data")
      APPDATA = os.getcwd() + os.sep + "multigeneblast_data"
  except:
    try:
      if os.path.exists(os.environ['HOME'] + os.sep + "multigeneblast_data"):
        APPDATA = os.getcwd() + os.sep + "multigeneblast_data"
      else:
        os.mkdir(os.environ['HOME'] + os.sep + "multigeneblast_data")
        APPDATA = os.environ['HOME'] + os.sep + "multigeneblast_data"
    except:
      print("No permission to write to installation folder. Please change user or save somewhere else.")
      sys.exit()
if sys.platform == ('darwin') or sys.platform == ('win32'):
  try:
    os.mkdir(APPDATA + os.sep + 'MultiGeneBlast')
    APPDATA = APPDATA + os.sep + 'MultiGeneBlast'
  except:
    if os.path.exists(APPDATA + os.sep + 'MultiGeneBlast'):
      APPDATA = APPDATA + os.sep + 'MultiGeneBlast'
#Find path to temporary files
if sys.platform == ('win32'):
  TEMP = os.environ['TEMP']
elif sys.platform == ('darwin'):
  TEMP = os.environ['TMPDIR']
else:
  try:
    os.mkdir(os.environ['HOME'] + os.sep + ".mgbtemp")
    TEMP = os.environ['HOME'] + os.sep + ".mgbtemp"
  except:
    TEMP = APPDATA
#Set other environment variables
os.environ['EXEC'] = MGBPATH + os.sep + "exec"
os.environ['PATH'] = os.environ['EXEC'] + os.pathsep + os.environ['PATH']
if sys.platform == ('win32'):
  USERDIR = os.environ['USERPROFILE'] + os.sep + 'Documents'
else:
  USERDIR = os.environ['HOME']


from tkinter import *
import tkinter.filedialog
from tkinter.messagebox import askyesno, showerror, showinfo
from tkinter.ttk import Frame, Label, Scale, Style
import webbrowser
from guilib import *
from read_input_gui import *
from multigeneblast import *
import time
import urllib.request, urllib.error, urllib.parse
import tarfile
import traceback

def main_mgb(frame, opts):
  global outbox
  global TEMP
  global GUI
  GUI = "y"
  currentdir = os.getcwd()
  os.chdir(TEMP)
  outbox = MessageBox(frame, "Running MultiGeneBlast...")
  starttime = time.time()
  error = "n"
  conscious = "n"
    
  frame.update()
 
  #Step 1: parse options
  #Skipped in GUI:   #parse_options(sys.argv, opts)
  outbox.text_insert("Step 1/11: Time since start: " + str((time.time() - starttime)) + "\n")
  frame.update()
  
  #Step 2: Read GBK / EMBL file, select genes from requested region and output FASTA file
  while True:
    try:
      try:
        proteins, genomic_accnr, dnaseqlength, nucname, querytags, names, seqs, seqdict, arch_search = read_input_file(opts.infile, opts.startpos, opts.endpos, opts.ingenes, opts.gui, outbox, frame)
      except(SystemExit):
        error = "y"
        conscious = "y"
        break
      except:
        errormessage = traceback.format_exc()
        outbox.text_insert("Error (" + str(sys.exc_info()[0]).rpartition(".")[2][:-2] + ") while reading input file\n")
        error = "y"
        break
      outbox.text_insert("Step 2/11: Time since start: " + str((time.time() - starttime)) + "\n")
      frame.update()

      #Step 3: Run internal BLAST
      try:
        internalhomologygroupsdict, seqlengths = internal_blast(opts.minseqcov, opts.minpercid, names, proteins, seqdict, opts.nrcpus)
      except(SystemExit):
        error = "y"
        conscious = "y"
        break
      except:
        errormessage = traceback.format_exc()
        outbox.text_insert("Error (" + str(sys.exc_info()[0]).rpartition(".")[2][:-2] + ") while running internal Blast\n")
        error = "y"
        break
      outbox.text_insert("Step 3/11: Time since start: " + str((time.time() - starttime)) + "\n")
      frame.update()

      #Step 4: Run BLAST on genbank_mf database
      try:
        blastoutput = db_blast(names, seqs, opts.db, opts.nrcpus, opts.hitspergene, opts.dbtype)
      except(SystemExit):
        error = "y"
        conscious = "y"
        break
      except:
        errormessage = traceback.format_exc()
        outbox.text_insert("Error (" + str(sys.exc_info()[0]).rpartition(".")[2][:-2] + ") while running Blast\n")
        error = "y"
        break
      outbox.text_insert("Step 4/11: Time since start: " + str((time.time() - starttime)) + "\n")
      frame.update()

      #Step 5: Parse BLAST output
      try:
        blastdict, querylist = parse_blast(blastoutput, opts.minseqcov, opts.minpercid, seqlengths, seqdict, opts.db, opts.dbtype)
      except(SystemExit):
        error = "y"
        conscious = "y"
        break
      except:
        errormessage = traceback.format_exc()
        outbox.text_insert("Error (" + str(sys.exc_info()[0]).rpartition(".")[2][:-2] + ") while parsing Blast outputs\n")
        error = "y"
        break
      outbox.text_insert("Step 5/11: Time since start: " + str((time.time() - starttime)) + "\n")

      #Step 6: Load genomic databases into memory
      try:
        nucdescriptions, nucdict, proteininfo = load_databases(querylist, blastdict, opts.nrcpus, opts.db, opts.dbtype)
      except(SystemExit):
        error = "y"
        conscious = "y"
        break
      except:
        errormessage = traceback.format_exc()
        outbox.text_insert("Error (" + str(sys.exc_info()[0]).rpartition(".")[2][:-2] + ") while loading database:\n")
        error = "y"
        break
      outbox.text_insert("Step 6/11: Time since start: " + str((time.time() - starttime)) + "\n")

      #Step 7: Locate Blast hits in genomes
      try:
        blastdict, geneposdict, hitclusters, clusters, multiplehitlist = find_genomic_loci(blastdict, nucdict, proteininfo, opts.distancekb, querylist, nucdescriptions, dnaseqlength)
      except(SystemExit):
        error = "y"
        conscious = "y"
        break
      except:
        errormessage = traceback.format_exc()
        outbox.text_insert("Error (" + str(sys.exc_info()[0]).rpartition(".")[2][:-2] + ") while locating Blast hits\n")
        outbox.text_insert(errormessage)
        error = "y"
        break
      outbox.text_insert("Step 7/11: Time since start: " + str((time.time() - starttime)) + "\n")

      #Step 8: Score Blast output on all loci
      try:
        opts.pages = score_blast_output(hitclusters, querylist, blastdict, multiplehitlist, proteins, proteininfo, querytags, opts.infile, clusters, nucdescriptions, opts.pages, arch_search, opts.syntenyweight)
      except(SystemExit):
        error = "y"
        conscious = "y"
        break
      except:
        errormessage = traceback.format_exc()
        outbox.text_insert("Error (" + str(sys.exc_info()[0]).rpartition(".")[2][:-2] + ") while scoring Blast outputs\n")
        error = "y"
        break
      outbox.text_insert("Step 8/11: Time since start: " + str((time.time() - starttime)) + "\n")

      #Output. From here, iterate for every page
      for page in [pagenr + 1 for pagenr in range(opts.pages)]:

        #Step 9: Write MultiGeneBlast SVGs
        try:
          queryclusterdata, colorschemedict, clusterblastpositiondata, blastdetails, mgb_scores = write_svgs(page, opts.screenwidth, internalhomologygroupsdict, arch_search)
        except(SystemExit):
          error = "y"
          conscious = "y"
          break
        except:
          errormessage = traceback.format_exc()
          outbox.text_insert("Error (" + str(sys.exc_info()[0]).rpartition(".")[2][:-2] + ") while writing SVGs\n")
          error = "y"
          break
        outbox.text_insert("Step 9/11, page " + str(page) + ": Time since start: " + str((time.time() - starttime)) + "\n")

        #Step 10: Create muscle alignments
        try:
          musclegroups = align_muscle(opts.muscle, colorschemedict, seqdict)
        except(SystemExit):
          error = "y"
          conscious = "y"
          break
        except:
          errormessage = traceback.format_exc()
          outbox.text_insert("Error (" + str(sys.exc_info()[0]).rpartition(".")[2][:-2] + ") while creating Muscle alignments\n")
          error = "y"
          break
        outbox.text_insert("Step 10/11, page " + str(page) + ": Time since start: " + str((time.time() - starttime)) + "\n")

        #Step 11: Create XHTML output file
        try:
          create_xhtml_file(queryclusterdata, clusters, clusterblastpositiondata, nucname, page, opts.pages, opts.screenwidth, blastdetails, mgb_scores, musclegroups, colorschemedict, opts.dbtype)
        except(SystemExit):
          error = "y"
          conscious = "y"
          break
        except:
          errormessage = traceback.format_exc()
          outbox.text_insert("Error while creating XHTML file\n")
          error = "y"
          break
        outbox.text_insert("Step 11/11, page " + str(page) + ": Time since start: " + str((time.time() - starttime)) + "\n")
        
    except(MemoryError):
      errormessage = traceback.format_exc()
      outbox.text_insert("Error: not enough memory to complete the operation.\n")
      error = "y"
      break
    except(SystemExit):
      conscious = "y"
      break
    except:
      errormessage = traceback.format_exc()
      outbox.text_insert("Unexpected error while running MultiGeneBlast\n")
      error = "y"
      break

    #Move all files to specified output folder
    try:
      move_outputfiles(opts.outputfolder, opts.pages)
    except:
      errormessage = traceback.format_exc()
      error = "y"
      outbox.text_insert("Error in copying files to specified folder\n")
      break

    #Open display page of first XHTML file
    if sys.platform == ('win32'):
      os.startfile(opts.outputfolder + os.sep + 'displaypage1.xhtml')
    else:
      try:
        firefox = webbrowser.get('firefox')
        firefox.open_new_tab("file://" + opts.outputfolder + os.sep + 'displaypage1.xhtml')
      except:
        try:
          safari = webbrowser.get('safari')
          safari.open_new_tab("file://" + opts.outputfolder + os.sep + 'displaypage1.xhtml')
        except:
          try:
            chrome = webbrowser.get('/usr/bin/google-chrome %s')
            chrome.open_new_tab("file://" + opts.outputfolder + os.sep + 'displaypage1.xhtml')
          except:
            try:
              webbrowser.open_new_tab("file://" + opts.outputfolder + os.sep + 'displaypage1.xhtml')
            except:
              pass
    outbox.text_insert("The output can be accessed by opening " + opts.outputfolder + os.sep + 'displaypage1.xhtml with a web browser\n')
    break

  
  #Add finishing message and add OK button
  if error == "n":
    outbox.text_insert("MultiGeneBlast successfully finished in " + str((time.time() - starttime)) + " seconds.\n")
  elif conscious == "y":
    pass
  else:
    outbox.text_insert("MultiGeneBlast experienced a problem. Please click the button below to send an error report, so we can solve the underlying problem.\n")
    outbox.change_errormessage(errormessage)
    outbox.add_error_button()
  outbox.add_ok_button()
  os.chdir(currentdir)



global runmgb
def runmgb(frame, infile, cstart, cend, genes, outfolder, database, cores, seqcov, pid, distancekb, pages, muscle, searchtype, hitspergene, syntenyscore):
    #Check for valid input file
    if database == "<No database selected>":
      showerror('Input error','No database selected.')
      return
    try:
      open(infile)
    except:
      showerror('Input error','No input file given.')
      return
      
    if outfolder == "<Select output folder>" or outfolder == "":
      showerror('Input error',"Please supply a folder name to store the output in.")
      return
        
    if muscle == 1:
      answer = askyesno('Warning', 'Generating MUSCLE alignments for all queries with all homologs may be time-consuming. Continue anyway?')
      if not answer:
        return

    #Check for correct start/end position or genes input
    if searchtype != "architecture":
      if not (cstart.isdigit() and cend.isdigit()) and not (cstart == "" and cend == "" and genes != ""):
        showerror('Input error',"No correct start or end position given.")
        return
      if genes == "<Select genes>" and cstart.isdigit() and cend.isdigit():
        genes = ""
        pass
      else:
        genes = checkgenes(genes)
        if genes == "<ERROR>":
          showerror('Input error',"Invalid characters found in gene names.")
          return
    
    #Execute MGB command #Can be used for calling multigeneblast.py externally
    #if searchtype != "architecture":
    #  if (cstart == "" and cend == "" and genes != ""):
    #    command = 'python multigeneblast.py -in "' + infile + '" -genes ' + genes + ' -db "' + database + '"'
    #  else:
    #    command = 'python multigeneblast.py -in "' + infile + '" -from ' + cstart + " -to " + cend + ' -db "' + database + '"'
    #else:
    #  command = 'python multigeneblast.py -in "' + infile + '"' + ' -db "' + database + '"'
    #
    #command = command + " -out " + outfolder + " -cores " + str(cores) + " -minseqcov " + str(seqcov) + " -minpercid " + str(pid) + " -distancekb " + str(distancekb) + " -outpages " + str(pages / 50)
    #if muscle == 1:
    #  command = command + " -muscle y"
    #else:
    #  command = command + " -muscle n"
    #print command
    #os.system(command)
    
    #Set options for MGB run
    opts = Options()
    opts.outputfolder = outfolder
    #infile = parse_absolute_paths(infile)
    opts.infile = infile
    global DBPATH
    DBPATH = database.rpartition("/")[0]
    if sys.platform == ('win32'):
      DBPATH = DBPATH.replace("/", os.sep)
    os.environ['BLASTDB'] = DBPATH
    database = parse_absolute_paths(database)
    if ".pal" in database:
      opts.dbtype = "prot"
    else:
      opts.dbtype = "nucl"
    database = database.partition(".pal")[0].partition(".nal")[0]
    opts.db = database
    opts.nrcpus = cores
    opts.minseqcov = int(seqcov)
    opts.syntenyweight = float(syntenyscore)
    opts.minpercid = int(pid)
    opts.screenwidth = 1024
    opts.distancekb = int(distancekb) * 1000
    if int(muscle) == 1:
      opts.muscle = "y"
    else:
      opts.muscle = "n"
    opts.startpos = int(cstart)
    opts.endpos = int(cend)
    opts.ingenes = [gene for gene in genes.split(";") if gene != ""]
    opts.hitspergene = int(hitspergene)
    if len(opts.ingenes) > 0 and ";" in genes:
      opts.startpos = "N/A"
      opts.endpos = "N/A"
    else:
      opts.ingenes = "N/A"
    opts.pages = pages / 50
    opts.gui = "y"
    
    main_mgb(frame, opts)
    

def file_open():
    global APPDATA
    if SearchPrefs.searchtype.get() == "homology":
      location = tkinter.filedialog.askopenfilename(filetypes=(("GenBank files",('*.gbk','*.gb','*.genbank')),("EMBL files",('*.embl','*.emb'))), initialdir=APPDATA)
      infile.set(location)
      genes, dnaseqlength = read_input_file_gui(infile.get())
      if dnaseqlength != 0:
        SearchPrefs.geneslist = genes
        SearchPrefs.dnaseqlength = dnaseqlength
        SearchPrefs.update_starts_ends(frame)
      else:
        infile.set("<No file selected>")
    else:
      location = tkinter.filedialog.askopenfilename(filetypes=(("FASTA files",('*.fasta','*.fas','*.fa')),("All file types",('*.*'))))
      if checkfasta(location) == 1:
        infile.set(location)
    #frame.update_idletasks()

def db_open():
    global APPDATA
    location = tkinter.filedialog.askopenfilename(initialdir = APPDATA, filetypes=(("MGB database files",('*.pal','*.nal')),("MGB database files",())))
    if len(location) > 0:
      folder = location.rpartition("/")[0]
      basename = location.rpartition("/")[2].rpartition(".")[0]
    if len(location) > 0 and len(folder) > 0 and len(basename) > 0:
      #if folder != os.getcwd() and folder.replace("/","\\") != os.getcwd():
      #  showerror('Database error',"Please select a database within the MultiGeneBlast working directory.")
      #  db_open()
      if ".pal" in location:
        letter = "p"
      else:
        letter = "n"
      if not (basename + "." + letter + "sq" in os.listdir(folder) or basename + ".00." + letter + "sq" in os.listdir(folder)):
        showerror('Database error',"No *.psq file found for this database.")
        db_open()
      elif not (basename + "." + letter + "hr" in os.listdir(folder) or basename + ".00." + letter + "hr" in os.listdir(folder)):
        showerror('Database error',"No *.phr file found for this database.")
        db_open()
      elif not (basename + "." + letter + "in" in os.listdir(folder) or basename + ".00." + letter + "in" in os.listdir(folder)):
        showerror('Database error',"No *.pin file found for this database.")
        db_open()
      elif not (basename + "_all_descrs.txt" in os.listdir(folder)):
        showerror('Database error',"No descriptions file found for this database.")
        db_open()
      elif ".pal" in location:
        if not (basename + ".cords.tar" in os.listdir(folder)):
          showerror('Database error',"No coordinates TAR file found for this database.")
          db_open()
        elif not (basename + ".pinfo.tar" in os.listdir(folder)):
          showerror('Database error',"No protein info TAR file found for this database.")
          db_open()
      database.set(location)


class SearchPreferences:

  #Response to searchtype radio buttons
  global defineSearchType
  def defineSearchType(self, frame, searchtype):
    if searchtype.get() == 1:
      self.searchtype.set("homology")
      infile.set("<No file selected>")
      if len(list(self.text2.grid_info().keys())) == 0:
        self.text2.grid(row=6,column=0)
        self.text3.grid(row=7,column=0)
        self.text4.grid(row=8,column=0)
        self.cstart = ScaleBar(self.master, self, [6,1], 0, self.dnaseqlength, 0, resetgenes="y")
        self.cend = ScaleBar(self.master, self, [7,1], 0, self.dnaseqlength, self.dnaseqlength, resetgenes="y")
        self.selectedgenes = GeneSelectionFrame(self.master, self, [8,1], self.geneslist)
    elif searchtype.get() == 2:
      self.searchtype.set("architecture")
      infile.set("<No file selected>")
      if len(list(self.text2.grid_info().keys())) > 0:
        self.text2.grid_forget()
        self.text3.grid_forget()
        self.text4.grid_forget()
        self.cstart.grid_remove()
        self.cend.grid_remove()
        self.selectedgenes.grid_remove()
        self.dnaseqlength = 0
        self.geneslist = []

  #Search type radio buttons
  global searchtype
  def searchtype(self, frame):
    space1 = Label(frame, text="\n")
    space1.grid(row=3,column=0, columnspan=2)
    searchtype = IntVar()
    homradio = Radiobutton(frame, text="Homology search", variable=searchtype, command=lambda:defineSearchType(self, frame, searchtype), value=1)
    homradio.grid(row=4,column=0)
    archradio = Radiobutton(frame, text="Architecture search", variable=searchtype, command=lambda:defineSearchType(self, frame, searchtype), value=2)
    archradio.grid(row=4,column=1)
    homradio.select()
    space2 = Label(frame, text="\n")
    space2.grid(row=5,column=0, columnspan=2)

  #Return functions
  def getgenes(self):
    return self.genes.get()
  def getfolder(self):
    return self.outputfolder.get()
  def getsearchtype(self):
    return self.searchtype.get()

  def update_starts_ends(self, master):
    if self.searchtype.get() != "architecture":
      self.cstart.grid_remove()
      self.cend.grid_remove()
      self.selectedgenes.grid_remove()
      self.cstart = ScaleBar(master, self, [6,1], 0, self.dnaseqlength, 0, resetgenes="y")
      self.cend = ScaleBar(master, self, [7,1], 0, self.dnaseqlength, self.dnaseqlength, resetgenes="y")
      self.selectedgenes = GeneSelectionFrame(master, self, [8,1], self.geneslist)

  #Class constructor with class widgets
  def __init__(self,master):

    searchtype(self, master)
    self.master = master
    self.searchtype = StringVar()
    self.searchtype.set("homology")
    self.geneslist = []
    self.dnaseqlength = 0
    #space = Label(master, text="\n")
    #space.pack()

    #Start nt selection
    self.text2 = Label(master, text="Nucleotide start position of fragment to search:")
    self.text2.grid(row=6,column=0, sticky=N)
    self.cstart = ScaleBar(master, self, [6,1], 0, self.dnaseqlength, 0, resetgenes="y")

    #End nt selection
    self.text3 = Label(master, text="Nucleotide end position of fragment to search:")
    self.text3.grid(row=7,column=0)
    self.cend = ScaleBar(master, self, [7,1], 0, self.dnaseqlength, self.dnaseqlength, resetgenes="y")

    #Genes selection
    self.text4 = Label(master, text="\nOr genes to search (locus tags / accession numbers):")
    self.text4.grid(row=8,column=0, sticky=N)
    self.selectedgenes = GeneSelectionFrame(master, self, [8,1], self.geneslist)
    
    #Output folder name
    text5 = Label(master, text="\nOutput folder name:")
    text5.grid(row=9,column=0, sticky=N)
    self.outputfolder = OutputFolderSelectionFrame(master, [9,1], USERDIR)
    
    #Number of cores to use
    text6 = Label(master, text="\nNumber of CPU cores to be used:")
    text6.grid(row=10,column=0, sticky=N)
    self.cores = ScaleBar(master, self, [10,1], 1, determine_cpu_nr("all"), determine_cpu_nr("all"))

    #Minimal sequence coverage of Blast hit
    text7 = Label(master, text="\nNumber of Blast hits per gene to be mapped:")
    text7.grid(row=11,column=0, sticky=N)
    self.hitspergene = ScaleBar(master, self, [11,1], 50, 1000, 250)

    #Weight of synteny conservation score
    text8 = Label(master, text="\nWeight of synteny conservation in hit sorting:")
    text8.grid(row=12,column=0, sticky=N)
    self.syntenyweight = ScaleBar(master, self, [12,1], 0.0, 2.0, 0.5, input_type="float")

    #Minimal sequence coverage of Blast hit
    text9 = Label(master, text="\nMinimal sequence coverage of BLAST hits:")
    text9.grid(row=13,column=0, sticky=N)
    self.seqcov = ScaleBar(master, self, [13,1], 0, 100, 25)
    
    #Minimal % ID of Blast hit
    text10 = Label(master, text="\nMinimal % identity of BLAST hits:")
    text10.grid(row=14,column=0, sticky=N)
    self.pid = ScaleBar(master, self, [14,1], 0, 100, 30)
    
    #Maximum allowed distance between hits in locus
    text11 = Label(master, text="\nMaximum distance between genes in locus (kb):")
    text11.grid(row=15,column=0, sticky=N)
    self.distancekb = ScaleBar(master, self, [15,1], 1, 100, 20)
    
    #Number of hit loci to show
    text12 = Label(master, text="\nNumber of hit loci to display:")
    text12.grid(row=16,column=0, sticky=N)
    self.pages = SpinBox(master, [16,1], 0, 500, 50, 250)
        
    #Generate Muscle alignments
    text13 = Label(master, text="\nMuscle alignment of homologs with queries:")
    text13.grid(row=17,column=0, sticky=N)
    self.muscle = CheckBox(master, [17,1], "")
    
    #Empty space
    space3 = Label(master, text="   ")
    space3.grid(row=18,column=0, columnspan=2)
    #space = Label(master, text="\n")
    #space.pack()

def file_download(frame):
  global APPDATA
  GenBankFileDownload(frame, APPDATA)

def cancel_extraction(tar, toplevel):
  tar.close()
  toplevel.destroy()

def db_download():
  file_url = 'http://gbic.webhosting.rug.nl/genbank_mf.tar.gz'
  #file_url = 'https://downloads.sourceforge.net/project/multigeneblast/db/genbank_mf_test.rar?use_mirror=switch'
  answer = askyesno('Confirmation', 'Database downloading will take a while, and the database will occuppy ~15 GB disk space.\n Are you sure?')
  if not answer:
    return
  global APPDATA
  currentdir = os.getcwd()
  os.chdir(APPDATA)
  #Check if database is not already present
  if "genbank_mf.pal" in os.listdir("."):
    answer = askyesno('File already present', 'Database appears to be present already. Overwrite?')
    if not answer:
      return
  #Check if sufficient disk space is available
  default_file_size = 357266739200
  dbfilesize = int(get_remote_filesize(file_url))
  if int(get_free_space(".")) < dbfilesize:
    showerror('Error', 'Insufficient disk space available.\n15GB needed.')
    return
  #Open 'loading...' message
  loading = Toplevel(frame, height=200, width=400)
  loading.title("Download progress")
  message = "Downloading: " + str(0) + "/ " + str(int(dbfilesize / 1024)) + " kb\nPlease wait..."
  a = Label(loading, text=message)
  a.grid(row=1,column=1, padx=50, pady=50)
  loading.bind("<Escape>", lambda e: "return")
  frame.update()
  #Download the MGB GenBank database
  try:
    req = urllib.request.urlopen(file_url)
  except:
    showerror('Error', 'File not found on server.\nPlease check your internet connection.')
    return
  loading.protocol('WM_DELETE_WINDOW', lambda: cancel_download(req, loading))
  CHUNK = 128 * 1024
  x = 0
  with open("genbank_mf.tar.gz", 'wb') as fp:
    while True:
      message = "Downloading: " + str(x) + "/ " + str(int(dbfilesize / 1024)) + " kb\nPlease wait..."
      try:
        a = Label(loading, text=message)
        a.grid_remove()
        frame.update()
        a.grid(row=1,column=1, padx=20, pady=20)
        frame.update()
        chunk = req.read(CHUNK)
        if not chunk:
          break
        fp.write(chunk)
        x += 128
        if x > int(dbfilesize / 1024):
          x = int(dbfilesize / 1024)
      except:
        showerror('Error', 'Download interrupted.')
        return
  #Report download success
  loading.destroy()
  showinfo("Download finished", "Download completed successfully.")
  #Extract the database TAR/GZ file
  extracting = Toplevel(frame, height=200, width=400)
  extracting.title("Database extraction")
  message = "Extracting database.\nPlease wait..."
  a = Label(extracting, text=message)
  a.grid(row=1,column=1, padx=50, pady=50)
  extracting.bind("<Escape>", lambda e: "return")
  frame.update()
  try:
    tar = tarfile.open("genbank_mf.tar.gz")
    extracting.protocol('WM_DELETE_WINDOW', lambda: cancel_download(tar, extracting))
    frame.update()
    tar.extractall()
    tar.close()
  except:
    showerror("Error","Error extracting database. Please try to extract it manually.")
    extracting.destroy()
    os.chdir(currentdir)
    return
  extracting.destroy()
  os.chdir(currentdir)
  showinfo("Extraction finished", "Database extraction finished.\nYou can now use this database by selecting 'genbank_mf.pal' under 'Select database' in the 'File' menu.")

def makedb_file(frame):
  global APPDATA
  MakeDatabase(frame, APPDATA)

def makedb_ncbi(frame):
  global APPDATA
  MakeOnlineDatabase(frame, APPDATA)

def makedb_gb():
  global APPDATA
  MakeGenBankDatabase(frame, APPDATA)

def menu(root):
  menu = Menu(root)
  root.config(menu=menu)

  filemenu = Menu(menu)
  menu.add_cascade(label="File", menu=filemenu)
  filemenu.add_command(label="Open input file", command=file_open)
  filemenu.add_command(label="Select database", command=db_open)
  filemenu.add_command(label="Exit", command=root.quit)

  downloadmenu = Menu(menu)
  menu.add_cascade(label="Download", menu=downloadmenu)
  downloadmenu.add_command(label="Download GenBank entry", command=lambda: file_download(frame))
  downloadmenu.add_command(label="Download MGB GenBank database", command=db_download)

  dbmenu = Menu(menu)
  menu.add_cascade(label="Database", menu=dbmenu)
  dbmenu.add_command(label="Create database from files", command=lambda: makedb_file(frame))
  dbmenu.add_command(label="Create database from online GenBank entries", command=lambda: makedb_ncbi(frame))
  dbmenu.add_command(label="Create database from GenBank subdivisions", command=makedb_gb)

  helpmenu = Menu(menu)
  menu.add_cascade(label="Help", menu=helpmenu)
  helpmenu.add_command(label="About...", command=about)

def maingui():

  os.environ['PYTHON'] = os.path.dirname(os.path.abspath(__file__)) + os.sep + "python"
  os.environ['EXEC'] = os.path.dirname(os.path.abspath(__file__)) + os.sep + "exec"
  os.environ['PATH'] = os.environ['EXEC'] + os.pathsep + os.environ['PYTHON'] + os.pathsep + os.environ['PATH']
  
  root = Tk()
  if sys.platform == ('win32'):
    try:
      root.iconbitmap(default='mgb.ico')
    except:
      pass
  root.title('MultiGeneBlast')
  root.geometry("%dx%d%+d%+d" % (850, 750, 0, 0))

  global frame
  frame = Frame()
  frame.pack()
  frame.grid_columnconfigure(0, minsize=300)
  frame.grid_columnconfigure(1, minsize=300)
  frame.grid_rowconfigure(8, minsize=32)
  frame.grid_rowconfigure(9, minsize=32)
  frame.grid_rowconfigure(10, minsize=34)
  frame.grid_rowconfigure(11, minsize=34)
  frame.grid_rowconfigure(12, minsize=34)
  frame.grid_rowconfigure(13, minsize=34)
  frame.grid_rowconfigure(14, minsize=34)
  frame.grid_rowconfigure(15, minsize=34)
  frame.grid_rowconfigure(16, minsize=34)

  htext1 = Label(frame, text="MultiGeneBlast input\n")
  htext1.grid(row=0,column=0)

  #Database selection
  global database
  database = StringVar()
  if "genbank_mf.pal" in os.listdir("."):
    database.set((os.getcwd() + os.sep + "genbank_mf").replace("\\","/"))
  else:
    database.set("<No database selected>")
  text0A = Label(frame, text = "Current database:")
  text0A.grid(row=1,column=0)
  text0B = Label(frame, textvariable = database)
  text0B.grid(row=1,column=1)

  #Input file selection
  global infile
  infile = StringVar()
  infile.set("<No file selected>")
  text1A = Label(frame, text = "Current input file: ")
  text1A.grid(row=2,column=0)
  text1B = Label(frame, textvariable = infile)
  text1B.grid(row=2,column=1)

  #Search preferences
  global SearchPrefs
  SearchPrefs = SearchPreferences(frame)

  #Run button to start analysis
  button = Button(frame, text="Run MultiGeneBlast", command=lambda : runmgb(frame, infile.get(), str(SearchPrefs.cstart.getval()), str(SearchPrefs.cend.getval()), SearchPrefs.selectedgenes.getval(), SearchPrefs.outputfolder.getfolder(), database.get(), SearchPrefs.cores.getval(), SearchPrefs.seqcov.getval(), SearchPrefs.pid.getval(), SearchPrefs.distancekb.getval(), SearchPrefs.pages.getval(), SearchPrefs.muscle.getval(), SearchPrefs.getsearchtype(), SearchPrefs.hitspergene.getval(), SearchPrefs.syntenyweight.getval()))
  button.grid(row=20,column=0)

  # create a menu
  menu(root)
  
if __name__ == '__main__':
  freeze_support()
  maingui()
  mainloop()
