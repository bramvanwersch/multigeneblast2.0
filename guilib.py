#!/usr/bin/env python
## Copyright (c) 2012 Marnix H. Medema
## Department of Microbial Physiology / Groningen Bioinformatics Centre
## University of Groningen
## License: GNU General Public License v3 or later
## A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

import multiprocessing
import webbrowser
import tkinter.filedialog
from tkinter.messagebox import askyesno, showerror, showwarning, showinfo
from tkinter.ttk import Frame, Label, Scale, Style
from tkinter import Tk, BOTH, IntVar, Checkbutton, Spinbox, Listbox, StringVar, Entry, Button, Toplevel, TOP, RIGHT, LEFT, BOTTOM, Y, END, EXTENDED, Scrollbar, Text, NORMAL, INSERT, DISABLED, SINGLE, S, N, TclError
import os
import sys
import shutil
import platform
import ctypes
import urllib.request, urllib.error, urllib.parse
from dblib.parse_gbk import parse_gbk_embl, nucparse_gbk_embl_fasta
from dblib.utils import fasta_names, make_blast_db, clean_up, get_accession
from dblib.generate_genecords_tar import generate_genecords_tar
from dblib.generate_proteininfo_tar import generate_proteininfo_tar

from utilities import ILLEGAL_CHARACTERS
# import makegbdb
# import makegbndb

#Find path to mgb files if run from another directory
# pathfolders = os.environ['PATH'].split(os.pathsep)
# pathfolders.append(os.getcwd())
# CURRENTDIR = os.getcwd()
# MGBPATH = ""
# for folder in pathfolders:
#   try:
#     if "read_input_gui.py" in os.listdir(folder) and "guilib.py" in os.listdir(folder) and "empty.xhtml" in os.listdir(folder) and "multigeneblast.py" in os.listdir(folder) and "mgb_gui.py" in os.listdir(folder):
#       MGBPATH = folder
#       break
#   except:
#     pass
# try:
#   if  MGBPATH == "" and os.sep in sys.argv[0] and "read_input_gui.py" in os.listdir(sys.argv[0].rpartition(os.sep)[0]) and "guilib.py" in os.listdir(sys.argv[0].rpartition(os.sep)[0]):
#     MGBPATH = sys.argv[0].rpartition(os.sep)[0]
#     os.chdir(MGBPATH)
# except:
#   pass
# if MGBPATH == "":
#   print("Error: Please add the MultiGeneBlast installation directory to your $PATH environment variable before running the executable from another folder.")
#   enter = eval(input())
#   sys.exit(1)
#
# #Find path to Application Data
# if sys.platform == ('win32'):
#   APPDATA = os.environ['ALLUSERSPROFILE'] + os.sep + 'Application Data'
# elif sys.platform == ('darwin'):
#   APPDATA = os.path.expanduser("~") + "/Library/Application Support"
# else:
#   try:
#     if os.path.exists(os.getcwd() + os.sep + "multigeneblast_data"):
#       APPDATA = os.getcwd() + os.sep + "multigeneblast_data"
#     else:
#       os.mkdir(os.getcwd() + os.sep + "multigeneblast_data")
#       APPDATA = os.getcwd() + os.sep + "multigeneblast_data"
#   except:
#     try:
#       if os.path.exists(os.environ['HOME'] + os.sep + "multigeneblast_data"):
#         APPDATA = os.getcwd() + os.sep + "multigeneblast_data"
#       else:
#         os.mkdir(os.environ['HOME'] + os.sep + "multigeneblast_data")
#         APPDATA = os.environ['HOME'] + os.sep + "multigeneblast_data"
#     except:
#       print("No permission to write to installation folder. Please change user or save somewhere else.")
#       sys.exit()
# if sys.platform == ('darwin') or sys.platform == ('win32'):
#   try:
#     os.mkdir(APPDATA + os.sep + 'MultiGeneBlast')
#     APPDATA = APPDATA + os.sep + 'MultiGeneBlast'
#   except:
#     if os.path.exists(APPDATA + os.sep + 'MultiGeneBlast'):
#       APPDATA = APPDATA + os.sep + 'MultiGeneBlast'
# #Find path to temporary files
# if sys.platform == ('win32'):
#   TEMP = os.environ['TEMP']
# elif sys.platform == ('darwin'):
#   TEMP = os.environ['TMPDIR']
# else:
#   try:
#     os.mkdir(os.environ['HOME'] + os.sep + ".mgbtemp")
#     TEMP = os.environ['HOME'] + os.sep + ".mgbtemp"
#   except:
#     TEMP = APPDATA
# #Set other environment variables
# os.environ['EXEC'] = MGBPATH + os.pathsep + "exec"
# os.environ['PATH'] = os.environ['EXEC'] + os.pathsep + os.environ['PATH']
# if sys.platform == ('win32'):
#   USERDIR = os.environ['USERPROFILE'] + os.sep + 'Documents'
# else:
#   USERDIR = os.environ['HOME']

def cancel_download(req, toplevel):
  req.close()
  toplevel.destroy()

def get_remote_filesize(url):
  try:
    usock = urllib.request.urlopen(url)
    dbfilesize = usock.info().get('Content-Length')
    if dbfilesize is None:
        dbfilesize = 0
  except:
    dbfilesize = 0
  dbfilesize = float(int(dbfilesize)) # db file size in bytes
  return dbfilesize

def get_free_space(folder):
    """ Return folder/drive free space (in bytes)
    """
    if platform.system() == 'Windows':
        free_bytes = ctypes.c_ulonglong(0)
        ctypes.windll.kernel32.GetDiskFreeSpaceExW(ctypes.c_wchar_p(folder), None, None, ctypes.pointer(free_bytes))
        return free_bytes.value
    else:
        return os.statvfs(folder).f_bfree * os.statvfs(folder).f_frsize

def determine_cpu_nr(cores):
    #Determine number of CPUs used
    if cores == "all":
        try:
            nrcpus = multiprocessing.cpu_count()
        except(IOError,OSError,NotImplementedError):
            nrcpus = 1
    else:
        try:
            nrcpus = multiprocessing.cpu_count()
        except(IOError,OSError,NotImplementedError):
            nrcpus = 1
        if cores < nrcpus:
            nrcpus = cores
    return nrcpus
    
global checkgenes
def checkgenes(genes):
    genes = genes.replace(", ",";")
    genes = genes.replace(",",";")
    genes = genes.replace(" ",";")
    genes = genes.replace("'","")
    genes = genes.replace('"','')
    forbiddencharacters = ["'",'"','=',':','[',']','>','<','|','\\',"/",'*','-',',','?',')','(','^','#','!','`','~','+','{','}','@','$','%','&']
    for char in genes:
      if char in forbiddencharacters:
        return "<ERROR>"
    return genes
    
global checkfoldername
def checkfoldername(foldername):
    forbiddencharacters = ["'",'"','=',':','[',']','>','<','|','\\',"/",'*','-',',','?',')','(','^','#','!','`','~','+','{','}','@','$','%','&']
    for char in foldername:
      if char in forbiddencharacters:
        return "<ERROR>"
    return foldername

def about():
    webbrowser.open('http://multigeneblast.sourceforge.net/')
    print("Navigating to MultiGeneBlast web page on Sourceforge.net")

def validname(string):
    forbiddencharacters = ["'",'"','=',':',';','-','[',']','>','<','|','\\',"/",'*',',','?',')','(','^','#','!','`','~','+','{','}','@','$','%','&']
    for char in forbiddencharacters:
      if char in string:
        return False
    return True

class ScaleBar(Frame):
  
    def __init__(self, parent, SpObj, positions, minimum, maximum, default="0", resetgenes="n", input_type="int"):
        Frame.__init__(self, parent)   
         
        self.parent = parent
        self.SpObj = SpObj
        self.positions = positions
        self.minimum = minimum
        self.maximum = maximum
        self.default = default
        self.resetgenes = resetgenes
        self.input_type = input_type
        self.initUI()
        
    def initUI(self):
      
        #self.parent.title("Scale")
        #self.style = Style()
        #self.style.theme_use("default")
        
        #self.pack(fill=BOTH, expand=1)
        
        self.grid(row=self.positions[0],column=self.positions[1], sticky=S)
        
        scale = Scale(self, from_=self.minimum, to=self.maximum, 
            command=self.onScale, length=200)
        scale.grid(row=0,column=1)
        self.scale = scale
        self.var = StringVar()
        self.var.set(self.default)
        self.lastvar = StringVar()
        self.lastvar.set(self.var.get())
        scale.set(self.default)
        self.entry = Entry(self, textvariable=self.var, width=8)
        self.entry.bind("<FocusOut>", self.OnValidate)
        self.entry.bind("<Return>", self.OnValidate)
        self.entry.grid(row=0,column=0)
        
    def OnValidate(self, val):
        if "-" in self.var.get():
          self.var.set(self.minimum)
        if not str(self.var.get()).isdigit():
          self.var.set(self.lastvar.get())
        if int(self.var.get()) <= int(self.minimum):
          self.var.set(self.minimum)
        if int(self.maximum) - int(self.var.get()) < 0:
          self.var.set(self.maximum)
        self.onScale(self.var.get())
        self.scale.set(self.var.get())
        self.lastvar.set(str(self.var.get()))

    def onScale(self, val):
        if self.input_type == "int":
          v = int(float(val))
        else:
          v = float("".join(str(val).partition(".")[0:2]) + str(val).partition(".")[2][:2])
        self.var.set(str(v))
        if hasattr(self.SpObj, 'selectedgenes') and self.resetgenes == "y":
          self.SpObj.selectedgenes.clear_selection()
        
    def setval(self, val):
        self.var.set(val)
        self.scale.set(val)

    def getval(self):
        return self.var.get()
        
class CheckBox(Frame):
  
    def __init__(self, parent, positions, description):
        Frame.__init__(self, parent)   
         
        self.parent = parent
        self.positions = positions
        self.description = description
        self.initUI()
        self.var.set(0)
        
    def initUI(self):
      
        #self.parent.title("Checkbutton")

        self.grid(row=self.positions[0],column=self.positions[1], sticky=S)
        self.var = IntVar()
        
        cb = Checkbutton(self, text=self.description,
            variable=self.var)
        cb.select()
        cb.grid(row=0,column=0)
            
    def getval(self):
        return self.var.get()
        
        
class SpinBox(Frame):
  
    def __init__(self, parent, positions, minimum, maximum, incrval, default=0):
        Frame.__init__(self, parent)   
         
        self.parent = parent
        self.positions = positions
        self.minimum = minimum
        self.maximum = maximum
        self.incrval = incrval
        self.default = default
        self.initUI()
        
    def initUI(self):
      
        #self.parent.title("Checkbutton")

        self.grid(row=self.positions[0],column=self.positions[1], sticky=S)
        self.var = IntVar()
        self.var.set(self.default)
        
        cb = Spinbox(self, from_=self.minimum, to=self.maximum, increment=self.incrval,
            textvariable=self.var)
        cb.grid(row=0,column=0)
            
    def getval(self):
        return self.var.get()
        
class ListBoxChoice(object):
    def __init__(self, master=None, title=None, message=None, list=[]):
        self.master = master
        self.value = None
        self.list = list[:]
        
        self.modalPane = Toplevel(self.master)

        self.modalPane.transient(self.master)
        self.modalPane.grab_set()

        self.modalPane.bind("<Return>", self._choose)
        self.modalPane.bind("<Escape>", self._cancel)

        if title:
            self.modalPane.title(title)

        if message:
            Label(self.modalPane, text=message).pack(padx=5, pady=5)

        listFrame = Frame(self.modalPane)
        listFrame.pack(side=TOP, padx=5, pady=5)
        
        scrollBar = Scrollbar(listFrame)
        scrollBar.pack(side=RIGHT, fill=Y)
        self.listBox = Listbox(listFrame, selectmode=EXTENDED)
        self.listBox.pack(side=LEFT, fill=Y)
        scrollBar.config(command=self.listBox.yview)
        self.listBox.config(yscrollcommand=scrollBar.set)
        self.list.sort()
        for item in self.list:
            self.listBox.insert(END, item)

        buttonFrame = Frame(self.modalPane)
        buttonFrame.pack(side=BOTTOM)

        chooseButton = Button(buttonFrame, text="Choose", command=self._choose)
        chooseButton.pack()

        cancelButton = Button(buttonFrame, text="Cancel", command=self._cancel)
        cancelButton.pack(side=RIGHT)

    def _choose(self, event=None):
        try:
            if len(self.listBox.curselection()) == 1:
              Selected = self.listBox.curselection()[0]
              self.value = self.list[int(Selected)]
            else:
              Selected = self.listBox.curselection()
              self.value = ";".join([self.list[int(idx)] for idx in Selected])
        except IndexError:
            self.value = None
        self.modalPane.destroy()

    def _cancel(self, event=None):
        self.modalPane.destroy()
        
    def returnValue(self):
        self.master.wait_window(self.modalPane)
        return self.value

class GenBankFileDownload(Frame):
  
    def __init__(self, master, directory):
        self.master = master
        self.directory = directory
        self.list = []
        self.keyword = StringVar()
        self.keyword.set("")
        self.organism = StringVar()
        self.organism.set("")
        self.accession = StringVar()
        self.accession.set("")
        self.initUI()

    def initUI(self):
        self.modalPane = Toplevel(self.master, height=200, width=400)

        self.modalPane.transient(self.master)
        self.modalPane.grab_set()

        self.modalPane.bind("<Return>", self._download)
        self.modalPane.bind("<Escape>", self._cancel)

        self.modalPane.title("Search GenBank")

        Label(self.modalPane, text="Find GenBank entries on NCBI server").pack(padx=5, pady=5)
        
        searchFrame = Frame(self.modalPane)
        searchFrame.pack(side=TOP)
        searchFrame.grid_columnconfigure(0, minsize=100)
        Label(searchFrame, text="Keyword: ").grid(row=0,column=0, pady=3)
        keywordentry = Entry(searchFrame, textvariable=self.keyword, width=25)
        keywordentry.grid(row=0,column=1, pady=3)
        Label(searchFrame, text="Organism: ").grid(row=1,column=0, pady=3)
        organismentry = Entry(searchFrame, textvariable=self.organism, width=25)
        organismentry.grid(row=1,column=1, pady=3)
        Label(searchFrame, text="Accession: ").grid(row=2,column=0, pady=3)
        accessionentry = Entry(searchFrame, textvariable=self.accession, width=25)
        accessionentry.grid(row=2,column=1, pady=3)
        searchButton = Button(searchFrame, text="Search", command=self._search)
        searchButton.grid(row=3,column=2, pady=3)

        listFrame = Frame(self.modalPane)
        listFrame.pack(side=TOP, padx=20, pady=20)
        
        scrollBar = Scrollbar(listFrame)
        scrollBar.pack(side=RIGHT, fill=Y)
        self.listBox = Listbox(listFrame, selectmode=SINGLE, width=100, height=20)
        self.listBox.pack(side=LEFT, fill=Y)
        scrollBar.config(command=self.listBox.yview)
        self.listBox.config(yscrollcommand=scrollBar.set)
        self.list.sort()
        for item in self.list:
            self.listBox.insert(END, item)

        buttonFrame = Frame(self.modalPane)
        buttonFrame.pack(side=BOTTOM)

        chooseButton = Button(buttonFrame, text="Download", command=self._download)
        chooseButton.pack(side=LEFT)

        cancelButton = Button(buttonFrame, text="Cancel", command=self._cancel)
        cancelButton.pack(side=RIGHT)

    def _search(self):
        #Redirect if no input given
        if self.organism.get() == "" and self.keyword.get() == "" and self.accession.get() == "":
          showerror("Input error", "Please specify a search term first.")
          return
        #Build URL
        url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term="
        if self.organism.get() != "":
          url += "%22" + self.organism.get().replace(" ","%20") + "%22[porgn]"
        if self.keyword.get() != "":
          if self.organism.get() != "":
            url += '%20AND%20' + "%22" + self.keyword.get().replace(" ","%20") + "%22"
          else:
            url += "%22" + self.keyword.get().replace(" ","%20") + "%22"
        if self.accession.get() != "":
          url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term=" + self.accession.get() + "[accession]"
        url += "&retmode=text&retmax=251&tool=multigeneblast&email=m.h.medema@rug.nl"
        #Search all IDs matching the query
        try:
          req = urllib.request.urlopen(url)
          xmltext = req.read()
          xmltext = xmltext.partition("<IdList>")[2].partition("</IdList>")[0]
          ids = [ID.partition("</Id>")[0] for ID in xmltext.split("<Id>")[1:]]
          if len(ids) == 0:
            showerror('Error', 'No matches found to your search query.')
            return
        except:
          showerror('Error', 'Could not connect to NCBI server.\nPlease check your internet connection.')
          return
        #Show message if more than 250 found, asking user to further specify his search
        if len(ids) > 250:
          showwarning("Warning", "More than 250 matches found for your search query.\nOnly top 250 matches will be shown.")
          ids = ids[:-1]
        self.ids = ids
        #Fetch titles for all IDs
        url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&id=" + ",".join(ids) + "&retmode=xml&retmax=251&tool=multigeneblast&email=m.h.medema@rug.nl"
        #Open 'loading...' message
        loading = Toplevel(self.modalPane, height=200, width=400)
        loading.title("Loading")
        message = "Loading. Please wait..."
        a = Label(loading, text=message)
        a.grid(row=1,column=1, padx=50, pady=50)
        loading.bind("<Escape>", lambda e: "return")
        self.modalPane.update()
        try:
          req = urllib.request.urlopen(url)
          xmltext = req.read()
          if "<DocSum>" not in xmltext or '<Item Name="Title" Type="String">' not in xmltext:
            showerror('Error', 'Error in fetching match descriptions from NCBI server.')
            return
          entries = xmltext.split("<DocSum>")[1:]
          descriptions = [entry.partition('<Item Name="Title" Type="String">')[2].partition('</Item>')[0] for entry in entries]
        except:
          showerror('Error', 'Could not connect to NCBI server.\nPlease check your internet connection.')
          return
        self._clear()
        loading.destroy()
        for description in descriptions:
          self._insert(description)

    def _clear(self):
        self.listBox.delete(0, END)

    def _insert(self, string):
        self.listBox.insert(END, string)

    def _download(self, event=None):
        try:
          Selected = self.listBox.curselection()[0]
        except IndexError:
          self.value = None
        currentdir = os.getcwd()
        os.chdir(self.directory)
        ID = self.ids[int(Selected)]
        url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=" + ID + "&rettype=gbwithparts&retmode=text&tool=multigeneblast&email=m.h.medema@rug.nl"
        filename = tkinter.filedialog.asksaveasfilename(defaultextension="gbk", filetypes = [('GenBank files', '.gbk'), ('all files', '.*')])
        #Open 'loading...' message
        loading = Toplevel(self.modalPane, height=200, width=400)
        loading.title("Download progress")
        message = "Downloading.\nPlease wait..."
        a = Label(loading, text=message)
        a.grid(row=1,column=1, padx=50, pady=50)
        loading.bind("<Escape>", lambda e: "return")
        self.modalPane.update()
        try:
          req = urllib.request.urlopen(url)
        except:
          showerror('Error', 'File not found on server.\nPlease check your internet connection.')
          loading.destroy()
          return
        loading.protocol('WM_DELETE_WINDOW', lambda: cancel_download(req, loading))
        CHUNK = 128 * 1024
        with open(filename, 'wb') as fp:
          while True:
            message = "Downloading.\nPlease wait..."
            try:
              a = Label(loading, text=message)
              a.grid_remove()
              a.grid(row=1,column=1, padx=20, pady=20)
              self.modalPane.update()
              chunk = req.read(CHUNK)
              if not chunk:
                break
              fp.write(chunk)
            except:
              showerror('Error', 'Download interrupted.')
              return
        #Report download success
        self.modalPane.destroy()
        os.chdir(currentdir)
        showinfo("Download finished", "Download completed successfully.")

    def _cancel(self, event=None):
        self.modalPane.destroy()


class GeneSelectionFrame(Frame):
  
    def __init__(self, master, SpObj, positions, entrieslist):
        Frame.__init__(self, master)   
         
        self.master = master
        self.SpObj = SpObj
        self.positions = positions
        self.entrieslist = entrieslist
        self.initUI()
      
    def select(self):
        selectedgenes = ListBoxChoice(self, "Select genes", "Select query genes:", self.entrieslist).returnValue()
        if selectedgenes != None:
          self.SpObj.cstart.setval(0)
          self.SpObj.cend.setval(0)
          self.selected.set(selectedgenes)

    def initUI(self):
        self.grid(row=self.positions[0],column=self.positions[1], sticky=S)
        self.selected = StringVar()
        self.selected.set("<Select genes>")
        self.genes_entry = Entry(self, textvariable=self.selected, width=25)
        self.genes_entry.grid(row=0,column=0)
        self.genes_entry.bind("<FocusOut>", self.OnValidate)
        self.genes_entry.bind("<Return>", self.OnValidate)
        selectionBox = Button(self, text="Select genes", command=self.select)
        selectionBox.grid(row=0,column=1)
        
    def OnValidate(self, val):
        if self.selected.get() != None and self.selected.get() != "<Select genes>" and len(self.selected.get().split(";")) > 1:
            selectedgenes = self.selected.get()
            self.SpObj.cstart.setval(0)
            self.SpObj.cend.setval(0)
            self.selected.set(selectedgenes)
        else:
            self.selected.set("<Select genes>")

    def clear_selection(self):
        self.selected.set("<Select genes>")

    def getval(self):
        return self.selected.get()

class MessageBox(object):
    def __init__(self, frame=None, master=None, title=None, message=None, list=[]):
        self.frame = frame
        self.master = master
        self.value = None
        self.errormessage = ""
        self.list = list[:]
        
        self.modalPane = Toplevel(self.frame, height=800, width=300)

        self.modalPane.transient(self.frame)
        self.modalPane.grab_set()

        self.modalPane.bind("<Escape>", self._ok)
        self.modalPane.protocol('WM_DELETE_WINDOW', self._ok)

        if title:
            self.modalPane.title(title)

        #define a new frame and put a text area in it
        self.messageFrame = Frame(self.modalPane)
        self.messageFrame.pack(side=TOP)
        self.messageFrame.text = Text(self.messageFrame,height=40,width=100,background='white',state=DISABLED)

        # put a scroll bar in the frame
        self.messageFrame.scroll = Scrollbar(self.messageFrame)
        self.messageFrame.text.configure(yscrollcommand=self.messageFrame.scroll.set)

        #pack everything
        self.messageFrame.text.pack(side=LEFT)
        self.messageFrame.scroll.pack(side=RIGHT,fill=Y)

        #listFrame = Frame(self.modalPane)
        #listFrame.pack(side=TOP, padx=5, pady=5)


    def text_insert(self, text):
        self.messageFrame.text.see(END)
        self.messageFrame.text.config(state=NORMAL)
        self.messageFrame.text.insert(INSERT, text)
        self.messageFrame.text.config(state=DISABLED)
        
    def add_ok_button(self):
        self.buttonFrame = Frame(self.modalPane)
        self.buttonFrame.pack(side=BOTTOM)

        self.okButton = Button(self.buttonFrame, text="OK", command=self._ok, width=10)
        self.okButton.pack()

    def add_error_button(self):
        self.buttonFrame = Frame(self.modalPane)
        self.buttonFrame.pack(side=BOTTOM)

        self.errorButton = Button(self.buttonFrame, text="Send Error Report", command=self.sendreport, width=20)
        self.errorButton.pack()

    def change_errormessage(self, errormessage):
      self.errormessage = errormessage

    def sendreport(self, event=None):
      try:
        webbrowser.open("mailto:multigeneblast@gmail.com?SUBJECT=ErrorReport&BODY=" + urllib.parse.quote(self.errormessage.encode("utf8")))
      except:
        webbrowser.open("sourceforge.net/tracker/?func=add&group_id=565495&atid=2293721")
      else:
        pass

    def _ok(self, event=None):
        #self.master._cancel() #Would destroy parent window as well
        self.modalPane.destroy()

class MakeDatabase(Frame):
  
    def __init__(self, master, directory):
        self.master = master
        self.files = []
        self.dbname = StringVar()
        self.directory = directory
        self.dbname.set("<Enter a name for your __database_file_label>")
        self.initUI()

    def initUI(self):
        self.modalPane = Toplevel(self.master, height=200, width=400)

        self.modalPane.transient(self.master)
        self.modalPane.grab_set()

        self.modalPane.bind("<Return>", self._makedb)
        self.modalPane.bind("<Escape>", self._cancel)

        self.modalPane.title("Make MGB __database_file_label from file")

        searchFrame = Frame(self.modalPane)
        searchFrame.pack(side=TOP)
        searchFrame.grid_columnconfigure(0, minsize=100)
        Label(searchFrame, text="Database name: ").grid(row=0,column=0, pady=3)
        self.dbname_entry = Entry(searchFrame, textvariable=self.dbname, width=50)
        self.dbname_entry.grid(row=0,column=1, pady=3)
        self.dbname_entry.bind("<FocusOut>", self.OnValidate)
        self.dbname_entry.bind("<Return>", self.OnValidate)

        text = Label(searchFrame, text="\nMake raw nucleotide __database_file_label for tblastn-searching:")
        text.grid(row=1,column=0, pady=5)
        self.dbtype = CheckBox(searchFrame, [1,1], "")

        searchButton = Button(searchFrame, text="Add files", command=self._addfile)
        searchButton.grid(row=2,column=0, pady=3)

        listFrame = Frame(self.modalPane)
        listFrame.pack(side=TOP, padx=20, pady=5)
        scrollBar = Scrollbar(listFrame)
        scrollBar.pack(side=RIGHT, fill=Y)
        self.listBox = Listbox(listFrame, selectmode=SINGLE, width=100, height=10)
        self.listBox.pack(side=LEFT, fill=Y)
        scrollBar.config(command=self.listBox.yview)
        self.listBox.config(yscrollcommand=scrollBar.set)
        self.listBox.bind("<Delete>", self._remove)
        self.files.sort()
        for item in self.files:
            self.listBox.insert(END, item)

        clearButton = Button(self.modalPane, text="Clear list", command=self._clear)
        clearButton.pack(side=LEFT)

        buttonFrame = Frame(self.modalPane)
        buttonFrame.pack(side=BOTTOM)
        chooseButton = Button(buttonFrame, text="Make __database_file_label", command=self._makedb)
        chooseButton.pack(side=LEFT)
        cancelButton = Button(buttonFrame, text="Cancel", command=self._cancel)
        cancelButton.pack(side=RIGHT)

    def OnValidate(self, val):
        forbiddencharacters = ["'",'"','=',':',';','-','[',']','>','<','|','\\',"/",'*',',','?',')','(','^','#','!','`','~','+','{','}','@','$','%','&']
        for char in forbiddencharacters:
          if self.dbname.get() != "<Enter a name for your __database_file_label>":
            if char in self.dbname.get():
              showerror("Error", "Forbidden character found in specified __database_file_label name: " + char)
              self.dbname.set("<Enter a name for your __database_file_label>")
              return
        if self.dbname.get() != None:
          if self.dbname.get() != "<Enter a name for your __database_file_label>":
            dbname = self.dbname.get().replace(" ","_")
            self.dbname.set(dbname)
        else:
          self.dbname.set("<Enter a name for your __database_file_label>")

    def _addfile(self):
        if self.dbtype.getval() == 0:
          filenames = tkinter.filedialog.askopenfilename(multiple=True, filetypes=(("GenBank files",('*.gbk','*.gb','*.genbank')),("EMBL files",('*.embl','*.emb'))), title="Select files", initialdir=self.directory)
        else:
          filenames = tkinter.filedialog.askopenfilename(multiple=True, filetypes=(("FASTA files",('*.fasta','*.fa','*.fas')), ("GenBank files",('*.gbk','*.gb','*.genbank')),("EMBL files",('*.embl','*.emb'))), title="Select files", initialdir=self.directory)
        if sys.platform == ('win32'):
          if "} {" in filenames:
            filenames = [str(filename) for filename in filenames[1:-1].split("} {")]
          else:
              if len(filenames) > 2 and filenames[2] == ":":
                filenames = [str(filenames[1:-1])]
              else:
                filenames = [str(filenames)]
        if sys.platform == ('darwin') or sys.platform == ('linux2'):
            filenames = list(filenames)
        for filename in filenames:
          if filename == "":
            return
          root, fext = os.path.splitext(filename)
          if self.dbtype.getval() == 0:
            extensions = [".gbk",".gb",".genbank",".embl",".emb"]
          else:
            extensions = [".gbk",".gb",".genbank",".embl",".emb",".fasta",".fa",".fas"]
          if fext.lower() not in extensions:
            if self.dbtype.getval() == 0:
              showerror("Error", filename + " not added, as it does not have a GenBank or EMBL file extension.")
            else:
              showerror("Error", filename + " not added, as it does not have a FASTA, GenBank or EMBL file extension.")
            continue
          if filename not in self.files:
            self.files.append(filename)
            self._insert(filename)

    def _clear(self):
        self.listBox.delete(0, END)
        self.files = []

    def _insert(self, string):
        self.listBox.insert(END, string)

    def _remove(self, idx):
        pos = self.listBox.curselection()
        del self.files[int(pos[0])]
        self.listBox.delete(pos, pos)

    def _makedb(self, event=None):
        global APPDATA
        #Final check for correct __database_file_label name
        if not validname(self.dbname.get()):
          showerror("Error", "Please select a valid __database_file_label name.")
          return
        if len(self.files) == 0:
          showerror("Error", "Please select input files for your __database_file_label first.")
          return
        #Remove spaces from dbname if this was not done automatically already; check for overwriting
        self.dbname.set(self.dbname.get().replace(" ","_"))
        if self.dbname.get() + ".pal" in os.listdir(APPDATA) or self.dbname.get() + ".nal" in os.listdir(APPDATA):
          answer = askyesno('Database name exists', 'A __database_file_label with this name already exists. Overwrite?')
          if not answer:
            return
        #Open up toplevel that reports progress
        try:
          self.makedb_outbox = MessageBox(self.master, self, "Creating Database...")
        except(TclError):
          pass
        
        #Create __database_file_label
        currentdir = os.getcwd()
        os.chdir(self.directory)
        if self.dbtype.getval() == 0:
          try:
            runmakedb(self.master, self.makedb_outbox, self.files, self.dbname.get())
          except(TclError):
            pass
          except:
            showerror("Failure", "Database creation failed.")
            os.chdir(currentdir)
            self.makedb_outbox.add_ok_button()
            self.modalPane.destroy()
            return
        else:
          try:
            runmakendb(self.master, self.makedb_outbox, self.files, self.dbname.get())
          except(TclError):
            pass
          except:
            showerror("Failure", "Database creation failed.")
            os.chdir(currentdir)
            self.makedb_outbox.add_ok_button()
            self.modalPane.destroy()
            return
        #Report success
        try:
          if self.dbtype.getval() == 0:
            self.makedb_outbox.text_insert("Database created.\nYou can now use this __database_file_label by selecting '" + self.dbname.get() + ".pal' under 'Select __database_file_label' in the 'File' menu.\n")
          else:
            self.makedb_outbox.text_insert("Database created.\nYou can now use this __database_file_label by selecting '" + self.dbname.get() + ".nal' under 'Select __database_file_label' in the 'File' menu.\n")
          self.makedb_outbox.add_ok_button()
          self.modalPane.destroy()
        except(TclError):
          pass

    def _cancel(self, event=None):
        self.modalPane.destroy()

def log(message, frame, makedb_outbox, exit=False, retcode=1):
    try:
      makedb_outbox.text_insert(message + "\n")
    except(TclError):
      pass
    frame.update()
    if exit:
      sys.exit(retcode)

def runmakedb(frame, makedb_outbox, inputfiles, dbname):
    if len(inputfiles) < 10:
      log("Creating __database_file_label. Please be patient, this can take a number of minutes...\n", frame, makedb_outbox)
    else:
      log("Creating __database_file_label. Please be patient, this can take a number of hours...\n", frame, makedb_outbox)
    
    #Create FASTA __database_file_label
    log("Creating FASTA file...", frame, makedb_outbox)
    descriptions = parse_gbk_embl(inputfiles, dbname, frame, makedb_outbox, GUI="y")

    #Create genbank_mf_all.txt file
    log("Generating " + dbname + "_all.txt file", frame, makedb_outbox)
    fasta_names(dbname + "_dbbuild.fasta", dbname + "_all.txt")

    #Create genbank_mf_all_descrs.txt file
    log("Generating " + dbname + "_all_descr.txt file", frame, makedb_outbox)
    descrfile = open(dbname + "_all_descrs.txt","w")
    descrkeys = list(descriptions.keys())
    descrkeys.sort()
    for key in descrkeys:
      descrfile.write(key + "\t" + descriptions[key] + "\n")
    descrfile.close()

    #Create Blast __database_file_label and TAR archives
    log("Creating BLAST __database_file_label...", frame, makedb_outbox)
    if "\n" not in open(dbname + "_dbbuild.fasta","r").read():
      log("Error making BLAST __database_file_label; no suitable sequences found in input.", frame, makedb_outbox, exit=True)
    make_blast_db(dbname, dbname + "_dbbuild.fasta", frame, GUI="n")
    log("Creating gene coordinates file...", frame, makedb_outbox)
    generate_genecords_tar(dbname, frame, makedb_outbox, GUI="y")
    log("Creating gene sources file...", frame, makedb_outbox)
    generate_proteininfo_tar(dbname, frame, makedb_outbox, GUI="y")

    #Clean up
    log("Cleaning up...\n", frame, makedb_outbox)
    clean_up(dbname, frame, makedb_outbox, GUI="y")

def runmakendb(frame, makedb_outbox, inputfiles, dbname):
    if len(inputfiles) < 10:
      log("Creating __database_file_label. Please be patient, this can take a number of minutes...\n", frame, makedb_outbox)
    else:
      log("Creating __database_file_label. Please be patient, this can take a number of hours...\n", frame, makedb_outbox)

    #Create FASTA __database_file_label
    log("Creating FASTA file...", frame, makedb_outbox)
    descriptions = nucparse_gbk_embl_fasta(inputfiles, dbname, frame, makedb_outbox, GUI="y")

    #Create genbank_mf_all_descr.txt file
    log("Generating _all.txt file", frame, makedb_outbox)
    outfile = open(dbname + "_all_descrs.txt","w")
    for key in list(descriptions.keys()):
      outfile.write(key + "\t" + descriptions[key] + "\n")
    outfile.close()

    #Create Blast __database_file_label
    if "\n" not in open(dbname + "_dbbuild.fasta","r").read():
      clean_up(dbname, dbtype="nucl")
      log("Error making BLAST __database_file_label; no suitable sequences found in input.", frame, makedb_outbox, exit=True)
    make_blast_db(dbname, dbname + "_dbbuild.fasta", frame, GUI="y", dbtype="nucl")

    #Clean up
    log("Cleaning up...\n", frame, makedb_outbox)
    clean_up(dbname, frame, makedb_outbox, GUI="y", dbtype="nucl")

class MakeOnlineDatabase(Frame):
  
    def __init__(self, master, directory):
        self.master = master
        self.directory = directory
        self.list = []
        self.selected = []
        self.keyword = StringVar()
        self.keyword.set("")
        self.organism = StringVar()
        self.organism.set("")
        self.accession = StringVar()
        self.accession.set("")
        self.dbname = StringVar()
        self.dbname.set("<Enter a name for your __database_file_label>")
        self.descriptions = []
        self.initUI()

    def initUI(self):
        self.modalPane = Toplevel(self.master, height=200, width=400)

        self.modalPane.transient(self.master)
        self.modalPane.grab_set()

        self.modalPane.bind("<Return>", self._download_makedb)
        self.modalPane.bind("<Escape>", self._cancel)

        self.modalPane.title("Search GenBank")

        Label(self.modalPane, text="Find GenBank entries on NCBI server").grid(row=0, column=0)
        
        searchFrame = Frame(self.modalPane)
        searchFrame.grid(row=1, column=0)
        searchFrame.grid_columnconfigure(0, minsize=100)
        Label(searchFrame, text="Keyword: ").grid(row=0,column=0, pady=3)
        keywordentry = Entry(searchFrame, textvariable=self.keyword, width=25)
        keywordentry.grid(row=0,column=1, pady=3)
        Label(searchFrame, text="Organism: ").grid(row=1,column=0, pady=3)
        organismentry = Entry(searchFrame, textvariable=self.organism, width=25)
        organismentry.grid(row=1,column=1, pady=3)
        Label(searchFrame, text="Accession: ").grid(row=2,column=0, pady=3)
        accessionentry = Entry(searchFrame, textvariable=self.accession, width=25)
        accessionentry.grid(row=2,column=1, pady=3)
        searchButton = Button(searchFrame, text="Search", command=self._search)
        searchButton.grid(row=3,column=2, pady=3)

        listFrame = Frame(self.modalPane)
        listFrame.grid(row=2, column=0, padx=20, pady=5)
        Label(listFrame, text="Search results:").pack(side=TOP)
        listInnerFrame = Frame(listFrame)
        listInnerFrame.pack(side=BOTTOM)
        scrollBar = Scrollbar(listInnerFrame)
        scrollBar.pack(side=RIGHT, fill=Y)
        self.listBox = Listbox(listInnerFrame, selectmode=EXTENDED, width=100, height=8)
        self.listBox.pack(side=LEFT, fill=Y)
        scrollBar.config(command=self.listBox.yview)
        self.listBox.config(yscrollcommand=scrollBar.set)
        self.list.sort()
        for item in self.list:
            self.listBox.insert(END, item)
        
        buttonFrame1 = Frame(self.modalPane, width=100)
        buttonFrame1.grid(row=3, column=0)
        selectButton = Button(buttonFrame1, text="Select", command=self._insert_selected)
        selectButton.pack(side=RIGHT)

        Label(self.modalPane, text=" ").grid(row=4, column=0, sticky=S)
        Label(self.modalPane, text="Selected entries:").grid(row=5, column=0, sticky=S)
        
        listFrame2 = Frame(self.modalPane)
        listFrame2.grid(row=6, column=0, padx=20, pady=5)
        scrollBar2 = Scrollbar(listFrame2)
        scrollBar2.pack(side=RIGHT, fill=Y)
        self.listBox2 = Listbox(listFrame2, selectmode=SINGLE, width=100, height=8)
        self.listBox2.pack(side=LEFT, fill=Y)
        scrollBar2.config(command=self.listBox2.yview)
        self.listBox2.config(yscrollcommand=scrollBar2.set)
        self.listBox2.bind("<Delete>", self._remove_selected)
        self.selected.sort()
        for item in self.selected:
            self.listBox2.insert(END, item)

        dbnameFrame = Frame(self.modalPane)
        dbnameFrame.grid(row=7,column=0, pady=5)
        Label(dbnameFrame, text="Database name: ").pack(side=LEFT)
        self.dbname_entry = Entry(dbnameFrame, textvariable=self.dbname, width=50)
        self.dbname_entry.pack(side=RIGHT)
        self.dbname_entry.bind("<FocusOut>", self.OnValidate)
        self.dbname_entry.bind("<Return>", self.OnValidate)

        checkboxFrame = Frame(self.modalPane)
        checkboxFrame.grid(row=8,column=0)
        text = Label(checkboxFrame, text="Make raw nucleotide __database_file_label for tblastn-searching:")
        text.pack(side=LEFT, pady=10)
        self.dbtype = CheckBox(checkboxFrame, [0,0], "")
        self.dbtype.pack(side=RIGHT)

        buttonFrame2 = Frame(self.modalPane)
        buttonFrame2.grid(row=9, column=0)

        clearButton = Button(buttonFrame2, text="Clear selection", command=self._clear_selected)
        clearButton.pack(side=LEFT)

        chooseButton = Button(buttonFrame2, text="Download and create __database_file_label", command=self._download_makedb)
        chooseButton.pack(side=LEFT)

        cancelButton = Button(buttonFrame2, text="Cancel", command=self._cancel)
        cancelButton.pack(side=RIGHT)

    def OnValidate(self, val):
        forbiddencharacters = ["'",'"','=',':',';','-','[',']','>','<','|','\\',"/",'*',',','?',')','(','^','#','!','`','~','+','{','}','@','$','%','&']
        for char in forbiddencharacters:
          if self.dbname.get() != "<Enter a name for your __database_file_label>":
            if char in self.dbname.get():
              showerror("Error", "Forbidden character found in specified __database_file_label name: " + char)
              self.dbname.set("<Enter a name for your __database_file_label>")
              return
        if self.dbname.get() != None:
          if self.dbname.get() != "<Enter a name for your __database_file_label>":
            dbname = self.dbname.get().replace(" ","_")
            self.dbname.set(dbname)
        else:
          self.dbname.set("<Enter a name for your __database_file_label>")

    def _search(self):
        #Redirect if no input given
        if self.organism.get() == "" and self.keyword.get() == "" and self.accession.get() == "":
          showerror("Input error", "Please specify a search term first.")
          return
        #Build URL
        url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term="
        if self.organism.get() != "":
          url += "%22" + self.organism.get().replace(" ","%20") + "%22[porgn]"
        if self.keyword.get() != "":
          if self.organism.get() != "":
            url += '%20AND%20' + "%22" + self.keyword.get().replace(" ","%20") + "%22"
          else:
            url += "%22" + self.keyword.get().replace(" ","%20") + "%22"
        if self.accession.get() != "":
          url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term=" + self.accession.get() + "[accession]"
        url += "&retmode=text&retmax=251&tool=multigeneblast&email=m.h.medema@rug.nl"
        #Search all IDs matching the query
        try:
          req = urllib.request.urlopen(url)
          xmltext = req.read()
          xmltext = xmltext.partition("<IdList>")[2].partition("</IdList>")[0]
          ids = [ID.partition("</Id>")[0] for ID in xmltext.split("<Id>")[1:]]
          if len(ids) == 0:
            showerror('Error', 'No matches found to your search query.')
            return
        except:
          showerror('Error', 'Could not connect to NCBI server.\nPlease check your internet connection.')
          return
        #Show message if more than 250 found, asking user to further specify his search
        warning250 = False
        if len(ids) > 250:
          warning250 = True
          ids = ids[:-1]
        self.ids = ids
        #Fetch titles for all IDs
        url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&id=" + ",".join(ids) + "&retmode=xml&retmax=251&tool=multigeneblast&email=m.h.medema@rug.nl"
        #Open 'loading...' message
        loading = Toplevel(self.modalPane, height=200, width=400)
        loading.title("Loading")
        message = "Loading. Please wait..."
        a = Label(loading, text=message)
        a.grid(row=1,column=1, padx=50, pady=50)
        loading.bind("<Escape>", lambda e: "return")
        self.modalPane.update()
        try:
          req = urllib.request.urlopen(url)
          xmltext = req.read()
          if "<DocSum>" not in xmltext or '<Item Name="Title" Type="String">' not in xmltext:
            showerror('Error', 'Error in fetching match descriptions from NCBI server.')
            return
          entries = xmltext.split("<DocSum>")[1:]
          self.descriptions = [entry.partition('<Item Name="Title" Type="String">')[2].partition('</Item>')[0] for entry in entries]
        except:
          showerror('Error', 'Could not connect to NCBI server.\nPlease check your internet connection.')
          return
        self._clear()
        loading.destroy()
        if warning250:
          showwarning("Warning", "More than 250 matches found for your search query.\nOnly top 250 matches will be shown.")
        for description in self.descriptions:
          self._insert(description)

    def _clear(self):
        self.listBox.delete(0, END)

    def _insert(self, string):
        self.listBox.insert(END, string)

    def _clear_selected(self):
        self.selected = []
        self.listBox2.delete(0, END)

    def _remove_selected(self, idx):
        pos = self.listBox2.curselection()
        self.listBox2.delete(pos, pos)
        ID = self.ids[int(pos)]
        del self.selected[self.selected.index(ID)]

    def _insert_selected(self):
        try:
          Selected = self.listBox.curselection()
        except IndexError:
          return
        for Sel in Selected:
          ID = self.ids[int(Sel)]
          if ID not in self.selected:
            self.listBox2.insert(END, self.descriptions[int(Sel)])
            self.selected.append(ID)     

    def _download_makedb(self, event=None):
        global APPDATA
        currentdir = os.getcwd()
        os.chdir(self.directory)
        if len(self.selected) == 0:
          showerror("Error", "No Entries selected.")
          return
        #Final check for correct __database_file_label name
        if not validname(self.dbname.get()):
          showerror("Error", "Please select a valid __database_file_label name.")
          return
        #Remove spaces from dbname if this was not done automatically already
        self.dbname.set(self.dbname.get().replace(" ","_"))
        #Check for overwriting
        if self.dbname.get() + ".pal" in os.listdir(APPDATA) or self.dbname.get() + ".nal" in os.listdir(APPDATA):
          answer = askyesno('Database name exists', 'A __database_file_label with this name already exists. Overwrite?')
          if not answer:
            return
        IDs = ",".join(self.selected)
        url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=" + IDs + "&rettype=gbwithparts&retmode=text&tool=multigeneblast&email=m.h.medema@rug.nl"
        filename = self.dbname.get() + ".gbk"
        #Open 'loading...' message
        loading = Toplevel(self.modalPane, height=200, width=400)
        loading.title("Download progress")
        message = "Downloading.\nPlease wait..."
        a = Label(loading, text=message)
        a.grid(row=1,column=1, padx=50, pady=50)
        loading.bind("<Escape>", lambda e: "return")
        self.modalPane.update()
        try:
          req = urllib.request.urlopen(url)
        except:
          showerror('Error', 'File not found on server.\nPlease check your internet connection.')
          loading.destroy()
          return
        loading.protocol('WM_DELETE_WINDOW', lambda: cancel_download(req, loading))
        CHUNK = 128 * 1024
        with open(filename, 'wb') as fp:
          while True:
            message = "Downloading.\nPlease wait..."
            try:
              a = Label(loading, text=message)
              a.grid_remove()
              a.grid(row=1,column=1, padx=20, pady=20)
              self.modalPane.update()
              chunk = req.read(CHUNK)
              if not chunk:
                break
              fp.write(chunk)
            except:
              showerror('Error', 'Download interrupted.')
              return
        #Report download success
        self.modalPane.destroy()
        showinfo("Download finished", "Download completed successfully.\nProceeding with creating the __database_file_label...\n")

        #Open up toplevel that reports progress
        try:
          self.makedb2_outbox = MessageBox(self.master, self, "Creating Database...")
        except(TclError):
          pass
        
        #Create __database_file_label
        if self.dbtype.getval() == 0:
          try:
            runmakedb(self.master, self.makedb2_outbox, [self.dbname.get() + ".gbk"], self.dbname.get())
          except(TclError):
            pass
          except:
            showerror("Failure", "Database creation failed.")
            os.chdir(currentdir)
            self.makedb2_outbox.add_ok_button()
            self.modalPane.destroy()
            return
        else:
          try:
            runmakendb(self.master, self.makedb2_outbox, [self.dbname.get() + ".gbk"], self.dbname.get())
          except(TclError):
            pass
          except:
            showerror("Failure", "Database creation failed.")
            os.chdir(currentdir)
            self.makedb2_outbox.add_ok_button()
            self.modalPane.destroy()
            return
        #Report success
        try:
          if self.dbtype.getval() == 0:
            self.makedb2_outbox.text_insert("Database created.\nYou can now use this __database_file_label by selecting '" + self.dbname.get() + ".pal' under 'Select __database_file_label' in the 'File' menu.\n")
          else:
            self.makedb2_outbox.text_insert("Database created.\nYou can now use this __database_file_label by selecting '" + self.dbname.get() + ".nal' under 'Select __database_file_label' in the 'File' menu.\n")
          self.makedb2_outbox.add_ok_button()
          self.modalPane.destroy()
        except(TclError):
          pass
        os.chdir(currentdir)
        
    def _cancel(self, event=None):
        self.modalPane.destroy()


class MakeGenBankDatabase(Frame):
  
    def __init__(self, master, directory):
        self.master = master
        self.directory = directory
        self.gbdivisions = ["BCT: Bacteria", "PLN: Plants, Fungi and Algae", "PHG: Bacteriophages", "SYN: Synthetic sequences", "PAT: Patent sequences", "ENV: Environmental sequences", "WGS: whole genome shotgun sequencing data", "CON: Unfinished genomes"]
        self.selected = []
        self.dbname = StringVar()
        self.dbname.set("<Enter a name for your __database_file_label>")
        self.initUI()

    def initUI(self):
        self.modalPane = Toplevel(self.master, height=200, width=400)

        self.modalPane.transient(self.master)
        self.modalPane.grab_set()

        self.modalPane.bind("<Return>", self._makedb)
        self.modalPane.bind("<Escape>", self._cancel)

        self.modalPane.title("Make MGB __database_file_label from GenBank divisions")

        searchFrame = Frame(self.modalPane)
        searchFrame.pack(side=TOP)
        searchFrame.grid_columnconfigure(0, minsize=100)
        Label(searchFrame, text="Database name: ").grid(row=0,column=0, pady=3)
        self.dbname_entry = Entry(searchFrame, textvariable=self.dbname, width=30)
        self.dbname_entry.grid(row=0,column=1, pady=3)
        self.dbname_entry.bind("<FocusOut>", self.OnValidate)
        self.dbname_entry.bind("<Return>", self.OnValidate)

        listFrame = Frame(self.modalPane)
        listFrame.pack(side=TOP, padx=20, pady=20)
        scrollBar = Scrollbar(listFrame)
        scrollBar.pack(side=RIGHT, fill=Y)
        self.listBox = Listbox(listFrame, selectmode=EXTENDED, width=45, height=8)
        self.listBox.pack(side=LEFT, fill=Y)
        scrollBar.config(command=self.listBox.yview)
        self.listBox.config(yscrollcommand=scrollBar.set)
        for item in self.gbdivisions:
            self.listBox.insert(END, item)

        buttonFrame = Frame(self.modalPane)
        buttonFrame.pack(side=BOTTOM)
        chooseButton = Button(buttonFrame, text="Make __database_file_label", command=self._makedb)
        chooseButton.pack(side=LEFT)
        cancelButton = Button(buttonFrame, text="Cancel", command=self._cancel)
        cancelButton.pack(side=RIGHT)

        checkboxFrame = Frame(self.modalPane)
        checkboxFrame.pack(side=BOTTOM)
        text = Label(checkboxFrame, text="Make raw nucleotide __database_file_label for tblastn-searching:")
        text.pack(side=LEFT, pady=5)
        self.dbtype = CheckBox(checkboxFrame, [0,0], "")
        self.dbtype.pack(side=RIGHT)


    def OnValidate(self, val):
        forbiddencharacters = ["'",'"','=',':',';','-','[',']','>','<','|','\\',"/",'*',',','?',')','(','^','#','!','`','~','+','{','}','@','$','%','&']
        for char in forbiddencharacters:
          if self.dbname.get() != "<Enter a name for your __database_file_label>":
            if char in self.dbname.get():
              showerror("Error", "Forbidden character found in specified __database_file_label name: " + char)
              self.dbname.set("<Enter a name for your __database_file_label>")
              return
        if self.dbname.get() != None:
          if self.dbname.get() != "<Enter a name for your __database_file_label>":
            dbname = self.dbname.get().replace(" ","_")
            self.dbname.set(dbname)
        else:
          self.dbname.set("<Enter a name for your __database_file_label>")

    def _makedb(self, event=None):
        global APPDATA
        #Final check for correct __database_file_label name
        selected = [self.gbdivisions[int(idx)].partition(":")[0] for idx in self.listBox.curselection()]
        if not validname(self.dbname.get()):
          showerror("Error", "Please select a valid __database_file_label name.")
          return
        if len(selected) == 0:
          showerror("Error", "Please select input files for your __database_file_label first.")
          return
        #Remove spaces from dbname if this was not done automatically already; check for overwriting
        self.dbname.set(self.dbname.get().replace(" ","_"))
        if self.dbname.get() + ".pal" in os.listdir(APPDATA) or self.dbname.get() + ".nal" in os.listdir(APPDATA):
          answer = askyesno('Database name exists', 'A __database_file_label with this name already exists. Overwrite?')
          if not answer:
            return
        #Check that user is really ok with doing this
        answer = askyesno('Time', 'Database downloading and creation will take minutes to hours, depending on the divisions selected.\nAre you sure you want to continue?')
        if not answer:
          return
        #Check to prevent overwriting
        currentdir = os.getcwd()
        os.chdir(self.directory)
        if "genbank" in os.listdir(APPDATA):
          user_response = askyesno("Warning", "A folder named 'genbank' (which is used by MultiGeneBlast to store downloaded GenBank files) is already present,\n probably from an earlier search that was interrupted.\n Overwrite?")
          if not user_response:
            return
          else:
            shutil.rmtree("genbank")
        if "WGS" in selected and "wgs" in os.listdir(APPDATA):
          user_response = askyesno("Warning", "A folder named 'wgs' (which is used by MultiGeneBlast to store downloaded WGS GenBank files) is already present,\n probably from an earlier search that was interrupted.\n Overwrite?")
          if not user_response:
            return
          else:
            shutil.rmtree("wgs")
        #Open up toplevel that reports progress
        try:
          self.makegbdb_outbox = MessageBox(self.master, self, "Creating Database...")
          self.makegbdb_outbox.text_insert("GenBank divisions to be downloaded: " + ",".join(selected) + "\n")

          #Create __database_file_label
          if self.dbtype.getval() == 0:
            makegbdb.maingui(selected, self.dbname.get(), self.master, self.makegbdb_outbox)
          else:
            makegbndb.maingui(selected, self.dbname.get(), self.master, self.makegbdb_outbox)
        except(TclError):
          pass
        except:
          os.chdir(currentdir)
          self.makegbdb_outbox.add_ok_button()
          self.modalPane.destroy()
          return
        if not (self.dbname.get() + ".pal" in os.listdir(APPDATA) or self.dbname.get() + ".nal" in os.listdir(APPDATA)):
          self.makegbdb_outbox.text_insert("Error: __database_file_label creation failed due to unknown error.\n")
          os.chdir(currentdir)
          self.makegbdb_outbox.add_ok_button()
          self.modalPane.destroy()
          return
        #Report success
        if self.dbtype.getval() == 0:
          self.makegbdb_outbox.text_insert("Database created.\nYou can now use this __database_file_label by selecting '" + self.dbname.get() + ".pal' under 'Select __database_file_label' in the 'File' menu.\n")
        else:
          self.makegbdb_outbox.text_insert("Database created.\nYou can now use this __database_file_label by selecting '" + self.dbname.get() + ".nal' under 'Select __database_file_label' in the 'File' menu.\n")
        self.makegbdb_outbox.add_ok_button()
        self.modalPane.destroy()
        os.chdir(currentdir)

    def _cancel(self, event=None):
        self.modalPane.destroy()
