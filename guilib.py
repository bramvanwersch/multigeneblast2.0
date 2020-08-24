#!/usr/bin/env python
## Copyright (c) 2012 Marnix H. Medema
## Department of Microbial Physiology / Groningen Bioinformatics Centre
## University of Groningen
## License: GNU General Public License v3 or later
## A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

import subprocess
import webbrowser
import time
import tkinter.filedialog
from tkinter.messagebox import askyesno, showerror, showwarning, showinfo
from tkinter.ttk import Frame, Label, Scale, Style
from tkinter import Tk, BOTH, IntVar, Checkbutton, Spinbox, Listbox, StringVar, Entry, Button, Toplevel, TOP, RIGHT, LEFT, BOTTOM, Y, END, EXTENDED, Scrollbar, Text, NORMAL, INSERT, DISABLED, SINGLE, S, N, W, TclError
import os
import urllib.request, urllib.error, urllib.parse

from utilities import ILLEGAL_CHARACTERS
from constants import APPDATA, GENBANK_EXTENSIONS, EMBL_EXTENSIONS, get_mgb_path, CHUNK, OUT_FOLDER_NAME

MGBPATH = get_mgb_path()


def about():
    webbrowser.open('http://multigeneblast.sourceforge.net/')


class ScaleBar(Frame):
  
    def __init__(self, parent, SpObj, positions, minimum, maximum, default="0", scale_command=None, input_type="int"):
        Frame.__init__(self, parent)   
         
        self.parent = parent
        self.SpObj = SpObj
        self.positions = positions
        self.minimum = minimum
        self.maximum = maximum
        self.default = default
        self.command = scale_command
        self.input_type = input_type
        self.initUI()
        
    def initUI(self):
      
        #self.parent.title("Scale")
        #self.style = Style()
        #self.style.theme_use("default")
        
        #self.pack(fill=BOTH, expand=1)
        
        self.grid(row=self.positions[0],column=self.positions[1], sticky=W)
        
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

    def set_scale(self, start, end, value):
        self.scale.configure(from_=start)
        self.scale.configure(to=end)
        self.var.set(str(value))
        self.scale.set(value)
        
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
        if self.command != None:
            self.command()
        
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

        self.grid(row=self.positions[0],column=self.positions[1], sticky=W, padx=25)
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

        self.grid(row=self.positions[0],column=self.positions[1], sticky=W)
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

    def set_selectable_genes(self, genes):
        self.entrieslist = genes
        self.clear_selection()

    def initUI(self):
        self.grid(row=self.positions[0],column=self.positions[1], sticky=W)
        self.selected = StringVar()
        self.selected.set("<Select genes>")
        self.genes_entry = Entry(self, textvariable=self.selected, width=25)
        self.genes_entry.grid(row=0,column=0)
        self.genes_entry.bind("<FocusOut>", self.OnValidate)
        self.genes_entry.bind("<Return>", self.OnValidate)
        selectionBox = Button(self, text="Select genes", command=self.select)
        selectionBox.grid(row=0,column=1, padx=10)
        
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

class MessageBox:
    def __init__(self, frame=None, master=None, title=None, message=None, list=[], outdir=None):
        self.frame = frame
        self.master = master
        self.value = None
        self.errormessage = ""
        self.list = list[:]
        self._log_loc = outdir
        if outdir != None:
            self._log_loc = outdir + os.sep + OUT_FOLDER_NAME + os.sep + "run.log"
        
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

        #listFrame = Frame(self.modal_pane)
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

        self.errorButton = Button(self.buttonFrame, text="Send Error Report",
                                  command=self.sendreport, width=20)
        self.errorButton.pack()

    def change_errormessage(self, errormessage):
        self.errormessage = errormessage

    def sendreport(self, event=None):
        #TODO make this more proper by sending automatically
        try:

            webbrowser.open("mailto:multigeneblast@gmail.com?SUBJECT=ErrorReport&BODY=" + urllib.parse.quote(self.errormessage.encode("utf8")))
        except:
            webbrowser.open("sourceforge.net/tracker/?func=add&group_id=565495&atid=2293721")
        else:
            pass

    def _ok(self, event=None):
        self.modalPane.destroy()

class GenBankFileDownload(Frame):

    def __init__(self, master):
        self.master = master
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

        self.modalPane.bind("<Escape>", self._cancel)

        self.modalPane.title("Search GenBank")

        Label(self.modalPane, text="Find GenBank entries on NCBI server").pack(padx=5, pady=5)

        searchFrame = Frame(self.modalPane)
        searchFrame.pack(side=TOP)
        searchFrame.grid_columnconfigure(0, minsize=100)
        Label(searchFrame, text="Keyword: ").grid(row=0, column=0, pady=3)
        keywordentry = Entry(searchFrame, textvariable=self.keyword, width=25)
        keywordentry.grid(row=0, column=1, pady=3)
        Label(searchFrame, text="Organism: ").grid(row=1, column=0, pady=3)
        organismentry = Entry(searchFrame, textvariable=self.organism, width=25)
        organismentry.grid(row=1, column=1, pady=3)
        Label(searchFrame, text="Accession: ").grid(row=2, column=0, pady=3)
        accessionentry = Entry(searchFrame, textvariable=self.accession, width=25)
        accessionentry.grid(row=2, column=1, pady=3)
        searchButton = Button(searchFrame, text="Search", command=self._search)
        searchButton.grid(row=3, column=1, pady=3)

        listFrame = Frame(self.modalPane)
        listFrame.pack(side=TOP, padx=20, pady=20)

        scrollBar = Scrollbar(listFrame)
        scrollBar.pack(side=RIGHT, fill=Y)
        self.listBox = Listbox(listFrame, selectmode="multiple", width=100, height=20)
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

        if self.organism.get() == "" and self.keyword.get() == "" and self.accession.get() == "":
            showerror("Input error", "Please specify a search term first.")
            return
        # Build URL
        url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term="
        if self.organism.get() != "":
            url += "%22" + self.organism.get().replace(" ", "%20") + "%22[porgn]"
        if self.keyword.get() != "":
            if self.organism.get() != "":
                url += '%20AND%20' + "%22" + self.keyword.get().replace(" ", "%20") + "%22"
            else:
                url += "%22" + self.keyword.get().replace(" ", "%20") + "%22"
        if self.accession.get() != "":
            url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term=" + self.accession.get() + "[accession]"
        url += "&retmode=text&retmax=251&tool=multigeneblast&email=m.h.medema@rug.nl"
        # Search all IDs matching the query
        try:
            with urllib.request.urlopen(url) as req:
                xmltext = req.read().decode('utf-8')
                xmltext = xmltext.partition("<IdList>")[2].partition("</IdList>")[0]
                ids = [ID.partition("</Id>")[0] for ID in xmltext.split("<Id>")[1:]]
                if len(ids) == 0:
                    showerror('Error', 'No matches found for your search query.')
                    return
        except Exception:
            showerror('Error', 'Could not connect to NCBI server.\nPlease check your internet connection.')
            return
        # Show message if more than 250 found, asking user to further specify his search
        if len(ids) > 250:
            showwarning("Warning",
                        "More than 250 matches found for your search query.\nOnly top 250 matches will be shown.")
            ids = ids[:-1]
        self.ids = ids
        # Fetch titles for all IDs
        url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&id=" + ",".join(
            ids) + "&retmode=xml&retmax=251&tool=multigeneblast&email=m.h.medema@rug.nl"
        # Open 'loading...' message
        loading_frame = Toplevel(self.modalPane, height=200, width=400)
        loading_frame.title("Loading")
        message = "Loading. Please wait..."
        loading_lablel = Label(loading_frame, text=message)
        loading_lablel.grid(row=1, column=1, padx=50, pady=50)
        loading_frame.bind("<Escape>", lambda e: "return")
        self.modalPane.update()
        try:
            with urllib.request.urlopen(url) as req:
                xmltext = req.read().decode('utf-8')
                if "<DocSum>" not in xmltext or '<Item Name="Title" Type="String">' not in xmltext:
                    showerror('Error', 'Error in fetching match descriptions from NCBI server.')
                    return
                entries = xmltext.split("<DocSum>")[1:]
                descriptions = [entry.partition('<Item Name="Title" Type="String">')[2].partition('</Item>')[0]
                                for
                                entry in entries]
        except:
            showerror('Error', 'Could not connect to NCBI server.\nPlease check your internet connection.')
            return
        self._clear()
        loading_frame.destroy()
        for description in descriptions:
            self._insert(description)

    def _clear(self):
        self.listBox.delete(0, END)

    def _insert(self, string):
        self.listBox.insert(END, string)

    def _download(self, event=None):
        ids = [self.ids[indx] for indx in self.listBox.curselection()]
        outdir = tkinter.filedialog.askdirectory(mustexist=False)
        if outdir == "":
            return

        # label to show download progress
        loading_frame = Toplevel(self.modalPane, height=200, width=400)
        loading_frame.title("Download progress")
        message = StringVar()
        message.set("")
        loading_label = Label(loading_frame, textvariable=message)
        loading_label.grid(row=1, column=1, padx=50, pady=50)
        loading_frame.bind("<Escape>", lambda e: "return")
        loading_frame.protocol('WM_DELETE_WINDOW', lambda: cancel_download(req, loading_frame))

        failed_download_IDs = []
        for index, id in enumerate(ids):
            if not loading_frame.winfo_exists():
                return
            message.set("Downloading file {} out of {}.\nPlease wait...".format(index + 1, len(ids)))
            url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=" + id + \
                  "&rettype=gbwithparts&retmode=text&tool=multigeneblast&email=m.h.medema@rug.nl"

            self.modalPane.update()
            try:
                with urllib.request.urlopen(url) as req:
                    accumulated_size = 0
                    with open(outdir + os.sep + id + ".gb", 'wb') as fp:
                        while True:
                            try:
                                # close when the loading window is destroyed for wathever reason
                                if not loading_frame.winfo_exists():
                                    failed_download_IDs.append(id)
                                    break
                                self.modalPane.update()
                                chunk = req.read(CHUNK)
                                if not chunk:
                                    break
                                fp.write(chunk)
                            except:
                                failed_download_IDs.append(id)
                                break

                    # make sure to not spam request faster then allowed
                    time.sleep(0.5)
            except:
                failed_download_IDs.append(id)
                continue

        loading_frame.destroy()
        # Report download success
        showinfo("Download finished",
                 "Download completed successfully. Files can be found at {}.The following IDs could not be downloaded: {}".format(
                     outdir, ", ".join(failed_download_IDs)))

    def _cancel(self, event=None):
        self.modalPane.destroy()

class MakeDatabase(Frame):
  
    def __init__(self, master):
        self.master = master
        self.files = []
        self.dbname = StringVar()
        self.dbname.set("<Enter a name for your database>")
        self.__outdir_label = StringVar()
        self.__outdir_path = APPDATA
        self.initUI()

    def initUI(self):
        self.modal_pane = Toplevel(self.master, height=200, width=400)

        self.modal_pane.transient(self.master)
        self.modal_pane.grab_set()

        self.modal_pane.bind("<Escape>", self._cancel)

        self.modal_pane.title("Make MGB database from file")

        searchFrame = Frame(self.modal_pane)
        searchFrame.pack(side=TOP)
        searchFrame.grid_columnconfigure(0, minsize=100)
        Label(searchFrame, text="Database name: ").grid(row=0,column=0, pady=3)
        self.dbname_entry = Entry(searchFrame, textvariable=self.dbname, width=50)
        self.dbname_entry.grid(row=0,column=1, pady=3)
        self.dbname_entry.bind("<FocusOut>", self.OnValidate)
        self.dbname_entry.bind("<Return>", self.OnValidate)

        # TODO: see if this feature is wanted and realistic considereing the new database
        #text = Label(searchFrame, text="\nMake raw nucleotide database for tblastn-searching:")
        # text.grid(row=1,column=0, pady=5)
        # self.dbtype = CheckBox(searchFrame, [1,1], "")

        self.__outdir_label.set(self.__outdir_path)
        outdir_lbl = Label(searchFrame, text="Output folder name:")
        outdir_lbl.grid(row=1, column=0, sticky=N)
        outdir_text = Label(searchFrame, textvariable=self.__outdir_label)
        outdir_text.grid(row=1, column=1, sticky=N)
        outdir_button = Button(searchFrame, text="Select the output folder", command=self.select_out_directory)
        outdir_button.grid(row=1, column=2, sticky=N)


        searchButton = Button(searchFrame, text="Add files", command=self._addfile)
        searchButton.grid(row=2,column=0, pady=3)

        listFrame = Frame(self.modal_pane)
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

        clearButton = Button(self.modal_pane, text="Clear list", command=self._clear)
        clearButton.pack(side=LEFT)

        buttonFrame = Frame(self.modal_pane)
        buttonFrame.pack(side=BOTTOM)
        chooseButton = Button(buttonFrame, text="Make database", command=self.make_database)
        chooseButton.pack(side=LEFT)
        cancelButton = Button(buttonFrame, text="Cancel", command=self._cancel)
        cancelButton.pack(side=RIGHT)

    def select_out_directory(self):
        selected = tkinter.filedialog.askdirectory(mustexist=False)
        if selected == "":
            return
        try:
            with open(selected + os.sep + "test.txt", "w") as f:
                f.write("test")
            os.remove(selected + os.sep + "test.txt")
        except PermissionError:
            showerror("Error", "No permission to write to this folder. "
                "Please choose a directory in which you have writing permissions.")
            self.select_out_directory()
        display_selected = selected
        if len(selected) > 50:
            display_selected = "{}...".format(selected[:47])
        self.__outdir_label.set(display_selected)
        self.__outdir_path = selected.replace("/", os.sep)

    def OnValidate(self, val):
        for char in ILLEGAL_CHARACTERS:
            if self.dbname.get() != "<Enter a name for your database>":
                if char in self.dbname.get():
                    showerror("Input Error", "Database name contains illegal character. {} is not allowed".format(char))
                    self.dbname.set("<Enter a name for your database>")
                return
        if self.dbname.get() != None:
            if self.dbname.get() != "<Enter a name for your database>":
                dbname = self.dbname.get().replace(" ","_")
                self.dbname.set(dbname)
        else:
            self.dbname.set("<Enter a name for your database>")

    def _addfile(self):
        #returns a tuple of names
        filenames = tkinter.filedialog.askopenfilename(multiple=True, filetypes=(("GenBank files",('*.gbk','*.gb','*.genbank')),("EMBL files",('*.embl','*.emb'))), title="Select files")
        # else:
        #     filenames = tkinter.filedialog.askopenfilename(multiple=True, filetypes=(("FASTA files",('*.fasta','*.fa','*.fas')), ("GenBank files",('*.gbk','*.gb','*.genbank')),("EMBL files",('*.embl','*.emb'))), title="Select files", initialdir=self.directory)
        filenames = list(filenames)
        for filename in filenames:
            #when clicking cancel
            if filename == "":
                return
            root, ext = os.path.splitext(filename)
            if ext.lower() not in GENBANK_EXTENSIONS + EMBL_EXTENSIONS:
                showerror("Error", filename + " not added, as it does not have a GenBank or EMBL file extension.")
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

    def make_database(self):
        if self.dbname.get() + ".pal" in os.listdir(self.__outdir_path):
            answer = askyesno('Database name exists',
                              'A database with this name already exists. Overwrite?')
            if not answer:
                return
        command = self.__make_db_command()
        if command != None:
            self.__run_db_command(command)

    def __make_db_command(self):
        if len(self.files) == 0:
            showerror("Invalid input Error", "No files where submitted to create a database.")
            return None
        for char in self.dbname.get():
            if char in ILLEGAL_CHARACTERS:
                showerror("Input Error", "Database name contains illegal character. {} is not allowed".format(char))
                return None
        base_command = "{}{}make_database.py ".format(MGBPATH, os.sep)
        required_information = '-o "{}" -n "{}" -i "{}" '.format(self.__outdir_path, self.dbname.get(), " ".join(self.files))
        full_command = base_command + required_information
        return full_command

    def __run_db_command(self, command):
        #create a place to push messages to
        outbox = MessageBox(frame=self.master, title="Creating database...")
        self.modal_pane.update()
        popen = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        lines_iterator = iter(popen.stdout.readline, b"")
        while popen.poll() is None:
            for line in lines_iterator:
                outbox.text_insert(line)
                self.modal_pane.update()
        #check exit message
        output, error = popen.communicate()
        if error != b'':
            outbox.text_insert(error.decode('utf-8') + "\n")
            outbox.text_insert("MultiGeneBlast experienced a problem while creating the database. Please click"
                " the button below to send an error report, so we can solve the underlying problem.\n")
            outbox.change_errormessage(error)
            outbox.add_ok_button()
            outbox.add_error_button()
            self.modal_pane.update()
        else:

            outbox.text_insert("Database created.\nYou can now use this database "
                "by selecting '{}.pal' by clicking 'open database file' in the main window.\n".format(self.dbname.get()))
            outbox.add_ok_button()
            self.modal_pane.update()

    def _cancel(self, event=None):
        self.modal_pane.destroy()


class MakeOnlineDatabase(Frame):
  
    def __init__(self, master):
        self.master = master
        self.list = []
        self.selected = []
        self.keyword = StringVar()
        self.keyword.set("")
        self.organism = StringVar()
        self.organism.set("")
        self.accession = StringVar()
        self.accession.set("")
        self.dbname = StringVar()
        self.dbname.set("<Enter a name for your database>")
        self.descriptions = []

        self.__search_warning = StringVar()
        self.__outdir_label = StringVar()
        self.__outdir_label.set(APPDATA)
        self.__outdir_path = APPDATA

        self.initUI()

    def initUI(self):
        self.modal_pane = Toplevel(self.master)

        self.modal_pane.transient(self.master)
        self.modal_pane.grab_set()

        self.modal_pane.bind("<Escape>", self._cancel)

        self.modal_pane.title("Search GenBank")

        Label(self.modal_pane, text="Find GenBank entries on NCBI server").grid(row=0, column=0)
        
        searchFrame = Frame(self.modal_pane)
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
        searchButton.grid(row=3,column=1, pady=3)

        listFrame = Frame(self.modal_pane)
        listFrame.grid(row=2, column=0, padx=20, pady=5)
        Label(listFrame, text="Search results:").pack(side=TOP)
        Label(listFrame, textvariable=self.__search_warning, foreground="red").pack(side=TOP)
        listInnerFrame = Frame(listFrame)
        listInnerFrame.pack(side=BOTTOM)
        scrollBar = Scrollbar(listInnerFrame)
        scrollBar.pack(side=RIGHT, fill=Y)
        self.search_listbox = Listbox(listInnerFrame, selectmode="multiple", width=100, height=8)
        self.search_listbox.pack(side=LEFT, fill=Y)
        scrollBar.config(command=self.search_listbox.yview)
        self.search_listbox.config(yscrollcommand=scrollBar.set)
        self.list.sort()
        for item in self.list:
            self.search_listbox.insert(END, item)
        
        buttonFrame1 = Frame(self.modal_pane, width=100)
        buttonFrame1.grid(row=3, column=0)
        selectButton = Button(buttonFrame1, text="Select", command=self._insert_selected)
        selectButton.pack(side=RIGHT)

        Label(self.modal_pane, text=" ").grid(row=4, column=0, sticky=S)
        Label(self.modal_pane, text="Selected entries:").grid(row=5, column=0, sticky=S)
        
        listFrame2 = Frame(self.modal_pane)
        listFrame2.grid(row=6, column=0, padx=20, pady=5)
        scrollBar2 = Scrollbar(listFrame2)
        scrollBar2.pack(side=RIGHT, fill=Y)
        self.select_listbox = Listbox(listFrame2, selectmode="multiple", width=100, height=8)
        self.select_listbox.pack(side=LEFT, fill=Y)
        scrollBar2.config(command=self.select_listbox.yview)
        self.select_listbox.config(yscrollcommand=scrollBar2.set)
        self.select_listbox.bind("<Delete>", self._remove_selected)
        self.selected.sort()
        for item in self.selected:
            self.select_listbox.insert(END, item)

        dbnameFrame = Frame(self.modal_pane)
        dbnameFrame.grid(row=7,column=0, pady=5)
        Label(dbnameFrame, text="Database name: ").grid(row=0, column=0, padx=3, pady=5)
        dbname_entry = Entry(dbnameFrame, textvariable=self.dbname, width=40)
        dbname_entry.grid(row=0, column=1, padx=3, pady=5)
        dbname_entry.bind("<FocusOut>", self.OnValidate)
        dbname_entry.bind("<Return>", self.OnValidate)

        Label(dbnameFrame, text="Output folder name:").grid(row=1, column=0, padx=3, pady=5)
        outdir_label = Label(dbnameFrame, textvariable=self.__outdir_label, width=40)
        outdir_label.grid(row=1, column=1, padx=3, pady=5)
        Button(dbnameFrame, text="Select output folder", command=self.select_out_directory).grid(row=1, column=2, padx=3, pady=5)


        buttonFrame2 = Frame(self.modal_pane)
        buttonFrame2.grid(row=9, column=0)

        clearButton = Button(buttonFrame2, text="Clear selection", command=self._clear_selected)
        clearButton.grid(row=0, column=0, padx=3, pady=5)

        chooseButton = Button(buttonFrame2, text="Download and create database", command=self.download_make_database)
        chooseButton.grid(row=0, column=1, padx=3, pady=5)

        cancelButton = Button(buttonFrame2, text="Cancel", command=self._cancel)
        cancelButton.grid(row=0, column=2, padx=3, pady=5)

    def OnValidate(self, val):
        for char in ILLEGAL_CHARACTERS:
            if self.dbname.get() != "<Enter a name for your database>":
                if char in self.dbname.get():
                    showerror("Input Error", "Database name contains illegal character. {} is not allowed".format(char))
                    self.dbname.set("<Enter a name for your database>")
                return
        if self.dbname.get() != None:
            if self.dbname.get() != "<Enter a name for your database>":
                dbname = self.dbname.get().replace(" ","_")
                self.dbname.set(dbname)
        else:
            self.dbname.set("<Enter a name for your database>")
            
    def select_out_directory(self):
        selected = tkinter.filedialog.askdirectory(mustexist=False)
        if selected == "":
            return
        try:
            with open(selected + os.sep + "test.txt", "w") as f:
                f.write("test")
            os.remove(selected + os.sep + "test.txt")
        except PermissionError:
            showerror("Error", "No permission to write to this folder. "
                "Please choose a directory in which you have writing permissions.")
            self.select_out_directory()
        display_selected = selected
        if len(selected) > 50:
            display_selected = "{}...".format(selected[:47])
        self.__outdir_label.set(display_selected)
        self.__outdir_path = selected.replace("/", os.sep)

    def _search(self):
        if self.organism.get() == "" and self.keyword.get() == "" and self.accession.get() == "":
            showerror("Input error", "Please specify a search term first.")
            return
        # Build URL
        url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term="
        if self.organism.get() != "":
            url += "%22" + self.organism.get().replace(" ",
                                                       "%20") + "%22[porgn]"
        if self.keyword.get() != "":
            if self.organism.get() != "":
                url += '%20AND%20' + "%22" + self.keyword.get().replace(" ",
                                                                        "%20") + "%22"
            else:
                url += "%22" + self.keyword.get().replace(" ", "%20") + "%22"
        if self.accession.get() != "":
            url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term=" + self.accession.get() + "[accession]"
        url += "&retmode=text&retmax=251&tool=multigeneblast&email=m.h.medema@rug.nl"
        # Search all IDs matching the query
        try:
            with urllib.request.urlopen(url) as req:
                xmltext = req.read().decode('utf-8')
                xmltext = \
                xmltext.partition("<IdList>")[2].partition("</IdList>")[0]
                ids = [ID.partition("</Id>")[0] for ID in
                       xmltext.split("<Id>")[1:]]
                if len(ids) == 0:
                    showerror('Error',
                              'No matches found for your search query.')
                    return
        except Exception:
            showerror('Error','Could not connect to NCBI server.\nPlease check your internet connection.')
            return
        # Show message if more than 250 found, asking user to further specify his search
        if len(ids) > 250:
            self.__search_warning.set("Warning. More than 250 matches found for your search query. Only top 250 matches will be shown.")
            ids = ids[:-1]
        self.ids = ids
        # Fetch titles for all IDs
        url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&id=" + ",".join(
            ids) + "&retmode=xml&retmax=251&tool=multigeneblast&email=m.h.medema@rug.nl"
        # Open 'loading...' message
        loading_frame = Toplevel(self.modal_pane, height=200, width=400)
        loading_frame.title("Loading")
        loading_lablel = Label(loading_frame, text="Loading. Please wait...")
        loading_lablel.grid(row=1, column=1, padx=50, pady=50)
        loading_frame.bind("<Escape>", lambda e: "return")
        self.modal_pane.update()
        try:
            with urllib.request.urlopen(url) as req:
                xmltext = req.read().decode('utf-8')
                if "<DocSum>" not in xmltext or '<Item Name="Title" Type="String">' not in xmltext:
                    showerror('Error', 'Error in fetching match descriptions from NCBI server.')
                    return
                entries = xmltext.split("<DocSum>")[1:]
                descriptions = [entry.partition('<Item Name="Title" Type="String">')[2].partition('</Item>')[0] for entry in entries]
        except:
            showerror('Error', 'Could not connect to NCBI server.\nPlease check your internet connection.')
            return
        loading_frame.destroy()
        self._insert_search_results(descriptions)

    def _insert_search_results(self, descriptions):
        #make sure to delete contents before putting new in
        self.search_listbox.delete(0, END)
        self.descriptions = []
        for desc in descriptions:
            self.search_listbox.insert(END, desc)
            self.descriptions.append(desc)

    def _clear_selected(self):
        self.selected = []
        self.select_listbox.delete(0, END)

    def _remove_selected(self, idx):
        pos = self.select_listbox.curselection()
        self.select_listbox.delete(pos, pos)
        ID = self.ids[int(pos)]
        del self.selected[self.selected.index(ID)]

    def _insert_selected(self):
        selected = self.search_listbox.curselection()
        for sel_indx in selected:
            ID = self.ids[sel_indx]
            if ID not in self.selected:
                self.select_listbox.insert(END, self.descriptions[sel_indx])
                self.selected.append(ID)

    def download_make_database(self):
        if len(self.selected) == 0:
            showerror("Error", "No Entries selected.")
            return
        #Final check for correct database name
        for char in self.dbname.get():
            if char in ILLEGAL_CHARACTERS:
                showerror("Input Error","Database name contains illegal character. {} is not allowed".format(char))
                return
        # create a place to push messages to
        outbox = MessageBox(frame=self.master, title="Creating database...")
        #if downloading was a succes proceed with creating the database
        if self.download_database(outbox):
            command = self.__make_db_command()
            outbox.text_insert("Running make_database.py:\n")
            self.__run_db_command(command, outbox)

    def download_database(self, outbox):

        ids = self.selected

        outbox.text_insert("Downloading {} files. This can take a while..\n".format(len(ids)))
        self.modal_pane.update()
        url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=" + ",".join(ids) +\
              "&rettype=gbwithparts&retmode=text&tool=multigeneblast&email=m.h.medema@rug.nl"
        start_time = time.time()
        last_message_time = start_time
        try:
            with urllib.request.urlopen(url) as req:
                with open(self.__outdir_path + os.sep + self.dbname.get() + ".gbk", 'wb') as fp:
                    while True:
                        time_passed = time.time() - last_message_time
                        chunk = req.read(CHUNK)
                        if not chunk:
                            break
                        fp.write(chunk)
                        if time_passed > 10:
                            last_message_time = time.time()
                            outbox.text_insert("Downloading in progress. Time since start {}s\n".format(time.time() - start_time))
                            self.modal_pane.update()
        except Exception as error:
            outbox.text_insert("Failed to download files with the following error: {}\n".format(error))
            outbox.add_ok_button()
            outbox.add_error_button()
            return False

        # Report download success
        outbox.text_insert("Download completed successfully.\n")
        return True

    def __make_db_command(self):
        base_command = "{}{}make_database.py ".format(MGBPATH, os.sep)
        required_information = "-o {} -n {} -i {} ".format(self.__outdir_path, self.dbname.get(),self.__outdir_path + os.sep + self.dbname.get() + ".gbk")
        full_command = base_command + required_information
        return full_command

    def __run_db_command(self, command, outbox):
        # create a place to push messages to
        self.modal_pane.update()
        popen = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        lines_iterator = iter(popen.stdout.readline, b"")
        while popen.poll() is None:
            for line in lines_iterator:
                outbox.text_insert(line)
                self.modal_pane.update()
        # check exit message
        output, error = popen.communicate()
        if error != b'':
            outbox.text_insert(error.decode('utf-8') + "\n")
            outbox.text_insert("MultiGeneBlast experienced a problem while creating the database. Please click"
                               " the button below to send an error report, so we can solve the underlying problem.\n")
            outbox.change_errormessage(error)
            outbox.add_ok_button()
            outbox.add_error_button()
            self.modal_pane.update()
        else:
            outbox.text_insert("Database created. You can now use this database "
                "by selecting '{}.pal' by clicking 'open database file' in the main window.\n".format(self.dbname.get()))
            outbox.add_ok_button()
            self.modal_pane.update()
        
    def _cancel(self, event=None):
        self.modal_pane.destroy()