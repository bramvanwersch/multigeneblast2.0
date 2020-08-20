#!/usr/bin/env python
## Copyright (c) 2012 Marnix H. Medema
## Department of Microbial Physiology / Groningen Bioinformatics Centre
## University of Groningen
## License: GNU General Public License v3 or later
## A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

import sys
import os
from multiprocessing import freeze_support
import subprocess

#Find path to mgb files if run from another directory
pathfolders = os.environ['PATH'].split(os.pathsep)
pathfolders.append(os.getcwd())

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
from multigeneblast import main
import time
import urllib.request, urllib.error, urllib.parse
import tarfile
import traceback

from databases import GenbankFile, embl_to_genbank
from constants import *

MGBPATH = get_mgb_path()

def read_input_file_gui(path):
    """
    Read the input query file provided by the user in the GUI.

    :param path: an abosulte string path to the file
    :return: a list of genes and the size of the query
    """
    #the user clicked cancel
    if len(path) == 0:
        return [], 0
    root, ext = os.path.splitext(path)
    genes = []
    dna_seq_lenght = 0
    if ext.lower() in GENBANK_EXTENSIONS:
        try:
            gbf = GenbankFile(path=path)
            genes = list(gbf.proteins.keys())
            dna_seq_lenght = gbf.lenght
        except Exception as e:
            showerror("Reading file error", str(e))

    elif ext.lower() in EMBL_EXTENSIONS:
        try:
            gb_file_text = embl_to_genbank(path)
            gbf = GenbankFile(file_text=gb_file_text)
            genes = list(gbf.proteins.keys())
            dna_seq_lenght = gbf.lenght
        except Exception as e:
            showerror("Reading file error", str(e))
    return genes, dna_seq_lenght

def checkfasta(infile):
    """
    Check a fasta input file provided for the user when running architecture
    mode

    :param infile: an absolute path to a fasta file
    :return: a boolean telling if the file is valid or not
    """
    with open(infile,"r") as f:
        content = f.read()
    entries = [">" + entry for entry in content.split(">")][1:]
    if len(entries) < 2:
        showerror("Fasta read error","File Error","Please provide a FASTA file "
                "with multiple entries, representing the amino acid sequences of your query genes.")
        return False
    for entry in entries:
        entryseq = entry.partition("\n")[2]
        if (entryseq.count('A') + entryseq.count('a') + entryseq.count('C') + entryseq.count('c') + entryseq.count('G') + entryseq.count('g') + entryseq.count('T') + entryseq.count('t')) > (0.5 * len(entryseq)):
            showerror("Fasta read error","Nucleotide FASTA sequences provided. "
                    "Please provide a multi-FASTA file with amino acid sequences.")
            return False
    return True


def file_download(frame):
    global APPDATA
    GenBankFileDownload(frame, APPDATA)

def cancel_extraction(tar, toplevel):
    tar.close()
    toplevel.destroy()


def makedb_file(frame):
    global APPDATA
    MakeDatabase(frame)

def makedb_ncbi(frame):
    global APPDATA
    MakeOnlineDatabase(frame, APPDATA)

def makedb_gb():
    global APPDATA
    MakeGenBankDatabase(frame, APPDATA)


class MainMultiGeneBlastGui(Frame):
    """
    The main gui class that holds all the widgets for the main window that is
    represented when opening the GUI
    """
    def __init__(self, master):
        Frame.__init__(self, master)
        self.master = master

        self.main_window_frame = Frame(master)
        self.grid()

        #setup the menu widgets
        self.__setup_menu(master)

        self.__database_file_label = StringVar()
        self.database_file_path = ""

        self.__input_file_label = StringVar()
        self.input_file_path = ""

        self.__outdir_label = StringVar()
        self.outdir_path = ""

        self.searchtype = StringVar()

        #values for homology search when reading a query
        self.geneslist = []
        self.dnaseqlength = 0

        #innitiate the widgets of the main window
        self.__innitialize_widgets()

    def __setup_menu(self, master):
        """
        Setup the menu at the top of the window

        :param master: the tk() componenet controlling the frame
        """
        menu = Menu(master)
        master.config(menu=menu)

        filemenu = Menu(menu)
        menu.add_cascade(label="File", menu=filemenu)
        # filemenu.add_command(label="Open input file", command=self.file_open)
        # filemenu.add_command(label="Select __database_file_label", command=db_open)
        filemenu.add_command(label="Exit", command=master.quit)

        downloadmenu = Menu(menu)
        menu.add_cascade(label="Download", menu=downloadmenu)
        downloadmenu.add_command(label="Download GenBank entry",
                                 command=lambda: file_download(frame))
        # downloadmenu.add_command(label="Download MGB GenBank database",
        #                          command=db_download)

        dbmenu = Menu(menu)
        menu.add_cascade(label="Database", menu=dbmenu)
        dbmenu.add_command(label="Create database from files",
                           command=lambda: makedb_file(self))
        dbmenu.add_command(label="Create database from online GenBank entries",
                           command=lambda: makedb_ncbi(frame))
        dbmenu.add_command(label="Create database from GenBank subdivisions",
                           command=makedb_gb)

        helpmenu = Menu(menu)
        menu.add_cascade(label="Help", menu=helpmenu)
        helpmenu.add_command(label="About...", command=about)

    def __innitialize_widgets(self):
        """
        Innitialize all widgets present in the innitial frame
        :return:
        """
        #frame to hold the widgets
        self.main_window_frame.grid()

        #minimum sizes
        self.main_window_frame.grid_columnconfigure(0, minsize=300)
        self.main_window_frame.grid_columnconfigure(1, minsize=300)
        self.main_window_frame.grid_rowconfigure(8, minsize=32)
        self.main_window_frame.grid_rowconfigure(9, minsize=32)
        self.main_window_frame.grid_rowconfigure(10, minsize=34)
        self.main_window_frame.grid_rowconfigure(11, minsize=34)
        self.main_window_frame.grid_rowconfigure(12, minsize=34)
        self.main_window_frame.grid_rowconfigure(13, minsize=34)
        self.main_window_frame.grid_rowconfigure(14, minsize=34)
        self.main_window_frame.grid_rowconfigure(15, minsize=34)
        self.main_window_frame.grid_rowconfigure(16, minsize=34)

        #header text before
        header_text = Label(self.main_window_frame, text="MultiGeneBlast input\n")
        header_text.grid(row=0,column=0)

        #Database selection
        if "genbank_mf.pal" in os.listdir("."):
            self.__database_file_label.set((os.getcwd() + os.sep + "genbank_mf").replace("\\", "/"))
        else:
            self.__database_file_label.set("<No database selected>")
        database_label = Label(self.main_window_frame, text = "Current database:")
        database_label.grid(row=1, column=0)
        database_text = Label(self.main_window_frame, textvariable = self.__database_file_label)
        database_text.grid(row=1, column=1)

        database_button = Button(self.main_window_frame, text = "Open database file", command=self.db_open)
        database_button.grid(row=1,column=2)

        #Input file selection
        self.__input_file_label.set("<No file selected>")
        infile_label = Label(self.main_window_frame, text = "Current input file: ")
        infile_label.grid(row=2,column=0)
        infile_text = Label(self.main_window_frame, textvariable = self.__input_file_label)
        infile_text.grid(row=2,column=1)

        infile_button = Button(self.main_window_frame, text = "Open input file", command=self.file_open)
        infile_button.grid(row=2,column=2)

        # Output folder selection
        self.__outdir_label.set("<No output directory selected>")
        outdir_lbl = Label(self.main_window_frame, text="Output folder name:")
        outdir_lbl.grid(row=3, column=0, sticky=N)
        outdir_text = Label(self.main_window_frame, textvariable=self.__outdir_label)
        outdir_text.grid(row=3, column=1, sticky=N)
        outdir_button = Button(self.main_window_frame, text="Select the output folder", command=self.select_out_directory)
        outdir_button.grid(row=3, column=2, sticky=N)

        self.searchtype.set("homology")

        # Search type radio buttons
        space1 = Label(self.main_window_frame, text="\n")
        space1.grid(row=4, column=0, columnspan=2)
        homradio = Radiobutton(self.main_window_frame, text="Homology search", variable=self.searchtype,
                               command=self.setup_search_type_widgets, value="homology")
        homradio.grid(row=5, column=0)
        archradio = Radiobutton(self.main_window_frame, text="Architecture search", variable=self.searchtype,
                                command=self.setup_search_type_widgets, value="architecture")
        archradio.grid(row=5, column=1)
        homradio.select()

        space2 = Label(self.main_window_frame, text="\n")
        space2.grid(row=6, column=0, columnspan=2)

        # Start nt selection
        self.nuc_start_lbl = Label(self.main_window_frame, text="Nucleotide start position of fragment to search:")
        self.nuc_start_lbl.grid(row=7, column=0, sticky=N)
        self.cstart = ScaleBar(self.main_window_frame, self, [7, 1], 0, self.dnaseqlength, 0,
                               resetgenes="y")

        # End nt selection
        self.nuc_end_lbl = Label(self.main_window_frame, text="Nucleotide end position of fragment to search:")
        self.nuc_end_lbl.grid(row=8, column=0)
        self.cend = ScaleBar(self.main_window_frame, self, [8, 1], 0, self.dnaseqlength,
                             self.dnaseqlength, resetgenes="y")

        # Genes selection
        self.gene_selection_lbl = Label(self.main_window_frame,
                                        text="\nOr genes to search (locus tags / accession numbers):")
        self.gene_selection_lbl.grid(row=9, column=0, sticky=N)
        self.selectedgenes = GeneSelectionFrame(self.main_window_frame, self, [9, 1],
                                                self.geneslist)

        # Number of cores to use
        nr_cores_lbl = Label(self.main_window_frame, text="\nNumber of CPU cores to be used:")
        nr_cores_lbl.grid(row=10, column=0, sticky=N)
        self.cores = ScaleBar(self.main_window_frame, self, [10, 1], 1,
                              determine_cpu_nr("all"), determine_cpu_nr("all"))

        # Minimal sequence coverage of Blast hit
        hitspergene_lbl = Label(self.main_window_frame,
                      text="\nNumber of Blast hits per gene to be mapped:")
        hitspergene_lbl.grid(row=11, column=0, sticky=N)
        self.hitspergene = ScaleBar(self.main_window_frame, self, [11, 1], 50, 1000, 250)

        # Weight of synteny conservation score
        syntentwheight_lbl = Label(self.main_window_frame,
                      text="\nWeight of synteny conservation in hit sorting:")
        syntentwheight_lbl.grid(row=12, column=0, sticky=N)
        self.syntenyweight = ScaleBar(self.main_window_frame, self, [12, 1], 0.0, 2.0, 0.5,
                                      input_type="float")

        # Minimal sequence coverage of Blast hit
        seqcov_lbl = Label(self.main_window_frame,
                      text="\nMinimal sequence coverage of BLAST hits:")
        seqcov_lbl.grid(row=13, column=0, sticky=N)
        self.seqcov = ScaleBar(self.main_window_frame, self, [13, 1], 0, 100, 25)

        # Minimal % ID of Blast hit
        percid_lbl = Label(self.main_window_frame, text="\nMinimal % identity of BLAST hits:")
        percid_lbl.grid(row=14, column=0, sticky=N)
        self.percid = ScaleBar(self.main_window_frame, self, [14, 1], 0, 100, 30)

        # Maximum allowed distance between hits in locus
        distancekb_lbl = Label(self.main_window_frame,
                       text="\nMaximum distance between genes in locus (kb):")
        distancekb_lbl.grid(row=15, column=0, sticky=N)
        self.distancekb = ScaleBar(self.main_window_frame, self, [15, 1], 1, 100, 20)

        # Number of hit loci to show
        hitsperloci_lbl = Label(self.main_window_frame, text="\nNumber of hit loci to display:")
        hitsperloci_lbl.grid(row=16, column=0, sticky=N)
        self.hitsperloci = SpinBox(self.main_window_frame, [16, 1], 0, 500, 50, 250)

        # Generate Muscle alignments
        muscle_lbl = Label(self.main_window_frame,
                       text="\nMuscle alignment of homologs with queries:")
        muscle_lbl.grid(row=17, column=0, sticky=N)
        self.muscle = CheckBox(self.main_window_frame, [17, 1], "")

        # Empty space
        space3 = Label(self.main_window_frame, text="   ")
        space3.grid(row=18, column=0, columnspan=2)

        #Run button to start analysis
        run_mgb_btn = Button(self.main_window_frame, text="Run MultiGeneBlast", command=self.run_multigeneblast)
        run_mgb_btn.grid(row=20,column=0)

    def setup_search_type_widgets(self):
        """
        Called when clicking one of the radio buttons. Adds some widgets in the
        case of homollogy mode and removes them for architecture mode
        """
        if self.searchtype.get() == "homology":
            self.__input_file_label.set("<No file selected>")
            if len(list(self.nuc_start_lbl.grid_info().keys())) == 0:
                self.nuc_start_lbl.grid(row=7, column=0)
                self.nuc_end_lbl.grid(row=8, column=0)
                self.gene_selection_lbl.grid(row=9, column=0)
                self.cstart = ScaleBar(self.main_window_frame, self, [7, 1], 0,
                                       self.dnaseqlength, 0, resetgenes="y")
                self.cend = ScaleBar(self.main_window_frame, self, [8, 1], 0,
                                     self.dnaseqlength, self.dnaseqlength,
                                     resetgenes="y")
                self.selectedgenes = GeneSelectionFrame(self.main_window_frame, self,
                                                        [9, 1], self.geneslist)
        elif self.searchtype.get() == "architecture":
            self.__input_file_label.set("<No file selected>")
            if len(list(self.nuc_start_lbl.grid_info().keys())) > 0:
                self.nuc_start_lbl.grid_forget()
                self.nuc_end_lbl.grid_forget()
                self.gene_selection_lbl.grid_forget()
                self.cstart.grid_remove()
                self.cend.grid_remove()
                self.selectedgenes.grid_remove()
                self.dnaseqlength = 0
                self.geneslist = []

    def update_starts_ends(self):
        """
        adjust the scalebars for the start and end location of the query
        selection as well as the genes that can be selected.
        """
        if self.searchtype.get() == "homology":
            self.cstart.grid_remove()
            self.cend.grid_remove()
            self.selectedgenes.grid_remove()
            self.cstart = ScaleBar(self.main_window_frame, self, [6,1], 0, self.dnaseqlength, 0, resetgenes="y")
            self.cend = ScaleBar(self.main_window_frame, self, [7,1], 0, self.dnaseqlength, self.dnaseqlength, resetgenes="y")
            self.selectedgenes = GeneSelectionFrame(self.main_window_frame, self, [8,1], self.geneslist)

    def file_open(self):
        """
        Select an input query file as well as load genes and the lenght of the
        file into the GUI
        """
        if self.searchtype.get() == "homology":
            location = tkinter.filedialog.askopenfilename(filetypes=(("GenBank files", ('*.gbk', '*.gb', '*.genbank')),("EMBL files", ('*.embl', '*.emb'))))

            #when canceling
            if location == "":
                return

            genes, dna_seq_length = read_input_file_gui(location)
            if dna_seq_length > 0 and len(genes) > 0:
                self.geneslist = genes
                self.dnaseqlength = dna_seq_length
                self.update_starts_ends()
            else:
                location = "<No file selected>"
        else:
            location = tkinter.filedialog.askopenfilename(filetypes=(("FASTA files", ('*.fasta', '*.fas', '*.fa')),
                                                                     ("All file types", ('*.*'))))
            # when canceling
            if location == "":
                return
            if not checkfasta(location):
                location = "<No file selected>"

        #set the location
        display_location = location
        if len(location) > 50:
            display_location = "{}...".format(location[:47])
        self.__input_file_label.set(display_location)
        self.input_file_path = location.replace("/", os.sep)

    def db_open(self):
        """
        Select a database and check if all files are present
        """
        location = tkinter.filedialog.askopenfilename(filetypes=(("MGB database files",('*.pal','*.nal')),("MGB database",())))
        if location == "":
            return
        root, ext = os.path.splitext(location)
        to_path, db_file = os.path.split(location)
        dbname = db_file.split(".")[0]

        # make sure the db folder contains all files required
        db_folder = os.listdir(to_path)
        expected_files = [dbname + ext for ext in DATABASE_EXTENSIONS]
        for file in expected_files:
            if file not in db_folder:
                showerror("Missing database file", "The following file {} is missing"
                        " for the database with alias file {}".format(file, db_file))
                location = "<No database selected>"

        display_location = location
        if len(location) > 50:
            display_location = "{}...".format(location[:47])
        self.__database_file_label.set(display_location)
        self.database_file_path = location.replace("/", os.sep)

    def select_out_directory(self):
        """
        Select the output directory. If files are present a new directory is
        created in that directory to prevent unwanted overwriting or removal of
        files.
        """
        selected = tkinter.filedialog.askdirectory(mustexist=False)
        if selected == "":
            return
        try:
            dir_files = os.listdir(selected)
        except FileNotFoundError:
            showerror("Error", "Folder does not exist. Select a different one")
            self.select_out_directory()
        #if files in the directory make a directory in the directory
        if len(dir_files) > 0:
            #generate a unique name
            count = 1
            name = "MultiGeneBlast_out{}".format(count)
            while name in dir_files:
                count += 1
                name = "MultiGeneBlast_out{}".format(count)
            try:
                os.mkdir(selected + os.sep + name)
            except:
                showerror("Error", "No permission to write to this folder. "
                    "Please choose a directory in which you have writing permissions.")
                self.select_out_directory()
            selected = selected + "/" + name
        #if not files are present test if you can write in the folder
        else:
            #easier to ask forgiveniss then permission
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
        self.outdir_path = selected.replace("/", os.sep)


    def run_multigeneblast(self):
        """
        Called when clicking the button to run multigeneblast
        """
        command = self.__construct_mgb_command()
        if command != None:
            self.__run_mgb_command(command)


    def __construct_mgb_command(self):
        """
        Construct a command from the data input by the user
        :return:
        """
        #first check if in database and out are defined
        if self.input_file_path == "" or self.database_file_path == "" or  self.outdir_path == "":
            showerror("Input incomplete Error", "Please fill in the input query file, the database file and the output directory")
            return None

        #if there are genes selected use genes, otherwise use locations
        filter_options = ""
        if self.searchtype.get() == "homology":
            genes = self.selectedgenes.getval().split(";")
            if len(genes) > 1:
                for gene in genes:
                    if gene not in self.geneslist:
                        showerror("Input Error", "Invalid gene name {} in the provided list of genes".format(gene))
                        return None
                filter_options = "-g {} ".format(" ".join(genes))
            else:
                start = self.cstart.getval()
                end = self.cend.getval()
                if end <= start:
                    showerror("Input Error", "Start location has to be smaller then end location")
                    return None
                filter_options = "-f {} -t {} ".format(start, end)

        base_command = "{}{}multigeneblast.py ".format(MGBPATH, os.sep)

        #change some values to adhere to the correct format of the command line tool
        if self.muscle.getval() == 0:
            muscle_val = "n"
        else:
            muscle_val = "y"
        nr_pages = self.hitsperloci.getval() / HITS_PER_PAGE
        if self.hitsperloci.getval() % HITS_PER_PAGE == 0:
            nr_pages = int(nr_pages)
        else:
            nr_pages = int(nr_pages) + 1

        general_options = "-in {} -db {} -out {} -c {} -hpg {} -msc {} -dkb {} -mpi {} -m {} -sw {} -op {}".format(self.input_file_path,
                                                                                                                   self.database_file_path, self.outdir_path, self.cores.getval(), self.hitspergene.getval(), self.seqcov.getval(),
                                                                                                                   self.distancekb.getval(), self.percid.getval(), muscle_val, self.syntenyweight.getval(), nr_pages)
        full_command = base_command + filter_options + general_options
        return full_command

    def __run_mgb_command(self, command):
        """
        Run multigeneblast using subprocess, while printing the sdtout to a
        MessageBox object

        :param command: a string that represents a valid mgb command
        """
        #create a place to push messages to
        outbox = MessageBox(frame=self, title="Running MultiGeneBlast...")
        self.update()
        popen = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        lines_iterator = iter(popen.stdout.readline, b"")
        while popen.poll() is None:
            for line in lines_iterator:
                outbox.text_insert(line)
                self.update()
        #check exit message
        output, error = popen.communicate()
        if error != b'':
            outbox.text_insert(str(error))
            outbox.text_insert("MultiGeneBlast experienced a problem. Please click"
                " the button below to send an error report, so we can solve the underlying problem.\n")
            outbox.change_errormessage(error)
            outbox.add_ok_button()
            outbox.add_error_button()
            self.update()
        else:
            self.open_final_webpage()

            outbox.text_insert('\nVisual output can be accessed by opening: '
                '"file://' + self.outdir_path + os.sep + 'visual' + os.sep +
                               'displaypage1.xhtml" with a web browser\n')
            outbox.add_ok_button()
            self.update()

    def open_final_webpage(self):
        """
        Try to open an xhtml file in a browser
        """
        display_page_1 = self.outdir_path + os.sep + 'visual' + os.sep + 'displaypage1.xhtml'
        if sys.platform == ('win32'):
            os.startfile(display_page_1)
        else:
            try:
                firefox = webbrowser.get('firefox')
                firefox.open_new_tab("file://" + display_page_1)
            except:
                pass
            try:
                safari = webbrowser.get('safari')
                safari.open_new_tab("file://" + display_page_1)
            except:
                pass
            try:
                chrome = webbrowser.get('/usr/bin/google-chrome %s')
                chrome.open_new_tab("file://" + display_page_1)
            except:
                pass
            try:
                webbrowser.open_new_tab("file://" + display_page_1)
            except:
                pass


def maingui():

    #to ensure that the MGB directory can be located
    # os.environ['PATH'] = os.getcwd() + os.pathsep + os.environ['PATH']

    root = Tk()
    if sys.platform == ('win32'):
        try:
            root.iconbitmap(default='mgb.ico')
        except:
            pass
    root.title('MultiGeneBlast')
    root.geometry("%dx%d%+d%+d" % (850, 750, 0, 0))

    mg = MainMultiGeneBlastGui(root)


if __name__ == '__main__':
    freeze_support()
    maingui()
    mainloop()
