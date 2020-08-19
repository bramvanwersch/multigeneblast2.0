#!/usr/bin/env python
## Copyright (c) 2012 Marnix H. Medema
## Department of Microbial Physiology / Groningen Bioinformatics Centre
## University of Groningen
## License: GNU General Public License v3 or later
## A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.

import sys
import os
from multiprocessing import freeze_support

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
from multigeneblast import *
import time
import urllib.request, urllib.error, urllib.parse
import tarfile
import traceback
from constants import *


def read_input_file_gui(path):
    #TODO move these files to more appropriate place
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
        except MultiGeneBlastException as e:
            showerror("Reading file error", str(e))

    elif ext.lower() in EMBL_EXTENSIONS:
        try:
            gb_file_text = embl_to_genbank(path)
            gbf = GenbankFile(file_text=gb_file_text)
            genes = list(gbf.proteins.keys())
            dna_seq_lenght = gbf.lenght
        except MultiGeneBlastException as e:
            showerror("Reading file error", str(e))
    return genes, dna_seq_lenght

def checkfasta(infile):
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

            #Step 4: Run BLAST on genbank_mf __database_file_label
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
                outbox.text_insert("Error (" + str(sys.exc_info()[0]).rpartition(".")[2][:-2] + ") while loading __database_file_label:\n")
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
    #    command = 'python multigeneblast.py -in "' + __input_file_label + '" -genes ' + genes + ' -db "' + __database_file_label + '"'
    #  else:
    #    command = 'python multigeneblast.py -in "' + __input_file_label + '" -from ' + cstart + " -to " + cend + ' -db "' + __database_file_label + '"'
    #else:
    #  command = 'python multigeneblast.py -in "' + __input_file_label + '"' + ' -db "' + __database_file_label + '"'
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
    #__input_file_label = parse_absolute_paths(__input_file_label)
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


    #frame.update_idletasks()


def file_download(frame):
    global APPDATA
    GenBankFileDownload(frame, APPDATA)

def cancel_extraction(tar, toplevel):
    tar.close()
    toplevel.destroy()

def db_download():
    file_url = 'http://gbic.webhosting.rug.nl/genbank_mf.tar.gz'
    #file_url = 'https://downloads.sourceforge.net/project/multigeneblast/db/genbank_mf_test.rar?use_mirror=switch'
    answer = askyesno('Confirmation', 'Database downloading will take a while, and the __database_file_label will occuppy ~15 GB disk space.\n Are you sure?')
    if not answer:
        return
    global APPDATA
    currentdir = os.getcwd()
    os.chdir(APPDATA)
    #Check if __database_file_label is not already present
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
    #Download the MGB GenBank __database_file_label
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
    #Extract the __database_file_label TAR/GZ file
    extracting = Toplevel(frame, height=200, width=400)
    extracting.title("Database extraction")
    message = "Extracting __database_file_label.\nPlease wait..."
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
        showerror("Error","Error extracting __database_file_label. Please try to extract it manually.")
        extracting.destroy()
        os.chdir(currentdir)
        return
    extracting.destroy()
    os.chdir(currentdir)
    showinfo("Extraction finished", "Database extraction finished.\nYou can now use this __database_file_label by selecting 'genbank_mf.pal' under 'Select __database_file_label' in the 'File' menu.")

def makedb_file(frame):
    global APPDATA
    MakeDatabase(frame, APPDATA)

def makedb_ncbi(frame):
    global APPDATA
    MakeOnlineDatabase(frame, APPDATA)

def makedb_gb():
    global APPDATA
    MakeGenBankDatabase(frame, APPDATA)


class MainMultiGeneBlastGui(Frame):
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
        self.searchtype = StringVar()

        #innitiate the widgets of the main window
        self.__innitialize_widgets()

    def __setup_menu(self, master):
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
        downloadmenu.add_command(label="Download MGB GenBank database",
                                 command=db_download)

        dbmenu = Menu(menu)
        menu.add_cascade(label="Database", menu=dbmenu)
        dbmenu.add_command(label="Create database from files",
                           command=lambda: makedb_file(frame))
        dbmenu.add_command(label="Create database from online GenBank entries",
                           command=lambda: makedb_ncbi(frame))
        dbmenu.add_command(label="Create database from GenBank subdivisions",
                           command=makedb_gb)

        helpmenu = Menu(menu)
        menu.add_cascade(label="Help", menu=helpmenu)
        helpmenu.add_command(label="About...", command=about)

    def __innitialize_widgets(self):
        self.main_window_frame.grid()

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
            self.__database_file_label.set("<Nodatabase selected>")
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

        self.searchtype.set("homology")
        self.geneslist = []
        self.dnaseqlength = 0


        # Search type radio buttons
        space1 = Label(self.main_window_frame, text="\n")
        space1.grid(row=3, column=0, columnspan=2)
        homradio = Radiobutton(self.main_window_frame, text="Homology search", variable=self.searchtype,
                               command=self.setup_search_type_widgets, value="homology")
        homradio.grid(row=4, column=0)
        archradio = Radiobutton(self.main_window_frame, text="Architecture search", variable=self.searchtype,
                                command=self.setup_search_type_widgets, value="architecture")
        archradio.grid(row=4, column=1)
        homradio.select()

        space2 = Label(self.main_window_frame, text="\n")
        space2.grid(row=5, column=0, columnspan=2)

        # Start nt selection
        self.text2 = Label(self.main_window_frame,
                           text="Nucleotide start position of fragment to search:")
        self.text2.grid(row=6, column=0, sticky=N)
        self.cstart = ScaleBar(self.main_window_frame, self, [6, 1], 0, self.dnaseqlength, 0,
                               resetgenes="y")

        # End nt selection
        self.text3 = Label(self.main_window_frame,
                           text="Nucleotide end position of fragment to search:")
        self.text3.grid(row=7, column=0)
        self.cend = ScaleBar(self.main_window_frame, self, [7, 1], 0, self.dnaseqlength,
                             self.dnaseqlength, resetgenes="y")

        # Genes selection
        self.text4 = Label(self.main_window_frame,
                           text="\nOr genes to search (locus tags / accession numbers):")
        self.text4.grid(row=8, column=0, sticky=N)
        self.selectedgenes = GeneSelectionFrame(self.main_window_frame, self, [8, 1],
                                                self.geneslist)

        # Output folder name
        text5 = Label(self.main_window_frame, text="\nOutput folder name:")
        text5.grid(row=9, column=0, sticky=N)
        self.outputfolder = OutputFolderSelectionFrame(self.main_window_frame, [9, 1], USERDIR)

        # Number of cores to use
        text6 = Label(self.main_window_frame, text="\nNumber of CPU cores to be used:")
        text6.grid(row=10, column=0, sticky=N)
        self.cores = ScaleBar(self.main_window_frame, self, [10, 1], 1,
                              determine_cpu_nr("all"), determine_cpu_nr("all"))

        # Minimal sequence coverage of Blast hit
        text7 = Label(self.main_window_frame,
                      text="\nNumber of Blast hits per gene to be mapped:")
        text7.grid(row=11, column=0, sticky=N)
        self.hitspergene = ScaleBar(self.main_window_frame, self, [11, 1], 50, 1000, 250)

        # Weight of synteny conservation score
        text8 = Label(self.main_window_frame,
                      text="\nWeight of synteny conservation in hit sorting:")
        text8.grid(row=12, column=0, sticky=N)
        self.syntenyweight = ScaleBar(self.main_window_frame, self, [12, 1], 0.0, 2.0, 0.5,
                                      input_type="float")

        # Minimal sequence coverage of Blast hit
        text9 = Label(self.main_window_frame,
                      text="\nMinimal sequence coverage of BLAST hits:")
        text9.grid(row=13, column=0, sticky=N)
        self.seqcov = ScaleBar(self.main_window_frame, self, [13, 1], 0, 100, 25)

        # Minimal % ID of Blast hit
        text10 = Label(self.main_window_frame, text="\nMinimal % identity of BLAST hits:")
        text10.grid(row=14, column=0, sticky=N)
        self.pid = ScaleBar(self.main_window_frame, self, [14, 1], 0, 100, 30)

        # Maximum allowed distance between hits in locus
        text11 = Label(self.main_window_frame,
                       text="\nMaximum distance between genes in locus (kb):")
        text11.grid(row=15, column=0, sticky=N)
        self.distancekb = ScaleBar(self.main_window_frame, self, [15, 1], 1, 100, 20)

        # Number of hit loci to show
        text12 = Label(self.main_window_frame, text="\nNumber of hit loci to display:")
        text12.grid(row=16, column=0, sticky=N)
        self.pages = SpinBox(self.main_window_frame, [16, 1], 0, 500, 50, 250)

        # Generate Muscle alignments
        text13 = Label(self.main_window_frame,
                       text="\nMuscle alignment of homologs with queries:")
        text13.grid(row=17, column=0, sticky=N)
        self.muscle = CheckBox(self.main_window_frame, [17, 1], "")

        # Empty space
        space3 = Label(self.main_window_frame, text="   ")
        space3.grid(row=18, column=0, columnspan=2)


        #Run button to start analysis
        button = Button(self.main_window_frame, text="Run MultiGeneBlast",) #command=lambda : runmgb(frame, __input_file_label.get(), str(SearchPrefs.cstart.getval()), str(SearchPrefs.cend.getval()), SearchPrefs.selectedgenes.getval(), SearchPrefs.outputfolder.getfolder(), __database_file_label.get(), SearchPrefs.cores.getval(), SearchPrefs.seqcov.getval(), SearchPrefs.pid.getval(), SearchPrefs.distancekb.getval(), SearchPrefs.pages.getval(), SearchPrefs.muscle.getval(), SearchPrefs.getsearchtype(), SearchPrefs.hitspergene.getval(), SearchPrefs.syntenyweight.getval()))
        button.grid(row=20,column=0)

    def setup_search_type_widgets(self):
        if self.searchtype.get() == "homology":
            self.__input_file_label.set("<No file selected>")
            if len(list(self.text2.grid_info().keys())) == 0:
                self.text2.grid(row=6, column=0)
                self.text3.grid(row=7, column=0)
                self.text4.grid(row=8, column=0)
                self.cstart = ScaleBar(self.main_window_frame, self, [6, 1], 0,
                                       self.dnaseqlength, 0, resetgenes="y")
                self.cend = ScaleBar(self.main_window_frame, self, [7, 1], 0,
                                     self.dnaseqlength, self.dnaseqlength,
                                     resetgenes="y")
                self.selectedgenes = GeneSelectionFrame(self.main_window_frame, self,
                                                        [8, 1], self.geneslist)
        elif self.searchtype.get() == "architecture":
            self.__input_file_label.set("<No file selected>")
            if len(list(self.text2.grid_info().keys())) > 0:
                self.text2.grid_forget()
                self.text3.grid_forget()
                self.text4.grid_forget()
                self.cstart.grid_remove()
                self.cend.grid_remove()
                self.selectedgenes.grid_remove()
                self.dnaseqlength = 0
                self.geneslist = []

    def update_starts_ends(self):
        if self.searchtype.get() != "architecture":
            self.cstart.grid_remove()
            self.cend.grid_remove()
            self.selectedgenes.grid_remove()
            self.cstart = ScaleBar(self.main_window_frame, self, [6,1], 0, self.dnaseqlength, 0, resetgenes="y")
            self.cend = ScaleBar(self.main_window_frame, self, [7,1], 0, self.dnaseqlength, self.dnaseqlength, resetgenes="y")
            self.selectedgenes = GeneSelectionFrame(self.main_window_frame, self, [8,1], self.geneslist)

    def file_open(self):
        """
        Open the input query file
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
        self.input_file_path = location

    def db_open(self):
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
        self.database_file_path = location



def maingui():

    # os.environ['PYTHON'] = os.path.dirname(os.path.abspath(__file__)) + os.sep + "python"
    # os.environ['EXEC'] = os.path.dirname(os.path.abspath(__file__)) + os.sep + "exec"
    # os.environ['PATH'] = os.environ['EXEC'] + os.pathsep + os.environ['PYTHON'] + os.pathsep + os.environ['PATH']

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
