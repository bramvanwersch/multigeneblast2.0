#!/usr/bin/env python3

"""
Classes used by several classes from the gui.

Original creator: Marnix Medena
Recent contributor: Bram van Wersch

Copyright (c) 2012 Marnix H. Medema
License: GNU General Public License v3 or later
A copy of GNU GPL v3 should have been included in this software package in LICENSE.txt.
"""


import time
import tkinter.filedialog
from tkinter.messagebox import askyesno, showinfo
from tkinter import SINGLE, S, W
import urllib.request
import urllib.error
import urllib.parse

from utilities import ILLEGAL_CHARACTERS
from constants import APPDATA, GENBANK_EXTENSIONS, EMBL_EXTENSIONS, get_mgb_path, CHUNK, FASTA_EXTENSIONS
from gui_utility import *

MGBPATH = get_mgb_path()


class GeneSelectionFrame(Frame):
    """
    Widget combination of button and entry to allow the user to select a collection of genes from an input query file
    """
    def __init__(self, master, mgb_gui_i, entrieslist):
        """
        :param master: the master object
        :param mgb_gui_i: multigeneblast_gui instance
        :param entrieslist:
        """
        Frame.__init__(self, master)
        self.SpObj = mgb_gui_i
        self.entrieslist = entrieslist
        self.selected = StringVar()

        self.init_widgets()

    def init_widgets(self):
        """
        Innitialize the widgets
        """
        self.selected.set("<Select genes>")
        genes_entry = Entry(self, textvariable=self.selected, width=25)
        genes_entry.grid(row=0, column=0)
        genes_entry.bind("<FocusOut>", self.on_validate)
        genes_entry.bind("<Return>", self.on_validate)
        selection_box = Button(self, text="Select genes", command=self.select)
        selection_box.grid(row=0, column=1, padx=10)
      
    def select(self):
        """
        Open a ListBoxChoice object to select a number of genes
        """
        selectedgenes = ListBoxChoice(self, "Select genes", "Select query genes:", self.entrieslist).return_value()
        if selectedgenes is not None:
            self.SpObj.cstart.setval(0)
            self.SpObj.cend.setval(0)
            self.selected.set(selectedgenes)

    def set_selectable_genes(self, genes):
        """
        Change the self.entrielist property

        :param genes: a list of strings
        """
        self.entrieslist = genes
        self.clear_selection()
        
    def on_validate(self, *args):
        """
        Called when finished text to the entry widget

        :param args: optional arguments
        """
        if self.selected.get() is not None and self.selected.get() != "<Select genes>" and\
                len(self.selected.get().split(";")) > 1:
            selectedgenes = self.selected.get()
            self.SpObj.cstart.setval(0)
            self.SpObj.cend.setval(0)
            self.selected.set(selectedgenes)
        else:
            self.selected.set("<Select genes>")

    def clear_selection(self):
        """
        Clear the text in the entry widget
        """
        self.selected.set("<Select genes>")

    def getval(self):
        """
        Get the text in the entry widget

        :return: String text in the entry widget
        """
        return self.selected.get()


class SearchFrame(Frame):
    """
    Implements a Text widget and some entries to search the NCBI database for genbank files
    """
    def __init__(self, master):
        """
        :param master: Tk or Frame object that is the master window of this TopLevel
        """
        super().__init__(master)
        self.descriptions = []
        self.keyword = StringVar()
        self.keyword.set("")
        self.organism = StringVar()
        self.organism.set("")
        self.accession = StringVar()
        self.accession.set("")

        self.__search_warning = StringVar()
        self.search_listbox = None
        self.init_widgets()

    def init_widgets(self):
        """
        Innitialize the widgets
        """
        search_frame = Frame(self.master)
        search_frame.grid(row=0)
        Label(search_frame, text="Keyword: ").grid(row=0, column=0, pady=3)
        keywordentry = Entry(search_frame, textvariable=self.keyword, width=25)
        keywordentry.grid(row=0, column=1, pady=3)
        Label(search_frame, text="Organism: ").grid(row=1, column=0, pady=3)
        organismentry = Entry(search_frame, textvariable=self.organism, width=25)
        organismentry.grid(row=1, column=1, pady=3)
        Label(search_frame, text="Accession: ").grid(row=2, column=0, pady=3)
        accessionentry = Entry(search_frame, textvariable=self.accession, width=25)
        accessionentry.grid(row=2, column=1, pady=3)
        search_button = Button(search_frame, text="Search", command=self._search)
        search_button.grid(row=3, column=1, pady=3)

        list_frame = Frame(self.master)
        list_frame.grid(row=1, padx=20, pady=20)
        Label(list_frame, textvariable=self.__search_warning, foreground="red").pack(side=TOP)

        scroll_bar = Scrollbar(list_frame)
        scroll_bar.pack(side=RIGHT, fill=Y)
        self.search_listbox = Listbox(list_frame, selectmode="multiple", width=100, height=12)
        self.search_listbox.pack(side=LEFT, fill=Y)
        scroll_bar.config(command=self.search_listbox.yview)
        self.search_listbox.config(yscrollcommand=scroll_bar.set)

    def _search(self):
        """
        Function for searching the NCBI database for genbank files using species, keywords and organisms
        """
        self.__search_warning.set("")
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
            url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term=" + self.accession.get()\
                  + "[accession]"
        url += "&retmode=text&retmax=251&tool=multigeneblast&email=m.h.medema@rug.nl"
        # Search all IDs matching the query
        try:
            with urllib.request.urlopen(url) as req:
                xmltext = req.read().decode('utf-8')
                xmltext = xmltext.partition("<IdList>")[2].partition("</IdList>")[0]
                ids = [ID.partition("</Id>")[0] for ID in
                       xmltext.split("<Id>")[1:]]
                if len(ids) == 0:
                    showerror('Error',
                              'No matches found for your search query.')
                    return
        except Exception:
            showerror('Error', 'Could not connect to NCBI server.\nPlease check your internet connection.')
            return
        # Show message if more than 250 found, asking user to further specify his search
        if len(ids) > 250:
            self.__search_warning.set("Warning. More than 250 matches found for your search query."
                                      " Only top 250 matches will be shown.")
            ids = ids[:-1]
        self.ids = ids
        # Fetch titles for all IDs
        url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&id=" + ",".join(
            ids) + "&retmode=xml&retmax=251&tool=multigeneblast&email=m.h.medema@rug.nl"
        # Open 'loading...' message
        loading_frame = Toplevel(self.master, height=200, width=400)
        loading_frame.title("Loading")
        loading_lablel = Label(loading_frame, text="Loading. Please wait...")
        loading_lablel.grid(row=1, column=1, padx=50, pady=50)
        loading_frame.bind("<Escape>", lambda e: "return")
        self.master.update()
        try:
            with urllib.request.urlopen(url) as req:
                xmltext = req.read().decode('utf-8')
                if "<DocSum>" not in xmltext or '<Item Name="Title" Type="String">' not in xmltext:
                    showerror('Error', 'Error in fetching match descriptions from NCBI server.')
                    return
                entries = xmltext.split("<DocSum>")[1:]
                descriptions = [entry.partition('<Item Name="Title" Type="String">')[2].partition('</Item>')[0] for
                                entry in entries]
        except Exception:
            showerror('Error', 'Could not connect to NCBI server.\nPlease check your internet connection.')
            return
        loading_frame.destroy()
        self._insert_search_results(descriptions)

    def _clear(self):
        """
        Clear the list_box with search results
        """
        self.search_listbox.delete(0, END)

    def _insert_search_results(self, descriptions):
        """
        Insert a list of descriptions that are search results

        :param descriptions: a list of strings
        """
        # make sure to delete contents before putting new in
        self.search_listbox.delete(0, END)
        self.descriptions = []
        for desc in descriptions:
            self.search_listbox.insert(END, desc)
            self.descriptions.append(desc)

    def get_selected_ids(self):
        """
        Get the IDs of all search results that are selected in the list_box

        :return: a list of string IDs
        """
        return [self.ids[indx] for indx in self.search_listbox.curselection()]

    def get_selected_descriptions(self):
        """
        Get the descriptions of all search results that are selected in the list_box

        :return: a list of string descriptions
        """
        return [self.descriptions[indx] for indx in self.search_listbox.curselection()]


class GenBankFileDownload(Toplevel):
    """
    TopLevel window for downloading one or more genbank files
    """

    def __init__(self, master):
        """
        :param master: Tk or Frame object that is the master window of this TopLevel
        """
        super().__init__(master)
        self.grid()
        self.transient(self.master)
        self.grab_set()

        self.bind("<Escape>", self._cancel)
        self.title("Search GenBank")

        self.search_frame = None
        self.init_widgets()

    def init_widgets(self):
        """
        innitialize the widgets
        """
        self.search_frame = SearchFrame(self)
        self.search_frame.grid(row=0)

        button_frame = Frame(self)
        button_frame.grid(row=3)

        choose_button = Button(button_frame, text="Download", command=self.one_by_one_download)
        choose_button.pack(side=LEFT)

        cancel_button = Button(button_frame, text="Cancel", command=self._cancel)
        cancel_button.pack(side=RIGHT)

    def one_by_one_download(self):
        """
        Download all selected genbank files into separate files.
        """
        ids = self.search_frame.get_selected_ids()
        outdir = tkinter.filedialog.askdirectory(mustexist=False)
        if outdir == "":
            return

        # label to show download progress
        loading_frame = Toplevel(self, height=200, width=400)
        loading_frame.title("Download progress")
        message = StringVar()
        message.set("")
        loading_label = Label(loading_frame, textvariable=message)
        loading_label.grid(row=1, column=1, padx=50, pady=50)
        loading_frame.bind("<Escape>", lambda e: "return")

        failed_download_ids = []
        for index, id_ in enumerate(ids):
            if not loading_frame.winfo_exists():
                return
            message.set("Downloading file {} out of {}.\nPlease wait...".format(index + 1, len(ids)))
            url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=" + id_ + \
                  "&rettype=gbwithparts&retmode=text&tool=multigeneblast&email=m.h.medema@rug.nl"

            self.update()
            try:
                with urllib.request.urlopen(url) as req:
                    with open(outdir + os.sep + id_ + ".gb", 'wb') as fp:
                        while True:
                            try:
                                # close when the loading window is destroyed for wathever reason
                                if not loading_frame.winfo_exists():
                                    failed_download_ids.append(id_)
                                    break
                                self.update()
                                chunk = req.read(CHUNK)
                                if not chunk:
                                    break
                                fp.write(chunk)
                            except Exception:
                                failed_download_ids.append(id_)
                                break

                    # make sure to not spam request faster then allowed
                    time.sleep(0.5)
            except Exception:
                failed_download_ids.append(id_)
                continue

        loading_frame.destroy()
        # Report download success
        message = "Download completed successfully. Files can be found at {}".format(outdir)
        if len(failed_download_ids) > 0:
            message += "\n\nWARNING: The following IDs could not be downloaded: {}"\
                .format(", ".join(failed_download_ids))
        showinfo("Download finished", message)

    def _cancel(self):
        """
        Called when escape is pressed to kill the toplevel window
        """
        self.destroy()


class MakeDatabase(Toplevel):
    """
    TopLevel for creating a database from local genbank files
    """
    def __init__(self, master):
        """
        :param master: Tk or Frame object that is the master window of this TopLevel
        """
        super().__init__(master)
        self.grid()

        self.files = []
        self.dbname = StringVar()
        self.dbname.set("<Enter a name for your database>")
        self.__outdir_label = StringVar()
        self.__outdir_path = APPDATA

        self.dbname_entry = None
        self.dbtype = None
        self.file_list_box = None
        self.init_widgets()

    def init_widgets(self):
        """
        Innitialize the widgets
        """
        self.transient(self.master)
        self.grab_set()
        self.bind("<Escape>", self._cancel)
        self.title("Make MGB database from file")

        search_frame = Frame(self)
        search_frame.grid(row=0, padx=20)
        Label(search_frame, text="Database name: ").grid(row=0, column=0, padx=5, sticky=W)
        self.dbname_entry = Entry(search_frame, textvariable=self.dbname, width=50)
        self.dbname_entry.grid(row=0, column=1, padx=5)
        self.dbname_entry.bind("<FocusOut>", self.on_validate)
        self.dbname_entry.bind("<Return>", self.on_validate)

        self.__outdir_label.set(self.__outdir_path)
        outdir_lbl = Label(search_frame, text="Output folder name:")
        outdir_lbl.grid(row=1, column=0, sticky=W, padx=5)
        outdir_text = Label(search_frame, textvariable=self.__outdir_label)
        outdir_text.grid(row=1, column=1, padx=5)
        outdir_button = Button(search_frame, text="Select the output folder", command=self.select_out_directory)
        outdir_button.grid(row=1, column=2, sticky=W)

        text = Label(search_frame, text="Make raw nucleotide database for tblastn-searching:")
        text.grid(row=2, column=0, padx=5, sticky=W)
        self.dbtype = CheckBox(search_frame, "")
        self.dbtype.grid(row=2, column=1)

        search_button = Button(search_frame, text="Add files", command=self._add_files)
        search_button.grid(row=3, column=1)

        list_frame = Frame(self)
        list_frame.grid(row=1, pady=10, padx=20)
        scroll_bar = Scrollbar(list_frame)
        scroll_bar.pack(side=RIGHT, fill=Y)
        self.file_list_box = Listbox(list_frame, selectmode=SINGLE, width=120, height=10)
        self.file_list_box.pack(side=LEFT, fill=Y)
        scroll_bar.config(command=self.file_list_box.yview)
        self.file_list_box.config(yscrollcommand=scroll_bar.set)
        self.file_list_box.bind("<Delete>", self._remove)

        button_frame = Frame(self)
        button_frame.grid(row=2, pady=10)
        clear_button = Button(button_frame, text="Clear all", command=self._clear)
        clear_button.grid(column=0, row=0, padx=5)
        clear_selected_button = Button(button_frame, text="Clear selected", command=self._remove)
        clear_selected_button.grid(column=1, row=0, padx=5)
        make_database_button = Button(button_frame, text="Make database", command=self.make_database)
        make_database_button.grid(column=2, row=0, padx=5)
        cancel_button = Button(button_frame, text="Cancel", command=self._cancel)
        cancel_button.grid(column=3, row=0, padx=5)

    def select_out_directory(self):
        """
        Select the output directory.
        """
        path = select_out_directory()
        if path is None:
            return
        display_selected = path
        if len(path) > 50:
            display_selected = "{}...".format(path[:47])
        self.__outdir_label.set(display_selected)
        self.__outdir_path = path.replace("/", os.sep)

    def on_validate(self, *args):
        """
        Called when leaving the entry widget

        :param args: optional arguments
        """
        for char in ILLEGAL_CHARACTERS:
            if self.dbname.get() != "<Enter a name for your database>":
                if char in self.dbname.get():
                    showerror("Input Error", "Database name contains illegal character. {} is not allowed".format(char))
                    self.dbname.set("<Enter a name for your database>")
                    return

    def _add_files(self):
        """
        Add one or more files to the self.file_list_box ListBox object
        """
        # returns a tuple of names
        dbtype = self.dbtype.getval()
        if dbtype == 0:
            filenames = tkinter.filedialog.askopenfilename(multiple=True, filetypes=(("GenBank files",
                            ('*.gbk', '*.gb', '*.genbank')), ("EMBL files", ('*.embl', '*.emb'))), title="Select files")
        else:
            filenames = tkinter.filedialog.askopenfilename(multiple=True, filetypes=(("FASTA files",
                    ('*.fasta', '*.fa', '*.fas')), ("GenBank files", ('*.gbk', '*.gb', '*.genbank')),
                    ("EMBL files", ('*.embl', '*.emb'))), title="Select files")
        filenames = list(filenames)
        if dbtype == 0:
            allowed_extensions = GENBANK_EXTENSIONS + EMBL_EXTENSIONS
        else:
            allowed_extensions = GENBANK_EXTENSIONS + EMBL_EXTENSIONS + FASTA_EXTENSIONS
        for filename in filenames:
            # when clicking cancel
            if filename == "":
                return
            root, ext = os.path.splitext(filename)
            if ext.lower() not in allowed_extensions:
                showerror("Error", filename + " not added, as it does not have a GenBank or EMBL file extension.")
            if filename not in self.files:
                self.files.append(filename)
                self._insert(filename)

    def _clear(self):
        """
        Clear the the complete file_list_box
        """
        self.file_list_box.delete(0, END)
        self.files = []

    def _insert(self, item):
        """
        Insert a file into the file_list_box

        :param item: string file_path
        :return:
        """
        self.file_list_box.insert(END, item)

    def _remove(self):
        """
        Remove a single entry from the list box
        """
        pos = self.file_list_box.curselection()
        if pos:
            del self.files[int(pos[0])]
            self.file_list_box.delete(pos, pos)

    def make_database(self):
        """
        Function called when the make_database_button is pressed
        """
        if self.dbname.get() + ".dmnd" in os.listdir(self.__outdir_path) or\
                self.dbname.get() + ".nal" in os.listdir(self.__outdir_path):
            answer = askyesno('Database name exists',
                              'A database with this name already exists. Overwrite?')
            if not answer:
                return
        command = self.__make_db_command()
        if command is not None:
            self.__run_db_command(command)

    def __make_db_command(self):
        """
        Constructs a command line command to run make_database.py

        :return: None or the command as a string
        """
        if len(self.files) == 0:
            showerror("Invalid input Error", "No files where submitted to create a database.")
            return None
        for char in self.dbname.get():
            if char in ILLEGAL_CHARACTERS:
                showerror("Input Error", "Database name contains illegal character. {} is not allowed".format(char))
                return None
        base_command = "{}{}make_database.py ".format(MGBPATH, os.sep)
        if self.dbtype.getval() == 0:
            dbtype = "prot"
        else:
            dbtype = "nucl"
        # make sure to place '"' around inputs that can contain spaces
        required_information = '-o "{}" -n "{}" -i {} -t {}'.format(self.__outdir_path, self.dbname.get(),
                                                                    " ".join(['"{}"'.format(f) for f in self.files]),
                                                                    dbtype)
        full_command = base_command + required_information
        return full_command

    def __run_db_command(self, command):
        """
        Run the command on the command line using the run_external_command function. This function simply
        adds some additional information depending on the finishing error code

        :param command: a command as string
        """
        # create a place to push messages to
        outbox = MessageBox(master=self.master, title="Creating database...")
        exit_code, expected = run_extrenal_command(command, outbox, self)
        if exit_code == 0:
            if self.dbtype.getval() == 0:
                outbox.text_insert("Database created. You can now use this database by selecting '{}.dmnd' by clicking"
                                   " 'open database file' in the main window.\n".format(self.dbname.get()))
            else:
                outbox.text_insert("Database created. You can now use this database by selecting '{}.nal' by clicking "
                                   "'open database file' in the main window.\n".format(self.dbname.get()))
        elif not expected:
            outbox.text_insert("MultiGeneBlast experienced a problem while creating the database. Please click the"
                               " button below to send an error report, so we can solve the underlying problem.\n")

    def _cancel(self):
        """
        Called when escape is pressed to kill the toplevel window
        """
        self.destroy()


class MakeOnlineDatabase(Toplevel):
  
    def __init__(self, master):
        """
        :param master: Tk or Frame object that is the master window of this TopLevel
        """
        super().__init__(master)
        self.transient(self.master)
        self.grab_set()

        self.bind("<Escape>", self._cancel)

        self.title("Search GenBank")

        self.selected = []
        self.dbname = StringVar()
        self.dbname.set("<Enter a name for your database>")
        self.search_frame = None
        self.dbtype = None

        self.__outdir_label = StringVar()
        self.__outdir_label.set(APPDATA)
        self.__outdir_path = APPDATA

        self.select_listbox = None
        self.init_widgets()

    def init_widgets(self):
        """
        Innitialize the widgets
        """
        Label(self, text="Find GenBank entries on NCBI server").grid(row=0, column=0)

        self.search_frame = SearchFrame(self)
        self.search_frame.grid(row=1)

        button_frame1 = Frame(self, width=100)
        button_frame1.grid(row=3, column=0)
        select_button = Button(button_frame1, text="Select", command=self._insert_selected)
        select_button.pack(side=RIGHT)

        Label(self, text=" ").grid(row=4, column=0, sticky=S)
        Label(self, text="Selected entries:").grid(row=5, column=0, sticky=S)
        
        list_frame2 = Frame(self)
        list_frame2.grid(row=6, column=0, padx=20, pady=5)
        scroll_bar2 = Scrollbar(list_frame2)
        scroll_bar2.pack(side=RIGHT, fill=Y)
        self.select_listbox = Listbox(list_frame2, selectmode="multiple", width=100, height=12)
        self.select_listbox.pack(side=LEFT, fill=Y)
        scroll_bar2.config(command=self.select_listbox.yview)
        self.select_listbox.config(yscrollcommand=scroll_bar2.set)
        self.select_listbox.bind("<Delete>", self._remove_selected)

        dbname_frame = Frame(self)
        dbname_frame.grid(row=7, column=0, pady=5)
        Label(dbname_frame, text="Database name: ").grid(row=0, column=0, padx=3, pady=5)
        dbname_entry = Entry(dbname_frame, textvariable=self.dbname, width=40)
        dbname_entry.grid(row=0, column=1, padx=3, pady=5)
        dbname_entry.bind("<FocusOut>", self.on_validate)
        dbname_entry.bind("<Return>", self.on_validate)

        Label(dbname_frame, text="Output folder name:").grid(row=1, column=0, padx=3, pady=5)
        outdir_label = Label(dbname_frame, textvariable=self.__outdir_label, width=40)
        outdir_label.grid(row=1, column=1, padx=3, pady=5)
        Button(dbname_frame, text="Select output folder", command=self.select_out_directory).grid(row=1, column=2,
                                                                                                  padx=3, pady=5)
        text = Label(dbname_frame, text="Make raw nucleotide database for tblastn-searching:")
        text.grid(row=2, column=1, padx=5, sticky=W)
        self.dbtype = CheckBox(dbname_frame, "")
        self.dbtype.grid(row=2, column=2, padx=5, sticky=W)

        button_frame2 = Frame(self)
        button_frame2.grid(row=9, column=0)

        clear_all_button = Button(button_frame2, text="Clear all", command=self._clear_all)
        clear_all_button.grid(row=0, column=0, padx=3, pady=5)

        clear_selected_button = Button(button_frame2, text="Clear selected", command=self._remove_selected)
        clear_selected_button.grid(row=0, column=1, padx=3, pady=5)

        choose_button = Button(button_frame2, text="Download and create database", command=self.download_make_database)
        choose_button.grid(row=0, column=2, padx=3, pady=5)

        cancel_button = Button(button_frame2, text="Cancel", command=self._cancel)
        cancel_button.grid(row=0, column=3, padx=3, pady=5)

    def on_validate(self, *args):
        """
        Called when leaving the entry widgets

        :param args: optional arguments
        """
        for char in ILLEGAL_CHARACTERS:
            if self.dbname.get() != "<Enter a name for your database>":
                if char in self.dbname.get():
                    showerror("Input Error", "Database name contains illegal character. {} is not allowed".format(char))
                    self.dbname.set("<Enter a name for your database>")
                    return

    def select_out_directory(self):
        """
        Select the output directory.
        """
        path = select_out_directory()
        if path is None:
            return
        display_selected = path
        if len(path) > 50:
            display_selected = "{}...".format(path[:47])
        self.__outdir_label.set(display_selected)
        self.__outdir_path = path.replace("/", os.sep)

    def _clear_all(self):
        """
        Clear all values in the select_listbox
        """
        self.selected = []
        self.select_listbox.delete(0, END)

    def _remove_selected(self):
        """
        Remove the specifically selected entry in the select_listbox
        """
        selected_positions = list(self.select_listbox.curselection())
        selected_positions.sort(key=int, reverse=True)
        for pos in selected_positions:
            self.select_listbox.delete(pos, pos)
            del self.selected[pos]

    def _insert_selected(self):
        """
        Insert the selected lines from the searchframe listbox into the select_listbox
        """
        descriptions = self.search_frame.get_selected_descriptions()
        ids = self.search_frame.get_selected_ids()
        for index, id_ in enumerate(ids):
            if id_ not in self.selected:
                self.select_listbox.insert(END,  descriptions[index])
                self.selected.append(id_)

    def download_make_database(self):
        """
        Function called when the button download and create database is pressed
        """
        if len(self.selected) == 0:
            showerror("Input Error", "No Entries selected.")
            return
        # Final check for correct database name
        for char in self.dbname.get():
            if char in ILLEGAL_CHARACTERS:
                showerror("Input Error", "Database name contains illegal character. {} is not allowed".format(char))
                return
        # create a place to push messages to
        outbox = MessageBox(master=self.master, title="Creating database...")
        # if downloading was a succes proceed with creating the database
        if self.download_database(outbox):
            command = self.__make_db_command()
            outbox.text_insert("Running make_database.py:\n")
            self.__run_db_command(command, outbox)

    def download_database(self, outbox):
        """
        Download all the entries in the select_listbox into one genbank file

        :return a boolean that tells if the download was a succes or not
        """
        ids = self.selected

        outbox.text_insert("Downloading {} files. This can take a while..\n".format(len(ids)))
        self.update()
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
                            self.update()
        except Exception as error:
            outbox.text_insert("Failed to download files with the following error: {}\n".format(error))
            outbox.add_ok_button()
            outbox.add_report_button()
            return False

        # Report download success
        outbox.text_insert("Download completed successfully.\n")
        return True

    def __make_db_command(self):
        """
        Create the command to run make_database.py to create a database from the downloaded genbank file

        :return: a string that is the full command to run make_database.py
        """
        base_command = "{}{}make_database.py ".format(MGBPATH, os.sep)
        if self.dbtype.getval() == 0:
            dbtype = "prot"
        else:
            dbtype = "nucl"
        required_information = "-o {} -n {} -i {} -t {}".format(self.__outdir_path, self.dbname.get(),
                                                                self.__outdir_path + os.sep + self.dbname.get() +
                                                                ".gbk", dbtype)
        full_command = base_command + required_information
        return full_command

    def __run_db_command(self, command, outbox):
        """
        Run the database making command. Use a MessageBox object to notify the user what is happening

        :param command: a string that is a command to run make_database.py
        :param outbox: a MessageBox where messages are pushed to giving feedback to the user
        """
        # create a place to push messages to
        self.update()
        exit_code, expected = run_extrenal_command(command, outbox, self.master)
        if exit_code == 0:
            if self.dbtype.getval() == 0:
                outbox.text_insert("Database created. You can now use this database by selecting '{}.dmnd' by clicking"
                                   " 'open database file' in the main window.\n".format(self.dbname.get()))
            else:
                outbox.text_insert("Database created. You can now use this database by selecting '{}.nal' by clicking "
                                   "'open database file' in the main window.\n".format(self.dbname.get()))
        elif not expected:
            outbox.text_insert("MultiGeneBlast experienced a problem while creating the database. Please click"
                               " the button below to send an error report, so we can solve the underlying problem.\n")
        
    def _cancel(self):
        """
        Called when escape the window is close
        """
        self.destroy()
