#!/usr/bin/env python3

"""
Classes for disecting genbank and embl files and putting them in a database format that can be used by multigeneblast

Original creator: Marnix Medena
Recent contributor: Bram van Wersch
"""

# imports
import tarfile
import pickle
from collections import OrderedDict
from abc import ABC, abstractmethod

import urllib.request
import urllib.error
import urllib.parse

from urllib.error import URLError
import http.client

from utilities import *
from constants import GENBANK_EXTENSIONS, EMBL_EXTENSIONS, TEMP, FASTA_EXTENSIONS


def clean_dna_sequence(dna_seq):
    """
    Clean a DNA sequence that is formatted like;

          1 atgcggcggg aaatt etc.

    :param dna_seq: a DNA sequence in the afformentioned format as a string
    :return: a DNA sequence that only contains nucleotides
    """
    dna_seq = dna_seq.replace(" ", "")
    dna_seq = dna_seq.replace("\t", "")
    dna_seq = dna_seq.replace("\n", "")
    dna_seq = dna_seq.replace("0", "")
    dna_seq = dna_seq.replace("1", "")
    dna_seq = dna_seq.replace("2", "")
    dna_seq = dna_seq.replace("3", "")
    dna_seq = dna_seq.replace("4", "")
    dna_seq = dna_seq.replace("5", "")
    dna_seq = dna_seq.replace("6", "")
    dna_seq = dna_seq.replace("7", "")
    dna_seq = dna_seq.replace("8", "")
    dna_seq = dna_seq.replace("9", "")
    dna_seq = dna_seq.replace("/", "")
    return dna_seq


def embl_to_genbank(embl_filepath):
    """
    Convert an embl file in such a way that the GenbankFile object can read it
    and convert it appropriately

    :param embl_filepath: a path to an embl file.
    :return: a string that can be read by a GenbankFile object
    """
    logging.debug("Converting to make embl file {} readable for the GenbankFile object".format(embl_filepath))
    try:
        with open(embl_filepath, "r") as f:
            file_text = f.read()
    except Exception:
        logging.critical("Invalid embl file {}. Exiting...".format(embl_filepath))
        raise MultiGeneBlastException("Invalid embl file {}.".format(embl_filepath))

    # make sure to remove potential old occurances of \r. Acts like a \n
    file_text = file_text.replace("\r", "\n")

    # do a basic check to see if the embl file is valid
    if "FT   CDS " not in file_text or ("\nSQ" not in file_text):
        logging.critical("Embl file {} is not properly formatted or contains no sequences. Exiting..."
                         .format(embl_filepath))
        raise MultiGeneBlastException("Embl file {} is not properly formatted or contains no sequences"
                                      .format(embl_filepath))
    text_lines = file_text.split("\n")

    # change certain line starts
    line_count = 0
    while line_count < len(text_lines):
        line = text_lines[line_count]
        if line.startswith("FT"):
            text_lines[line_count] = line.replace("FT", "  ", 1)
        if line.startswith("SQ"):
            text_lines[line_count] = "ORIGIN"
        elif line.startswith("AC"):
            text_lines[line_count] = line.replace("AC   ", "ACCESSION   ", 1)
        elif line.startswith("DE"):
            # change all the definition lines to make them readable
            text_lines[line_count] = line.replace("DE   ", "DEFINITION  ", 1)
            line_count += 1
            def_line = text_lines[line_count]
            while not def_line.startswith("XX"):
                text_lines[line_count] = def_line.replace("DE   ", "            ", 1)
                line_count += 1
                def_line = text_lines[line_count]
        line_count += 1
    logging.debug("Embl file {} made readable for genbank file object.".format(embl_filepath))
    return "\n".join(text_lines)


class Database(ABC):
    """
    Abstraction level for tracking multiple genbank and embl files and writing
    the output into a format that can be used by MultiGeneBlast
    """
    def __init__(self, base_path, paths):
        logging.info("Started creating database...")
        # care this is a generator object not a list
        self._files = self._read_files(base_path, paths)
        if len(self._files) == 0:
            logging.critical("Failed to load any of the provided database files.")
            raise MultiGeneBlastException("Failed to load any of the provided database files.")

    @abstractmethod
    def _read_files(self, base_path, paths):
        """
        Method for reading input files, each subclass of Database should have
        this method

        :param base_path: an absolute path to the directory make_database.py was
        executed from
        :param paths: a absolute or relative path to all the database files from
        the make_database.py directory
        :return: a list
        """
        return []

    @abstractmethod
    def get_fasta(self):
        """
        Return a fasta of the files in the database. This method is required
        to make Blast+ databases

        :return: a string in fasta format
        """
        return ""

    def create(self, outdir, dbname):
        """
        Create the files for the database using the self._files property

        :param outdir: an output directory for the database
        :param dbname: the name of all files in the output directory
        """
        self._write_files(outdir, dbname)
        self.__create_tar_file(outdir, dbname)

    @abstractmethod
    def _write_files(self, outdir, dbname):
        """
        Write the contigs in self._files as pickled objects.

        :param outdir: an output directory for the database
        :param dbname: the name of all files in the output directory
        """
        pass

    def __create_tar_file(self, outdir, dbname):
        """
        Create gzipped tar file to save all the pickled contigs into.

        :param outdir: an output directory for the database
        :param dbname: the name of all files in the output directory
        """
        logging.info("Writing to tar file...")
        with tarfile.open("{}{}{}_contigs.tar.gz".format(outdir, os.sep, dbname), "w:gz") as tar:
            tar.add("{}{}pickles".format(TEMP, os.sep), arcname=os.path.basename(dbname))

    def _convert_wgs_master_record(self, text, file_name):
        """
        Convert WGS master record file into multiple files containing contigs

        :param text: the file text of the WGS file
        :param file_name: the name of the WGS file, it ensures that the records
        get read into unique files
        :return: a list of new gb files that need to be parsed in addition to
        the existing list
        """
        # If seq_record is a WGS master record, parse out contig accession numbers and download these as
        # separate seq_records
        contig_ranges = None
        if "WGS_SCAFLD  " in text:
            contig_ranges = text.split("WGS_SCAFLD  ")[1].split("\n")[0]
        elif "WGS         " in text:
            contig_ranges = text.split("WGS         ")[1].split("\n")[0]

        # failsafe
        if contig_ranges is None:
            return []

        # extract all the ranges
        if ";" in contig_ranges:
            ranges = []
            for contig_range in contig_ranges.split(";"):
                if "-" in contig_range:
                    ranges.append(contig_range.split("-"))
                else:
                    ranges.append([contig_range])
            contig_ranges = ranges
        elif "-" in contig_ranges:
            contig_ranges = [contig_ranges.split("-")]
        else:
            contig_ranges = [[contig_ranges]]

        # extract all contigs for each range
        allcontigs = []
        for contig_range in contig_ranges:
            if len(contig_range) == 1:
                allcontigs.extend(contig_range)
                continue
            startnumber, endnumber = '', ''
            alpha_tag = ''
            for char in contig_range[0].partition(".")[0]:
                if char.isdigit():
                    startnumber += char
                else:
                    alpha_tag += char
            for char in contig_range[1].partition(".")[0]:
                if char.isdigit():
                    endnumber += char
            nrzeros = 0
            # certain entries rely on a number of zeroes. This can be at the start of a number
            for char in startnumber:
                if char == "0":
                    nrzeros += 1
                else:
                    break
            contig_range = [alpha_tag + nrzeros * "0" + str(number) for number in
                            range(int(startnumber), int(endnumber))]
            allcontigs.extend(contig_range)

        # Create contig groups of 50 (reasonable download size per download)
        contig_groups = [allcontigs[i:i+50] for i in range(0, len(allcontigs), 50)]

        # Download contigs and parse into seq_record objects
        new_files = self.__fetch_urls(contig_groups, '&rettype=gbwithparts&retmode=text', file_name)
        return new_files

    def __fetch_urls(self, contig_groups, url_end, file_name):
        """
        Fetch information from the NCBI database

        :param contig_groups: a list of lists of identifiers that are submitted
        in the url requests
        :param url_end: the end part of the url defining what type of data and
        return type are requested
        :param file_name: the name of the WGS file, it ensures that the records
        get read into unique files
        :return: a list of file names that contain the requested queries.
        """

        new_files = []
        for index, contig_group in enumerate(contig_groups):
            logging.debug("Fetching url information {}/{}".format(index + 1, len(contig_groups)))
            efetch_url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id='
            efetch_url = efetch_url + ",".join(contig_group) + url_end

            # fetch the url
            url_finished = False
            nrtries = 0
            output = ""
            while not url_finished or nrtries < 4:
                try:
                    nrtries += 1
                    time.sleep(1)
                    with urllib.request.urlopen(efetch_url) as response:
                        output = response.read()
                        if len(output) > 5:
                            url_finished = True
                except (IOError, http.client.BadStatusLine, URLError, http.client.HTTPException):
                    logging.debug("Entry fetching from NCBI failed. Waiting for connection...")
                    time.sleep(5)
            # write the file to appdata to save it.
            if output != "":
                file_loc = TEMP + os.sep + file_name + str(index) + ".gb"
                with open(file_loc, "wb") as f:
                    f.write(output)
                new_files.append(file_loc)
            else:
                logging.debug("Failed to retrieve files. NCBI servers are temporarily unavailable")
        return new_files


class NucleotideDataBase(Database):
    """
    Raw nucleotide database implementation
    """

    def _read_files(self, base_path, paths):
        """
        Read all files into NucleotideContigs. These are contigs with minimal
        information on sequence and a accession number

        :See super._read_files()
        :return a list of NucleotideContig objects
        """
        contigs = []
        while len(paths) > 0:
            path = paths.pop()
            # allow relative paths to b
            path = os.path.join(base_path, path)
            root, ext = os.path.splitext(path)
            if ext in GENBANK_EXTENSIONS:
                new_contigs, new_paths = self.__read_genbank_file(path)
                paths.extend(new_paths)
            elif ext in EMBL_EXTENSIONS:
                new_contigs = self.__read_embl_file(path)
            elif ext in FASTA_EXTENSIONS:
                new_contigs = self.__read_fasta_file(path)
            else:
                logging.warning("Invalid extension {}. Skipping file {}.".format(ext, path))
                continue
            contigs.extend(new_contigs)
        return contigs

    def get_fasta(self):
        """
        :See super.get_fasta()
        """
        fasta_text = ""
        for contig in self._files:
            for index, frame in enumerate(contig.reading_frames):
                fasta_text += ">{}_rf{}\n".format(contig.accession, index + 1)
                fasta_text += "{}\n".format(frame)
        return fasta_text

    def __read_fasta_file(self, file):
        """
        Read the sequence(s) from a fasta file and make sure that the file
        contains DNA code.

        :param file: a string that is the path to a fasta file
        :return: a list of NucleotideContig objects
        """
        fasta_dict = fasta_to_dict(file)
        contigs = []
        for key in fasta_dict:
            if not is_dna(fasta_dict[key]):
                logging.warning("Ignoring fasta entry {} from {}, because it does not contain a DNA sequence."
                                .format(key, file))
                continue
            contigs.append(NucleotideContig(sequence=fasta_dict[key], accession=key, definition="From fasta"))
        return contigs

    def __read_embl_file(self, file):
        """
        Read the sequence(s) from a embl file.

        :param file: a string that is the path to a fasta file
        :return: a list of NucleotideContig objects
        """
        file_text = embl_to_genbank(file)
        contigs_text = file_text.split("//\n")[:-1]
        contigs = []
        for contig_text in contigs_text:
            contigs.append(NucleotideContig(contig_text))
        return contigs

    def __read_genbank_file(self, file):
        """
        Read the sequence(s) from a genbank file. Take into account certain
        parts of the file that can contain references to other records

        :param file: a string that is the path to a fasta file
        :return: a list of NucleotideContig objects and an optional list of
        additional files to be read.
        """
        try:
            with open(file, "r") as f:
                file_text = f.read()
        except Exception:
            logging.critical("Invalid genbank file {}. Exiting...".format(file))
            raise MultiGeneBlastException("Invalid genbank file {}.".format(file))

        if "WGS_SCAFLD  " in file_text or "WGS         " in file_text:
            logging.debug("Encountered WGS master file. Extracting all contigs...")
            root, file_name = os.path.split(file)
            new_file_paths = self._convert_wgs_master_record(file_text, file_name.rsplit(".", 1)[0])
            logging.debug("{} new genbank file(s) where added containing the contigs.".format(len(new_file_paths)))
            return [], new_file_paths
        # TODO handle supercontig records. Example of how to handle can be found in the old parse_gbk.py
        #  file in dblib folder
        elif "CONTIG      " in file_text:
            logging.warning("Supercontig files are not supported at the moment. Skipping...")
            return [], []
        else:
            contigs_text = file_text.split("//\n")[:-1]
            contigs = []
            for contig_text in contigs_text:
                contigs.append(NucleotideContig(contig_text))
            return contigs, []

    def _write_files(self, outdir, dbname):
        """
        Write the contigs of the genbank objects to pickle files.

        :param outdir: an output directory for the database
        :param dbname: the name of all files in the output directory
        """
        index_list = []

        # create a directory for the pickle files, if the directory is there clean it
        try:
            os.mkdir(TEMP + os.sep + "pickles")
        except FileExistsError:
            shutil.rmtree(TEMP + os.sep + "pickles")
            os.mkdir(TEMP + os.sep + "pickles")

        # pickle each contig
        total_contigs = 0
        covered_accessions = set()
        for index, contig in enumerate(self._files):
            if contig.accession in covered_accessions:
                logging.warning("Double accession {} encountered. This is not "
                                "allowed. Skipping writing the second file...".format(contig.accession))
                continue
            with open("{}{}pickles{}{}.pickle".format(TEMP, os.sep, os.sep, index), "wb") as f:
                pickle.dump(contig, f)
            total_contigs += 1
            index_list.append(contig.accession)
        logging.debug("Created {} contig pickles in total.".format(total_contigs))

        # write pickle index file
        with open("{}{}{}_database_index.pickle".format(outdir, os.sep, dbname), "wb") as f:
            pickle.dump(index_list, f)


class ProteinDataBase(Database):
    """
    Protein database implementation
    """

    def get_fasta(self):
        """
        :See super.get_fasta()
        """
        full_fasta = ""
        for gb_file in self._files:
            full_fasta += gb_file.fasta_text()
        return full_fasta

    def _read_files(self, base_path, paths):
        """
        Read a list of files, either genbank or embl and make them into Genbnak
        objects

        :See super._read_files()
        :return: a list of GenbankFile objects
        """
        files = []
        while len(paths) > 0:
            path = paths.pop()
            # allow relative paths to b
            path = os.path.join(base_path, path)
            root, ext = os.path.splitext(path)
            file = None
            if ext in GENBANK_EXTENSIONS:
                file_text, new_paths = self.__read_genbank_file(path)
                paths.extend(new_paths)
            elif ext in EMBL_EXTENSIONS:
                file_text = embl_to_genbank(path)
            else:
                logging.warning("Invalid extension {}. Skipping file {}.".format(ext, path))
                continue
            try:
                if file_text != "":
                    file = GenbankFile(file_text=file_text)
            except MultiGeneBlastException as error:
                logging.debug(error)
                logging.warning("Failed to process the file {}. Skipping...".format(path))
            else:
                if file:
                    files.append(file)
        return files

    def __read_genbank_file(self, file):
        """
        Read a genbank file and check preemptively for invalid files

        :param file: a file path
        :return: a large string that is the text in the geanbank file
        """
        try:
            with open(file, "r") as f:
                file_text = f.read()
        except Exception:
            logging.critical("Invalid genbank file {}. Exiting...".format(file))
            raise MultiGeneBlastException("Invalid genbank file {}.".format(file))

        new_file_paths = []
        # make sure to remove potential old occurances of \r. Acts like a \n
        file_text = file_text.replace("\r", "\n")
        # if the file is a WGS master file disect it
        if "WGS_SCAFLD  " in file_text or "WGS         " in file_text:
            logging.debug("Encountered WGS master file. Extracting all contigs...")
            root, file_name = os.path.split(file)
            new_file_paths = self._convert_wgs_master_record(file_text, file_name.rsplit(".", 1)[0])
            logging.debug("{} new genbank file(s) where added containing the contigs.".format(len(new_file_paths)))
            file_text = ""
        elif "CONTIG      " in file_text:
            file_text = ""
            logging.warning("Supercontig files are not supported at the moment. Skipping...")
        # do a basic check to see if the genbank file is valid
        elif "     CDS             " not in file_text or "\nORIGIN" not in file_text:
            logging.warning("Genbank file {} is not properly formatted or contains no sequences. Skipping..."
                            .format(file))
        return file_text, new_file_paths

    def _write_files(self, outdir, dbname):
        """
        :See super._write_files
        """
        index_list = []

        # create a directory for the pickle files, if the directory is there clean it
        try:
            os.mkdir(TEMP + os.sep + "pickles")
        except FileExistsError:
            shutil.rmtree(TEMP + os.sep + "pickles")
            os.mkdir(TEMP + os.sep + "pickles")

        # pickle each contig
        contig_count = 0
        covered_accessions = set()
        for gb_file in self._files:
            for contig_name in gb_file.contigs:
                if contig_name in covered_accessions:
                    logging.warning("Double accession {} encountered. This is not "
                                    "allowed. Skipping writing the second file...".format(contig_name))
                    continue
                with open("{}{}pickles{}{}.pickle".format(TEMP, os.sep, os.sep, contig_count), "wb") as f:
                    pickle.dump(gb_file.contigs[contig_name], f)
                contig_count += 1
                covered_accessions.add(contig_name)
                index_list.append(set(gb_file.contigs[contig_name].proteins.keys()))
        logging.debug("Created {} contig pickles in total.".format(contig_count))

        # write pickle index file
        with open("{}{}{}_database_index.pickle".format(outdir, os.sep, dbname), "wb") as f:
            pickle.dump(index_list, f)


class GenbankFile:
    """
    An Object for easily disecting a genbank file into the individual proteins
    """
    def __init__(self, path=None, file_text=None, contigs=None, protein_range=None, allowed_proteins=None):
        """
        :param path: a path to a genbank file
        :param file_text: an optional parameter that will be used as the file text
        instead of the read version of the file.
        :param contigs: a optional list of contigs to skip reading the file
        :param protein_range: an optional range in which proteins can be selected
        from the file instead of all proteins
        :param allowed_proteins: a list of protein IDs that can be selected from
        the genbank file instead of all proteins
        """
        if protein_range is None:
            protein_range = [None for _ in range(2)]
        if contigs is not None:
            self.contigs = {cont.accession: cont for cont in contigs}
        elif file_text is not None:
            self.contigs = self.__create_contigs(file_text, protein_range, allowed_proteins)
        elif path is not None:
            file_text = self.__read_genbank_file(path)
            self.contigs = self.__create_contigs(file_text, protein_range, allowed_proteins)
        else:
            raise MultiGeneBlastException("One of the three options; path, file_text,"
                                          " or contigs need to be selected")
        # an OrderedDict.
        self.proteins = self.__list_proteins()
        self.lenght = sum(con.lenght for con in self.contigs.values())

    def fasta_text(self):
        """
        Give a fatsa text that is the combined fasta of all proteins in the
        database

        :return: a String containing the combined fasta of all unique proteins
        in the genbank file.
        """
        text = ""
        for protein in self.proteins.values():
            text += protein.fasta_text()
        return text

    def __list_proteins(self):
        """
        List all unique proteins from all contigs in a OrderedDict

        :return: an OrderedDict
        """
        protein_dict = OrderedDict()
        for contig in self.contigs.values():
            prots = contig.proteins.values()
            for prot in prots:
                if prot.name in protein_dict:
                    logging.warning("Double fasta entry '{}'. Skipping...".format(prot.name))
                else:
                    protein_dict[prot.name] = prot
        return protein_dict

    def __read_genbank_file(self, file):
        """
        Read a genbank file and check preemptively for invalid files

        :param file: a file path
        :return: a large string that is the text in the geanbank file
        """
        try:
            with open(file, "r") as f:
                file_text = f.read()
        except Exception:
            logging.critical("Invalid genbank file {}. Exiting...".format(file))
            raise MultiGeneBlastException("Invalid genbank file {}.".format(file))

        # make sure to remove potential old occurances of \r. Acts like a \n
        file_text = file_text.replace("\r", "\n")

        # if there is a WGS in the genbank record make sure to raise an error.
        if "WGS_SCAFLD  " in file_text or "WGS         " in file_text:
            logging.critical("WGS master files cannot be used as query input files.")
            raise MultiGeneBlastException("WGS master files cannot be used as query input files.")
        # do a basic check to see if the genbank file is valid
        if "     CDS             " not in file_text or "\nORIGIN" not in file_text:
            logging.critical("Genbank file {} is not properly formatted or contains no sequences".format(file))
            raise MultiGeneBlastException("Genbank file {} is not properly formatted or contains no sequences"
                                          .format(file))
        return file_text

    def __create_contigs(self, text, protein_range, allowed_proteins):
        """
        Create an unordered dictionary of ProteinContig objects

        :param text: Text that contains one contig
        :param protein_range: an optional range in which proteins can be selected
        from the file instead of all proteins
        :param allowed_proteins: a list of protein IDs that can be selected from
        the genbank file instead of all proteins
        :return: a dictionary of ProteinContig objects
        """
        contigs = {}
        for cont in text.split("//\n")[:-1]:
            contig = ProteinContig(cont, protein_range, allowed_proteins)
            contigs[contig.accession] = contig
        return contigs


class Contig(ABC):
    """
    Base abstraction class for contigs. Reads the header.
    """
    COUNTER = 0

    def __init__(self, dna_sequence, gene_information):
        """
        :param dna_sequence: The full DNA sequence in string form present in the genbank file for a certain contig
        :param gene_information: a list of large pieces of text containing the information of individual genes
        """
        self.lenght = len(dna_sequence)

        # read general information
        self.accession = ""
        self.definition = ""
        self.__extract_header(gene_information[0])

    def __extract_header(self, header):
        """
        Disect the header of the contig and extract the accesion and
        definition

        :param header: a string that contains all the information in the header
        of the genbank file
        """
        header_lines = header.split("\n")
        for index, line in enumerate(header_lines):
            if line.lower().startswith("accession"):
                # extract the accesion
                self.accession = remove_illegal_characters(line.replace("ACCESSION   ", ""))
            if line.lower().startswith("definition"):
                self.definition += line.replace("DEFINITION  ", "").strip()
                for def_line in header_lines[index + 1:]:
                    # end of definition line
                    if not def_line.startswith("    "):
                        break
                    self.definition += def_line.strip()
        # Test if accession number is probably real GenBank/RefSeq acc nr
        if not is_valid_accession(self.accession):
            logging.debug("Probably invalid GenBank/Refseq accesion {} found.".format(self.accession))
            self.accession = ""
        # auto generate a unique accession when no valid one is found
        if self.accession == "":
            self.accession = "accession{}".format(self.COUNTER)
            Contig.COUNTER += 1
            logging.debug("The following accession was auto generated: {}. Because no valid accession was found."
                          .format(self.accession))
        if self.definition == "":
            logging.debug("No definition found for contig {}".format(self.accession))


class NucleotideContig(Contig):
    """
    Abstraction for nucleotide contigs, saves a sequence and header information
    """
    def __init__(self, file_text=None, sequence=None, accession="", definition=""):
        """
        :param file_text: string to be read.
        :param sequence: Instead of providing a file_text
        :param accession: Optional name of an accession. If no name is provided the accession is auto generated
        :param definition: Optional definition to give some extra information
        """
        if file_text is not None:
            gene_defenitions, dna_sequence = file_text.split("\nORIGIN")
            # Extract DNA sequence and calculate complement of it
            dna_sequence = clean_dna_sequence(dna_sequence)

            # Extract gene information
            gene_information = gene_defenitions.split("     CDS             ")

            super().__init__(dna_sequence, gene_information)
        elif sequence is not None:
            dna_sequence = sequence
            self.lenght = len(dna_sequence)
            if accession == "":
                self.accession = "accession{}".format(self.COUNTER)
                Contig.COUNTER += 1
            else:
                self.accession = accession
            self.definition = definition
        else:
            raise MultiGeneBlastException("When creating a NucleotideContig you either have to supply a "
                                          "file text or sequence.")
        self.reading_frames = self.__get_reading_frames(dna_sequence)

    def __get_reading_frames(self, dna_seq):
        rf1 = translate(dna_seq)
        rf2 = translate(dna_seq[1:])
        rf3 = translate(dna_seq[2:])
        reverse_seq = dna_seq[::-1]
        rf4 = translate(reverse_seq)
        rf5 = translate(reverse_seq[1:])
        rf6 = translate(reverse_seq[2:])
        return rf1, rf2, rf3, rf4, rf5, rf6


class ProteinContig(Contig):
    """
    Acts as abstraction for a contig present in a genbank file
    """
    def __init__(self, file_text, protein_range, allowed_proteins):
        """
        :param file_text: string to be read.

        These last 2 options are mainly here for allowing filtering of the
        query contig.

        :param protein_range: a range of nucleotides wherein proteins are
        considered part of the contig to read.
        :param allowed_proteins: a list of protein names that are allowed
         to be in the contig
        """
        gene_defenitions, dna_sequence = file_text.split("\nORIGIN")
        # Extract DNA sequence and calculate complement of it
        dna_sequence = clean_dna_sequence(dna_sequence)

        # Extract gene information
        gene_information = gene_defenitions.split("     CDS             ")

        super().__init__(dna_sequence, gene_information)

        c_dna_sequence = complement(dna_sequence)

        # exract all the entries from the genbank file
        self.entries = self.__extract_entries(gene_information[1:], dna_sequence, c_dna_sequence, protein_range,
                                              allowed_proteins)
        self.entries.sort(key=lambda x: (x.protein.start, x.protein.stop))
        self.proteins = self.__get_unique_proteins()

    def __get_unique_proteins(self):
        """
        Make sure that you retrieve a unique list of named proteins after names
        have been cleaned for illegal characters. Also filter out entries that
        mis sequences.

        :return: a dictionary of protein objects
        """
        protein_dict = OrderedDict()
        for entry in self.entries:
            if entry.protein.name in protein_dict:
                logging.warning("Double protein entry '{}' for contig '{}'. Skipping...".format(entry.protein.name,
                                                                                                self.accession))
            elif len(entry.protein.sequence) == 0:
                logging.warning("Cannot find or reconstruct sequence for entry {}. Skipping..."
                                .format(entry.protein.name))
            else:
                protein_dict[entry.protein.name] = entry.protein
        return protein_dict

    def __extract_entries(self, entries, dna_seq, c_dna_seq, protein_range, allowed_proteins):
        """
        Extract genbank protein entries from a list of strings that contain
        these entries

        :param entries: a list of strings that contain one CDS entry each
         :param dna_seq: the dna sequence present in the genbank file
        :param c_dna_seq: the complement sequence present in the genbank file
        :param protein_range: an optional range in which proteins can be selected
        from the file instead of all proteins
        :param allowed_proteins: a list of protein IDs that can be selected from
        the genbank file instead of all proteins
        :return: a list of GenbakEnrty objects
        """
        genbank_entries = []

        # number to check if all allowed proteins have been found or not
        # this allows the function to stop prematurely if all proteins are found
        found_proteins = 0
        # ignore the firs hit which is the start of the genbnk file
        for index, gene in enumerate(entries):
            g = GenbankEntry(gene, index + 1, dna_seq, c_dna_seq, self.accession, self.definition)
            # skip if no valid protein is found
            if g.protein is None:
                continue
            elif protein_range[0] is not None and g.protein.start >= protein_range[0] and \
                    g.protein.stop <= protein_range[1]:
                genbank_entries.append(g)
            elif allowed_proteins is not None and (g.protein.protein_id in allowed_proteins or g.protein.name in
                                                   allowed_proteins or g.protein.locus_tag in allowed_proteins):
                genbank_entries.append(g)
                found_proteins += 1
                if found_proteins >= len(allowed_proteins):
                    break
            elif allowed_proteins is None and protein_range[0] is None:
                genbank_entries.append(g)
        return genbank_entries


class GenbankEntry:
    """
    Disect a genbank entry and substract allot of potential values.
    """
    def __init__(self, entry_string, nr, dna_seq, c_dna_seq, contig_id, contig_description):
        """
        :param entry_string: a string that contains a CDS genbank entry
        :param nr: a unique integer that can be used to create a genename is no
        name is available
        :param dna_seq: the dna sequence present in the genbank file
        :param c_dna_seq: the complement sequence present in the genbank file
        :param contig_id: a generated id for the contig that the entry belongs to
        :param contig_description: a description of the contig the entry belongs to
        """
        self.location = None
        self.contig_id = contig_id
        self.contig_description = contig_description

        self.__codon_start = 1
        self.locus_tag = ""
        self.gene_name = ""
        self.protein_id = ""
        self.protein_code = ""
        self.annotation = ""

        self.__disect_string(entry_string, nr)
        self.protein = self.__create_protein(dna_seq, c_dna_seq)

    def __disect_string(self, entry_string, nr):
        """
        Parse the entry_string and try to extract allot of different potentialy
        available values

        :param entry_string: a string that contains a CDS genbank entry
        :param nr: a unique integer that can be used to create a genename is no
        name is available
        """
        # create a list of lines that are part of the coding sequences (CDS)
        cds_part = self.__cds_lines(entry_string)
        for index, line in enumerate(cds_part[1:]):
            line = line.replace("/", "").strip()
            if line.lower().startswith("codon_start="):
                self.__codon_start = int(self.__clean_line(line))
            elif line.lower().startswith("protein_id="):
                self.protein_id = self.__clean_line(line)
            elif line.lower().startswith("locus_tag="):
                self.locus_tag = self.__clean_line(line)
            elif line.lower().startswith("gene="):
                self.gene_name = self.__clean_line(line)
            elif line.lower().startswith("translation="):
                self.protein_code += self.__clean_line(line)
                for prot_line in cds_part[index + 2:]:
                    self.protein_code += prot_line.strip().replace('"', "")
                break
            elif line.lower().startswith("product="):
                self.annotation = self.__clean_line(line)

        # configure the genename with this order, locus tag, protein_id, gene_name auto generated name
        if self.locus_tag:
            self.gene_name = self.locus_tag
        elif self.protein_id:
            self.gene_name = self.protein_id
        elif self.gene_name:
            pass
        else:
            self.gene_name = "orf {}".format(nr)

        # extract the locations
        self.locations = self.__get_locations(cds_part[0].strip())

    def __cds_lines(self, string):
        """
        Extract all the CDS location qualifiers form the string. This separates
        them from other potential qualifiers

        :param string: a string that contains a CDS genbank entry
        :return: a list of lines that are part of the CDS location qualifier
        """
        potential_lines = string.split("\n")
        # include the first line always
        lines = [potential_lines[0]]
        for line in potential_lines[1:]:
            if not line.startswith("                     "):
                break
            lines.append(line)
        return lines

    def __clean_line(self, line):
        """
        Clean a line from the genbank file by removing unwanted characters
        as well as the name associated with the value

        :param line: a string containing no newlines
        :return: the cleaned version of the input line
        """
        # remove unwanted characters and only inlcude the value not the name
        try:
            return line.replace("\\", "_").replace("/", "_").replace('"', '').replace(" ", "_").split("=")[1]
        except IndexError:
            logging.warning("{} does not have the right format, as a consequence will be ignored".format(line))

    def __get_locations(self, string):
        """
        Find the location in a location string. This can be tricky because of
        complement and join modifiers.
        Note: this script does not take the order modifier into account, and i
        cannot see why that would be neccesairy

        :param string: a string that should contain a location of minimal format
        start..stop
        :return: a list of lists that holds a list of locations (to take join
        into account) and a boolean that tells if the locations are on the
        complement string or not. A location has the format [start, stop] but
        can be None when no valid location was found
        """
        if string.startswith("complement"):
            if "join" in string:
                return self.__disect_join(string)
            else:
                # split the string without the complement text
                return [self.__convert_location(string[11:-1])], True

        elif string.startswith("join"):
            return self.__disect_join(string)
        else:
            return [self.__convert_location(string)], False

    def __disect_join(self, string):
        """
        Disect a location that contains a join statement

        :param string: a location in string format that has a join format
        :return:
        """
        # remove the join tag
        comp = False
        if string.startswith("complement"):
            comp = True
            string = string[11:-1]
        string = string[5:-1]
        individual_strings = string.split(",")
        locations = []
        for s in individual_strings:
            locations.append(self.__convert_location(s))
        return locations, comp

    def __convert_location(self, string):
        """
        Convert the location string s of the format start..stop into two integers
        while removign potential present characters

        :param string: a string containing a location
        :return: a list of lenght 2 with 2 integers
        """
        if ".." not in string:
            logging.warning("No valid location specifier found for genbank entry {}".format(self.gene_name))
            return None
        try:
            str_locations = string.replace("<", "").replace(">", "").strip().split("..")
            return list(map(int, str_locations))
        # in case the str locations cannot be converted to integers
        except AttributeError:
            logging.warning("No valid location specifier found for genbank entry {}".format(self.gene_name))
            return None

    def __create_protein(self, dna_seq, c_dna_seq):
        """
        Create a Protein object using class variables extracted from the entry
        string

        :param dna_seq: the dna sequence present in the genbank file
        :param c_dna_seq: the complement sequence present in the genbank file
        :return: a Protein object
        """
        # loc is a list of start, stop and True or false for reverse complement or not
        locations, reverse_complement = self.locations
        strand = "+"
        if reverse_complement:
            strand = "-"
        # make sure locations are sorted
        if None in locations:
            logging.warning("Can not process {}. No valid start stop found. Skipping...".format(self.gene_name))
            return
        locations.sort()

        # if a sequence is available use that, otherwise extract it
        if self.protein_code:
            sequence = self.protein_code
        else:
            sequence = ""
            for loc in locations:
                start, end = loc
                if reverse_complement:
                    nc_seq = c_dna_seq[max(0, (start - 1) - self.__codon_start - 1): end]
                else:
                    nc_seq = dna_seq[max(0, (start - 1) - self.__codon_start - 1): end]
                sequence += nc_seq
            sequence = translate(sequence)

        # the first location and then the start
        start = locations[0][0]
        return Protein(sequence, start, self.gene_name, strand,
                       locus_tag=self.locus_tag, annotation=self.annotation,
                       contig_id=self.contig_id, protein_id=self.protein_id,
                       contig_description=self.contig_description)


class Protein:
    """
    A representation of a protein, containing some optional fields that can or
    cannot be specified.
    """
    def __init__(self, sequence, start, name, strand, end=None, annotation="",
                 locus_tag="", contig_id="", contig_description="", protein_id=""):
        """
        :param sequence: a string that is the amino acid sequence of the protein
        :param start: an integer that is the start coordinate in the the protein
        in the larger contig it is located
        :param name: the gene name of the protein
        :param strand: the + or - strand the protein can be located on
        :param end: optional parameter to set the stop of the protein manually
        :param annotation: the optional annotation of the protein, description of
        function
        :param locus_tag: an optional locus tag of the protein
        :param contig_id: The ID of the contig as defined in the genbank file
        :param contig_description: the description of the contig as defined in
        the genbank file
        :param protein_id: an optional id of the protein
        """
        self.sequence = sequence
        self.strand = strand
        self.aa_lenght = len(self.sequence)
        self.nt_lenght = self.aa_lenght * 3

        # both in nucleotides
        self.start = start
        if end:
            self.stop = end
        else:
            self.stop = self.start + self.nt_lenght

        # gene name, this name should be used when comparing proteins.
        self.name = remove_illegal_characters(name)
        self.annotation = remove_illegal_characters(annotation)

        self.locus_tag = locus_tag
        self.protein_id = protein_id
        self.contig_id = contig_id
        self.contig_description = contig_description

    def summary(self):
        """
        Give a short string that is a summary of the protein

        :return: a string
        """
        lt = self.locus_tag
        if self.locus_tag == "":
            lt = "No locus tag"
        return "{}\t{}\t{}\t{}\t{}\t{}".format(self.name, self.start, self.stop,
                                               self.strand, self.annotation, lt)

    def fasta_text(self):
        """
        Full fasta entry for writing into a fasta file. Uses the pre computed
        header and the sequence

        :return: a string that is a complete fasta entry
        """
        text = ">{}\n{}\n".format(self.name, self.sequence)
        return text
