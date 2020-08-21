
# imports
import logging
import tarfile
import pickle
import os
import shutil
from string import ascii_letters
from collections import OrderedDict

import urllib.request, urllib.error, urllib.parse
from urllib.request import Request, urlopen
from urllib.error import URLError, HTTPError
import http.client
from http.client import BadStatusLine,HTTPException

from utilities import *
from constants import GENBANK_EXTENSIONS, EMBL_EXTENSIONS, TEMP


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
    except Exception as e:
        logging.critical("Invalid embl file {}. Exiting...".format(embl_filepath))
        raise MultiGeneBlastException("Invalid embl file {}.".format(embl_filepath))

    # make sure to remove potential old occurances of \r. Acts like a \n
    file_text = file_text.replace("\r", "\n")

    # do a basic check to see if the embl file is valid
    if "FT   CDS " not in file_text or ("\nSQ" not in file_text):
        logging.critical("Embl file {} is not properly formatted or contains no sequences. Exiting...".format(embl_filepath))
        raise MultiGeneBlastException("Embl file {} is not properly formatted or contains no sequences".format(embl_filepath))
    text_lines = file_text.split("\n")

    #change certain line starts
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
            #change all the definition lines to make them readable
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


class DataBase:
    """
    Tracks multiple genbank and embl files and can write the output into a format
    that can be used by MultiGeneBlast
    """
    def __init__(self, base_path, paths):
        logging.info("Started creating database...")
        #care this is a generator object not a list
        self.__gb_files = self.__read_files(base_path, paths)
        if len(self.__gb_files) == 0:
            logging.critical("Failed to load any of the provided database files.")
            raise MultiGeneBlastException("Failed to load any of the provided database files.")

    def create(self, outdir, dbname):
        """
        Create the files for the database using the self.__gb_files genbank
        objects

        :param outdir: an output directory for the database
        :param dbname: the name of all files in the output directory
        """
        self.__write_files(outdir, dbname)
        self.__create_tar_file(outdir, dbname)

    def get_fasta(self):
        """
        Convenience method for getting the fasta strings from all Genbank Objects
        :return: a String
        """
        full_fasta = ""
        for gb_file in self.__gb_files:
            full_fasta += gb_file.fasta_text()
        return full_fasta

    def __read_files(self, base_path, paths):
        """
        Read a list of files, either genbank or embl and make them into Genbnak
        objects

        :param base_path: an absolute path to the directory make_database.py was
        executed from
        :param paths: a absolute or relative path to all the database files from
        the make_database.py directory
        :return: a list of Genbnak objects
        """
        files = []
        while len(paths) > 0:
            path = paths.pop()
            #allow relative paths to b
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
                if file_text:
                    file = GenbankFile(file_text=file_text)
            except MultiGeneBlastException as error:
                print(error)
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
        except Exception as e:
            logging.critical("Invalid genbank file {}. Exiting...".format(file))
            raise MultiGeneBlastException("Invalid genbank file {}.".format(file))

        new_file_paths = []
        # make sure to remove potential old occurances of \r. Acts like a \n
        file_text = file_text.replace("\r", "\n")
        #if the file is a WGS master file disect it
        #TODO consider checking if the user has an internet connection
        if "WGS_SCAFLD  " in file_text or "WGS         " in file_text:
            logging.debug("Encountered WGS master file. Extracting all contigs...")
            root, file_name = os.path.split(file)
            new_file_paths = self.__convert_wgs_master_record(file_text, file_name.rsplit(".", 1)[0])
            logging.debug("{} new genbank file(s) where added containing the contigs.".format(len(new_file_paths)))
            file_text = ""
        # do a basic check to see if the genbank file is valid
        elif "     CDS             " not in file_text or "\nORIGIN" not in file_text:
            logging.warning("Genbank file {} is not properly formatted or contains no sequences. Skipping...".format(file))
        return file_text, new_file_paths

    def __write_files(self, outdir, dbname):
        """
        Write the contigs of the genbank objects to pickle files.

        :param outdir: an output directory for the database
        :param dbname: the name of all files in the output directory
        """
        #TODO think about handling double contigs. I would think not to relevant
        index_list = []

        #create a directory for the pickle files, if the directory is there clean it
        try:
            os.mkdir(TEMP + os.sep + "pickles")
        except FileExistsError:
            shutil.rmtree(TEMP + os.sep + "pickles")
            os.mkdir(TEMP + os.sep + "pickles")

        #pickle each contig
        total_pickles = 0
        for gb_file in self.__gb_files:
            for index, contig in enumerate(gb_file.contigs):
                with open("{}{}pickles{}{}.pickle".format(TEMP, os.sep, os.sep, index), "wb") as f:
                    pickle.dump(gb_file.contigs[contig], f)
                    total_pickles += 1
                index_list.append(set(gb_file.contigs[contig].proteins.keys()))
        logging.debug("Created {} contig pickles in total.".format(total_pickles))

        #write pickle index file
        with open("{}{}{}_database_index.pickle".format(outdir, os.sep, dbname), "wb") as f:
            pickle.dump(index_list, f)

    def __create_tar_file(self, outdir, dbname):
        """
        Create gzipped tar file to save all the pickled contigs into.

        :param outdir: an output directory for the database
        :param dbname: the name of all files in the output directory
        """
        logging.info("Writing to tar file...")
        with tarfile.open("{}{}{}_contigs.tar.gz".format(outdir, os.sep, dbname), "w:gz") as tar:
            tar.add("{}{}pickles".format(TEMP, os.sep), arcname=os.path.basename(dbname))

    def __convert_wgs_master_record(self, text, file_name):
        #If seq_record is a WGS master record, parse out contig accession numbers and download these as separate seq_records
        if "WGS_SCAFLD  " in text:
            contig_ranges = text.split("WGS_SCAFLD  ")[1].split("\n")[0]
        elif "WGS         " in text:
            contig_ranges = text.split("WGS         ")[1].split("\n")[0]

        #extract all the ranges
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

        #extract all contigs for each range
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
            #certain entries rely on a number of zeroes. This can be at the start of a number
            for char in startnumber:
                if char == "0":
                    nrzeros += 1
                else:
                    break
            contig_range = [alpha_tag + nrzeros * "0" + str(number) for number in range(int(startnumber), int(endnumber))]
            allcontigs.extend(contig_range)

        #Create contig groups of 50 (reasonable download size per download)
        nr_groups = len(allcontigs) / 50 + 1
        contig_groups = [allcontigs[i:i+50] for i in range(0, len(allcontigs), 50)]

        #Download contigs and parse into seq_record objects
        new_files = self.__fetch_urls(contig_groups, '&rettype=gbwithparts&retmode=text', file_name)
        return new_files

    def __fetch_urls(self, contig_groups, url_end, file_name):
        new_files = []
        for index, contig_group in enumerate(contig_groups):
            logging.debug("Fetching url information {}/{}".format(index + 1, len(contig_groups)))
            efetch_url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id='
            efetch_url = efetch_url + ",".join(contig_group) + url_end

            #fetch the url
            url_finished = False
            nrtries = 0
            output = ""
            while not url_finished or nrtries < 4:
                try:
                    nrtries += 1
                    time.sleep(3)
                    req = urllib.request.Request(efetch_url)
                    response = urllib.request.urlopen(req)
                    output = response.read()
                    if len(output) > 5:
                        url_finished = True
                except (IOError, http.client.BadStatusLine, URLError,http.client.HTTPException):
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


class GenbankFile:
    """
    An Object for easily disecting a genbank file into the individual proteins
    """
    def __init__(self, path = None, file_text = None, contigs = None, protein_range=None, allowed_proteins=None):
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
        if contigs != None:
            self.contigs = {cont.accession:cont for cont in contigs}
        elif file_text != None:
            self.contigs = self.__create_contigs(file_text, protein_range, allowed_proteins)
        elif path != None:
            file_text = self.__read_genbank_file(path)
            self.contigs = self.__create_contigs(file_text, protein_range, allowed_proteins)
        else:
            raise MultiGeneBlastException("One of the three options; path, file_text,"
                                          " or contigs need to be selected")
        self.proteins = self.__list_proteins() #an OrderedDict.
        self.lenght = sum(con.lenght for con in self.contigs.values())


    def fasta_text(self):
        """
        Give a fatsa text that is the combined fasta of all proteins in the
        database
        :return:
        """
        text = ""
        for protein in self.proteins.values():
            text += protein.fasta_text()
        return text

    def __list_proteins(self):
        protein_dict = OrderedDict()
        for contig in self.contigs.values():
            prots = contig.proteins.values()
            for prot in prots:
                if prot.name in protein_dict:
                    logging.warning("Double fasta entry '{}' in file '{}'. Skipping...".format(entry.protein.name, self.file))
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
        except Exception as e:
            logging.critical("Invalid genbank file {}. Exiting...".format(file))
            raise MultiGeneBlastException("Invalid genbank file {}.".format(file))

        # make sure to remove potential old occurances of \r. Acts like a \n
        file_text = file_text.replace("\r", "\n")

        # do a basic check to see if the genbank file is valid
        if "WGS_SCAFLD  " in file_text or "WGS         " in file_text:
            root, file_name = os.path.split(file)
            file_text = convert_wgs_master_record(file_text, file_name)
        if "     CDS             " not in file_text or "\nORIGIN" not in file_text:
            logging.critical("Genbank file {} is not properly formatted or contains no sequences".format(file))
            raise MultiGeneBlastException("Genbank file {} is not properly formatted or contains no sequences".format(file))
        return file_text

    def __create_contigs(self, text, protein_range, allowed_proteins):
        contigs = {}
        for cont in text.split("//\n")[:-1]:
            contig = Contig(cont, protein_range, allowed_proteins)
            contigs[contig.accession] = contig
        return contigs


class Contig:
    """
    Acts as abstraction for a contig present in a genbank file
    """
    def __init__(self, file_text, protein_range=None, allowed_proteins=None):
        gene_defenitions, dna_sequence = file_text.split("\nORIGIN")
        # Extract DNA sequence and calculate complement of it
        dna_sequence = clean_dna_sequence(dna_sequence)
        c_dna_sequence = complement(dna_sequence)

        self.lenght = len(dna_sequence)

        # Extract gene information
        gene_information = gene_defenitions.split("     CDS             ")

        #read general information
        self.accession = ""
        self.definition = ""
        self.__extract_header(gene_information[0])

        #exract all the entries from the genbank file
        self.entries = self.__extract_entries(gene_information[1:], dna_sequence, c_dna_sequence, protein_range, allowed_proteins)
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
                logging.warning("Double fasta entry '{}' in file '{}'. Skipping...".format(entry.protein.name, self.file))
            elif len(entry.protein.sequence) == 0:
                logging.warning("Cannot find or reconstruct sequence for entry {}. Skipping...".format(entry.protein.name))
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

        #number to check if all allowed proteins have been found or not
        #this allows the function to stop prematurely if all proteins are found
        found_proteins = 0
        # ignore the firs hit which is the start of the genbnk file
        for index, gene in enumerate(entries):
            g = GenbankEntry(gene, index + 1, dna_seq, c_dna_seq, self.accession, self.definition)
            #skip if no valid protein is found
            if g.protein == None:
                continue
            elif protein_range != None and g.protein.start >= protein_range[0] and g.protein.stop <= protein_range[1]:
                genbank_entries.append(g)
            elif allowed_proteins != None and g.protein.protein_id in allowed_proteins:
                genbank_entries.append(g)
                found_proteins += 1
                if found_proteins >= len(allowed_proteins):
                    break
            elif allowed_proteins == None and protein_range == None:
                genbank_entries.append(g)
        return genbank_entries

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
                self.accession = line.replace("ACCESSION   ","")
            if line.lower().startswith("definition"):
                self.definition += line.replace("DEFINITION  ", "").strip()
                for def_line in header_lines[index + 1:]:
                    #end of definition line
                    if not def_line.startswith("    "):
                        break
                    self.definition += def_line.strip()
        # Test if accession number is probably real GenBank/RefSeq acc nr
        if not is_valid_accession(self.accession):
            logging.debug("Probably invalid GenBank/Refseq accesion {} found.".format(self.accession))
            self.accession = ""
        if self.accession == "":
            logging.debug("No valid accesion found for file {}".format(self.file))
        if self.definition == "":
            logging.debug("No definition found for file {}".format(self.file))


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
        :param file_accession: the genbank file accession this entry originated
        from
        """
        self.location = None
        self.contig_id = contig_id
        self.contig_description = contig_description

        #relevant attributes that can be present
        self.__codon_start = 1

        #optional parameters that can be found in the entry
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
        #create a list of lines that are part of the coding sequences (CDS)
        cds_part = self.__CDS_lines(entry_string)
        protein_id = locus_tag = gene_name = None
        for index, line in enumerate(cds_part[1:]):
            line = line.replace("/","").strip()
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

        #configure the genename with this order, locus tag, gene_name protein_id auto generated name
        if self.locus_tag:
            self.gene_name = self.locus_tag
        elif self.gene_name:
            pass
        elif self.protein_id:
            self.gene_name = self.protein_id
        else:
            self.gene_name = "orf {}".format(nr)

        #extract the locations
        self.locations = self.__get_locations(cds_part[0].strip())

    def __CDS_lines(self, string):
        """
        Extract all the CDS location qualifiers form the string. This separates
        them from other potential qualifiers

        :param string: a string that contains a CDS genbank entry
        :return: a list of lines that are part of the CDS location qualifier
        """
        potential_lines = string.split("\n")
        #include the first line always
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
        #remove unwanted characters and only inlcude the value not the name
        try:
            return line.replace("\\", "_").replace("/", "_").replace('"', '').replace(" ","_").split("=")[1]
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
                #split the string without the complement text
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
        #remove the join tag
        complement = False
        if string.startswith("complement"):
            complement = True
            string = string[11:-1]
        string = string[5:-1]
        individual_strings = string.split(",")
        locations = []
        for s in individual_strings:
            locations.append(self.__convert_location(s))
        return locations, complement

    def __convert_location(self, s):
        """
        Convert the location string s of the format start..stop into two integers
        while removign potential present characters

        :param s: a string containing a location
        :return: a list of lenght 2 with 2 integers
        """
        if ".." not in s:
            logging.warning("No valid location found for genbank entry {}".format(self.gene_name))
            return None
        try:
            str_locations = s.replace("<", "").replace(">", "").strip().split("..")
            return list(map(int, str_locations))
        #in case the str locations cannot be converted to integers
        except AttributeError:
            logging.warning("No valid location found for genbank entry {}".format(self.gene_name))
            return None

    def __create_protein(self, dna_seq, c_dna_seq):
        """
        Create a Protein object using class variables extracted from the entry
        string

        :param dna_seq: the dna sequence present in the genbank file
        :param c_dna_seq: the complement sequence present in the genbank file
        :return: a Protein object
        """
        #loc is a list of start, stop and True or false for reverse complement or not
        locations, reverse_complement = self.locations
        strand = "+"
        if reverse_complement:
            strand = "-"
        #make sure locations are sorted
        if None in locations:
            logging.warning("Can not process {}. No valid start stop found. Skipping...".format(self.gene_name))
            return
        locations.sort()

        #if a sequence is available use that, otherwise extract it
        if self.protein_code:
            sequence = self.protein_code
        else:
            sequence = ""
            for loc in locations:
                start, end = loc
                if reverse_complement:
                    nc_seq = c_dna_seq[(start - 1) - self.__codon_start - 1 : end]
                else:
                    nc_seq = dna_seq[(start - 1) - self.__codon_start - 1 : end]
                sequence += nc_seq
            sequence = translate(sequence)

        #the first location and then the start
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
    def __init__(self, sequence, start, name, strand, end = None, annotation = "",
                 locus_tag = "", contig_id = "", contig_description = "", protein_id = ""):
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
        :param start_header: the start of the fasta header. This can be used in
        the case of query proteins to allow certain identifyers
        """
        self.sequence = sequence
        self.strand = strand
        self.aa_lenght = len(self.sequence)
        self.nt_lenght = self.aa_lenght * 3

        #both in nucleotides
        self.start = start
        if end:
            self.stop = end
        else:
            self.stop = self.start + self.nt_lenght

        #gene name, this name should be used when comparing proteins.
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
