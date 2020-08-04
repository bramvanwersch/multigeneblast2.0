import logging
from string import ascii_letters
ILLEGAL_CHARACTERS = ["'",'"','=',';',':','[',']','>','<','|','\\',"/",'*','-','_','.',',','?',')','(','^','#','!','`','~','+','{','}','@','$','%','&']


class Protein:
    """
    A representation of a protein, containing some optional fields that can or
    cannot be specified.
    """
    def __init__(self, sequence, start, name, strand, end = None, annotation = "",
                 locus_tag = "", genbank_file = "", protein_id = "", start_header = ""):
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
        :param genbank_file: an optional genbank file accesion where the protein
        originated from
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

        #gene name
        self.name = remove_illegal_characters(name)
        self.annotation = annotation

        self.locus_tag = locus_tag
        self.protein_id = protein_id
        self.genbank_file = genbank_file
        self.fasta_header = self.__construct_fasta_header(start_header)

    def __construct_fasta_header(self, start_header):
        """
        Construct a fasta header using available class variables

        :param start_header: start of the header, can be ''
        :return: a string that is the text of the fasta header. Meaning the > is
        not included
        """
        header = "{}-{}|{}|{}|{}|{}|{}".format(self.start, self.stop, self.strand, self.name, self.annotation,self.protein_id, self.locus_tag)
        return "{}|{}".format(start_header, header)

    def fasta_text(self):
        """
        Full fasta entry for writing into a fasta file. Uses the pre computed
        header and the sequence

        :return: a string that is a complete fasta entry
        """
        text = ">{}\n{}\n".format(self.name, self.sequence)
        return text

class GenbankFile:
    """
    An Object for easily disecting a genbank file into the individual proteins
    """
    def __init__(self, file, file_text = None, protein_range=None, allowed_proteins=None):
        """
        :param file: a path to a genbank file
        :param file_text: an optional parameter that will be used as the file text
        instead of the read version of the file.
        :param protein_range: an optional range in which proteins can be selected
        from the file instead of all proteins
        :param allowed_proteins: a list of protein IDs that can be selected from
        the genbank file instead of all proteins
        """
        logging.debug("Started parsing genbank file {}.".format(file))
        if file_text == None:
            self.file = file
            file_text = self.__read_genbank_file(file)
        self.contigs = self.__create_contigs(file_text, protein_range, allowed_proteins)
        self.proteins = self.__list_proteins()

        logging.debug("Finished parsing genbank file {}.".format(file))

    def __list_proteins(self):
        protein_dict = {}
        for contig in self.contigs:
            prots = contig.proteins.values()
            for prot in prots:
                if prot.name in protein_dict:
                    logging.warning(
                        "Double fasta entry '{}' in file '{}'. Skipping...".format(
                            entry.protein.name, self.file))
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
            logging.critical("Invalid genbank file {}. Exiting...".format(self.file))
            raise MultiGeneBlastException("Invalid genbank file {}.".format(self.file))

        # make sure to remove potential old occurances of \r. Acts like a \n
        file_text = file_text.replace("\r", "\n")

        # do a basic check to see if the genbank file is valid
        if "     CDS             " not in file_text or "\nORIGIN" not in file_text:
            logging.critical("Genbank file {} is not properly formatted or contains no sequences".format(self.file))
            raise MultiGeneBlastException("Genbank file {} is not properly formatted or contains no sequences".format(self.file))
        return file_text

    def __create_contigs(self, text, protein_range, allowed_proteins):
        contigs = []
        for cont in text.split("//\n")[:-1]:
            contig = Contig(cont,protein_range, allowed_proteins)
            contigs.append(contig)
        return contigs



class Contig:
    def __init__(self, file_text, protein_range=None, allowed_proteins=None):
        gene_defenitions, dna_sequence = file_text.split("\nORIGIN")
        # Extract DNA sequence and calculate complement of it
        dna_sequence = clean_dna_sequence(dna_sequence)
        c_dna_sequence = complement(dna_sequence)

        # Extract gene information
        gene_information = gene_defenitions.split("     CDS             ")

        #read general information
        self.accession = ""
        self.definition = ""
        self.__extract_header(gene_information[0])

        #exract all the entries from the genbank file
        self.entries = self.__extract_entries(gene_information[1:], dna_sequence, c_dna_sequence, protein_range, allowed_proteins)
        self.proteins = self.__get_unique_proteins()


    def __get_unique_proteins(self):
        """
        Make sure that you retrieve a unique list of named proteins after names
        have been cleaned for illegal characters

        :return: a dictionary of protein objects
        """
        protein_dict = {}
        for entry in self.entries:
            if entry.protein.name in protein_dict:
                logging.warning("Double fasta entry '{}' in file '{}'. Skipping...".format(entry.protein.name, self.file))
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
            g = GenbankEntry(gene, index + 1, dna_seq, c_dna_seq, self.accession)
            if protein_range != None and g.protein.start >= protein_range[0] and g.protein.stop <= protein_range[1]:
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
        Disect the header of the genbank file and extract the accesion and
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
                for def_line in header_lines[index + 2:]:
                    #end of definition line
                    if not def_line.startswith("    "):
                        break
                    self.definition += def_line.strip()
        # Test if accession number is probably real GenBank/RefSeq acc nr
        if testaccession(self.accession) == "n":
            logging.debug("Probably invalid GenBank/Refseq accesion found for {}.".format(self.file))
            self.accession = ""
        if self.accession == "":
            logging.debug("No valid accesion found for file {}".format(self.file))
        if self.definition == "":
            logging.debug("No definition found for file {}".format(self.file))


class GenbankEntry:
    """
    Disect a genbank entry and substract allot of potential values.
    """
    def __init__(self, entry_string, nr, dna_seq, c_dna_seq, file_accession):
        """
        :param entry_string: a string that contains a CDS genbank entry
        :param nr: a unique integer that can be used to create a genename is no
        name is available
        :param dna_seq: the dna sequence present in the genbank file
        :param c_dna_seq: the complement sequence present in the genbank file
        :param file_accession: the genbank file accession this entry originated
        from
        """
        #a list of locations
        self.location = None
        self.file_accession = file_accession

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
            if line.lower().startswith("codon_start"):
                self.__codon_start = int(self.__clean_line(line))
            elif line.lower().startswith("protein_id"):
                self.protein_id = self.__clean_line(line)
            elif line.lower().startswith("locus_tag"):
                self.locus_tag = self.__clean_line(line)
            elif line.lower().startswith("gene"):
                self.gene_name = self.__clean_line(line)
            elif line.lower().startswith("translation"):
                self.protein_code += self.__clean_line(line)
                for prot_line in cds_part[index + 2:]:
                    self.protein_code += prot_line.strip().replace('"', "")
                break
            elif line.lower().startswith("product"):
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
        string = string[5:-1]
        individual_strings = string.split(",")
        locations = []
        complement = False
        for s in individual_strings:
            if "complement" in s:
                complement = True
            locations.append(self.__convert_location(s[11:-1]))
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
        locations.sort()

        #if a sequence is available use that, otherwise extract it
        if self.protein_code:
            sequence = self.protein_code
        else:
            sequence = ""
            for loc in self.locations:
                #in case of invalid location skip to the next one.
                if loc[0] == None:
                    continue
                start, end = loc
                if reverse_complement:
                    nc_seq = c_dna_seq[(start - 1) - self.__codon_start - 1 : end]
                else:
                    nc_seq = dna_seq[(start - 1) - self.__codon_start - 1 : end]
                sequence += nc_seq
            sequence = translate(sequence)

        #the first location and then the start
        start = locations[0][0]
        return Protein(sequence, start, self.gene_name, strand, start_header="input|c1",
                       locus_tag=self.locus_tag, annotation=self.annotation, genbank_file=self.file_accession,
                       protein_id=self.protein_id)


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

def complement(seq):
    """
    Create the complement strand of a sequence. Use replacing for speed

    :param seq: a DNA sequence as a string
    :return: the complement of that string
    """
    seq = seq.replace("a", "x").replace("A", "X")
    seq = seq.replace("t", "a").replace("T", "A")
    seq = seq.replace("x", "t").replace("X", "T")

    seq = seq.replace("c", "x").replace("C", "X")
    seq = seq.replace("g", "c").replace("G", "C")
    seq = seq.replace("x", "g").replace("X", "G")
    return seq

def testaccession(accession):
    #TODO look into this function at some time
    #Test if accession number is probably real GenBank/RefSeq acc nr
    numbers = list(range(0,10))
    letters = []
    for i in ascii_letters:
        letters.append(i)
    nrnumbers = 0
    nrletters = 0
    for i in accession:
        if i in letters:
            nrletters += 1
        try:
            j = int(i)
            if j in numbers:
                nrnumbers += 1
        except:
            pass
    test = "y"
    if nrnumbers < 3 or nrletters < 1:
        test = "n"
    return test

def remove_illegal_characters(string):
    """
    Remove any character from ILLEGAL_CHARACTERS from string

    :param string: a string
    :return: the string without illegal characters. If no illegal characters
    are found the originall string is returned
    """
    clean_string = ""
    for letter in string:
        if letter not in ILLEGAL_CHARACTERS:
            clean_string += letter
    return clean_string