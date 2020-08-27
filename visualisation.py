

#imports
import logging
from collections import OrderedDict
from pysvg.filter import *
from pysvg.gradient import *
from pysvg.linking import *
from pysvg.script import *
from pysvg.shape import *
from pysvg.structure import *
from pysvg.style import *
from pysvg.text import *
from pysvg.builders import *

from constants import *
from utilities import is_valid_accession
from utilities import MultiGeneBlastException

#constant specifically required for mbg
MGBPATH = get_mgb_path()

class ClusterCollectionSvg:
    """
    Contains a list of all the clusters in the collection as well as an image
    that is the collection of clusters.
    """
    def __init__(self ,query_cluster, clusters, gene_color_dict, user_options):
        """
        :param query_cluster: a Cluster object containing the user defined cluster
        :param clusters: a list of Cluster objects
        :param gene_color_dict: a dictionary linking protein names to rgb colors
        :param user_options: an option object with user defined options
        """
        self.width = SCREENWIDTH
        self.dpn = self.__calculate_distance_per_nucleotide \
            (clusters + [query_cluster])
        self.cluster_svgs = self.__configure_cluster_svgs(clusters, query_cluster, gene_color_dict, user_options)
        self.cluster_image = self.__create_image()

    @property
    def XML(self):
        """
        Create an XML string that can be used in an XHTML file for visuals

        :return: a string
        """
        return self.cluster_image.getXML()

    def __create_image(self):
        """
        Create an image from all the elements saved in the self.cluster_svgs
        value

        :return: an svg object
        """
        svg_image = svg(x=0, y=0, width=int(self.width * 0.75), height=2770)
        viewbox = "0 0 " + str(int(self.width * 0.8)) + " " + str(2680)
        svg_image.set_viewBox(viewbox)
        svg_image.set_preserveAspectRatio("none")
        for cluster_svg in self.cluster_svgs.values():
            for group in cluster_svg.cluster_image_groups:
                svg_image.addElement(group)
        return svg_image

    def __calculate_distance_per_nucleotide(self, all_clusters):
        """
        Calculate a distance per nucleotide for the biggest cluster in
        all_clusters. This is a fraction that needs to be applied to start and
        stop locations to scale them to the size of the display window

        :param all_clusters: a list of Cluster objects including the query cluster
        :return: a float
        """
        # biggest core size
        biggest_size = 0
        for cluster in all_clusters:
            start, stop = cluster.core_start_stop()
            size = (stop - start) + 2 * SVG_CORE_EXTENSION
            if size > biggest_size:
                biggest_size = size
        return (self.width * 0.75) / biggest_size

    def __configure_cluster_svgs(self, clusters, query_cluster, gene_color_dict, user_options):
        """
        Create ClusterSvg objects for all clusters and the query

        :param query_cluster: the cluster that is the query for the mgb run
        :param clusters: a list of clusters that have to be presented with the
        query cluster
        :param gene_color_dict: a dictionary that links proteins with blast hits to
        colors so all proteins with the same blast subject have the same color
        :param user_options: an Option object with user defined options
        :return: a dictionary of ClusterSvg objects
        """
        cluster_svgs = OrderedDict()

        line = True
        if user_options.architecture_mode:
            line = False
        cluster_svgs["query"] = ClusterSvg \
            (list(query_cluster.proteins.values()),
                                           self.dpn, 0, (self.width * 0.75),
                                           gene_color_dict, line=line)

        query_cluster_size = query_cluster.core_start_stop()[1] - query_cluster.core_start_stop()[0]
        for index, cluster in enumerate(clusters):
            # make sur the cluster is sorted. This should be the case but calling
            # this method does no harm when the cluster is sorted already
            cluster.sort()

            # take an a good size sample of proteins
            reduced_proteins = self.__reduced_cluster_proteins(cluster, query_cluster_size)

            # check if the orientation of a cluster is inverted compared to the query
            # if so reverse the proteins.
            reverse = False
            if not self.__equal_to_query_strand_orientation(cluster):
                reverse = True
            cluster_svgs[cluster.no] = ClusterSvg(reduced_proteins, self.dpn,
                                                  index + 1, (self.width * 0.75),
                                                  gene_color_dict, reverse=reverse)
        return cluster_svgs

    def __reduced_cluster_proteins(self, cluster, query_cluster_size):
        """
        Get a list of proteins of the cluster that contains all proteins with blast
        hits and not to many of the surrounding proteins. The size has to be minimal
        0.8 * that of the query if the size of the cluster allows that.

        :param cluster: a Cluster object
        :param query_cluster_size: the size of the query cluster
        :return: a list of Protein objects that are around the blast hits in the
        cluster
        """
        # reduce the amount of proteins in the cluster by ignoreing proteins that
        # are to far from blast hits
        cluster_proteins = list(cluster.proteins.values())
        core_start, core_stop = cluster.core_start_stop()
        core_size = core_stop - core_start
        if query_cluster_size * 0.8 > core_size:
            extra_distance = (query_cluster_size * 0.8 - core_size) / 2
            core_start = max(core_start - extra_distance,
                             cluster.start_stop()[0])
            core_stop = min(core_stop + extra_distance,
                            cluster.start_stop()[1])
        min_start = core_start - SVG_CORE_EXTENSION
        max_end = core_stop + SVG_CORE_EXTENSION

        reduced_proteins = []
        for prot in cluster_proteins:
            if prot.start > min_start and prot.stop < max_end:
                reduced_proteins.append(prot)
        return reduced_proteins

    def __equal_to_query_strand_orientation(self, cluster):
        """
        Check if the cluster orientation is the same as the query orientation. This
        is done by checking for each blast hit in the cluster if it has the same
        orientation as the protein in blasted against.

        :param cluster: a Cluster object
        :return: a Boolean
        """
        equal_score = 0
        for protein_name, result_set in cluster.blast_hit_proteins.items():
            for blast_result in result_set:
                # compare the query protein strand with that of the cluster
                if blast_result.query_protein.strand == cluster.get_protein(protein_name).strand:
                    equal_score += 1
                else:
                    equal_score -= 1
        return equal_score >= 0


class ClusterSvg:
    """
    An svg of an individual cluster containing one or more proteins
    """
    def __init__(self, proteins, dpn, index, screenwidth, color_dict, line = True, reverse=False):
        """
        :param proteins: a list of Protein objects
        :param dpn: stands for distance per nucleotide. This is the distance on the
        svg object per nucleotide on the protein
        :param index: the position of the cluster on the current page
        :param screenwidth: the width the screen
        :param color_dict: a dictionary that links proteins with blast hits to
        colors so all proteins with the same blast subject have the same color
        :param line: a boolean that tells if a lign should be drawn trough all the
        gene shapes
        :param reverse: if the cluster is reversed compared to the query
        """
        self.reversed = reverse
        self.dpn = dpn
        self.start = int(proteins[0].start * self.dpn)
        self.stop = int(proteins[-1].stop * self.dpn)

        self.start_height = int(40 + HITS_PER_PAGE * index)
        self.total_width = screenwidth
        if reverse:
            start_offset = - (self.total_width - (self.start - self.stop)) / 2
        else:
            start_offset = (self.total_width - (self.stop - self.start)) / 2

        # gene arrow object saved by protein names
        self.gene_arrows = self.__configure_gene_arrows(proteins, start_offset, color_dict)

        # the final complete image of all gene arrow objects
        self.cluster_image_groups = self._cluster_image(color_dict, line)

    def __configure_gene_arrows(self, proteins, start_offset, color_dict):
        """
        Create GeneArrow objects for all proteins.

        :param proteins: a list of Protein objects
        :param start_offset: an offset from the start of the image to center
        the cluster on the line
        :param color_dict:
        :return:
        """
        arrows = OrderedDict()
        for prot in proteins:
            if prot.name in color_dict:
                color = color_dict[prot.name]
            else:
                color = "#FFFFFF"
            arrows[prot.name] = GeneArrow(prot, start_offset, self.start, self.dpn, color, self.start_height, reverse=self.reversed)
        return arrows

    def _cluster_image(self, color_dict, line):
        """
        Create groups of shapes that together form a protein cluster represenataion

        :param color_dict: a dictionary that links proteins with blast hits to
        colors so all proteins with the same blast subject have the same color
        :param line: a boolean that tells if a lign should be drawn trough all the
        gene shapes
        :return: a list of groups that form an image of a cluster
        """
        groups = []
        builder = ShapeBuilder()
        group = g()
        if line:
            group.addElement \
                (builder.createLine(10, self.start_height - 5, self.total_width,
                                                self.start_height - 5, strokewidth=1, stroke="grey"))
        else:
            # create an invisible line to allow the xhtml to properly work but for
            # it not to show
            group.addElement \
                (builder.createLine(10, self.start_height - 5, self.total_width,
                                                self.start_height - 5, strokewidth=1
                                   ,stroke="white"))
        groups.append(group)

        for gene_arrow in self.gene_arrows.values():
            groups.append(gene_arrow.arrow)

        return groups


class GeneArrow:
    """
    A polygon that represents an arrow that has a proportional size to the
    window size
    """
    ARROW_HEIGHT = 10
    def __init__(self, protein, start_offset, cluster_start, dpn, color, start_height, reverse = False):
        """
        :param protein: a Protein object that the arrow should be made of
        :param start_offset: a x offset
        :param cluster_start: the scaled start of the first protein in the cluster
        :param dpn: stands for distance per nucleotide. This is the distance on the
        svg object per nucleotide on the protein
        :param color: the rgb color of the protein
        :param start_height: a y offset
        :param reverse: if the cluster of the gene is reversed or not
        """
        self.start = int \
            (abs(start_offset + (protein.start * dpn - cluster_start)))
        self.stop = int \
            (abs(start_offset + (protein.stop * dpn - cluster_start)))

        # configure strand
        self.strand = self.__configure_strand(protein.strand, reverse)
        self.color = color

        self.arrow = self._gene_arrow(start_height, self.ARROW_HEIGHT)

    def __configure_strand(self, prot_strand, reverse):
        """
        Flips the strand of if reverse is True

        :param prot_strand: a strand orientation '+' or '-'
        :param reverse: a boolean
        :return: a strand orientation '+' or '-'
        """
        if reverse:
            if prot_strand == "+":
                return "-"
            else:
                return "+"
        return prot_strand

    def _gene_arrow(self, start_height, height):
        """
        Create a gene arrow using a polygon

        :param start_height:
        :param height: height of the arrow
        :return:
        """
        # TODO look to clean up this function if there is time left
        halfheight = height / 2
        if self.start > self.stop:
            start2 = self.stop
            stop2 = self.start
            self.start = start2
            self.stop = stop2
        builder = ShapeBuilder()
        if (self.stop - self.start) < halfheight:
            if (self.strand == "+"):
                points_as_tuples = [(self.start, start_height),
                                    (self.stop, start_height - halfheight),
                                    (self.start, start_height - height),
                                    (self.start, start_height)
                                    ]
            if (self.strand == "-"):
                points_as_tuples = [(self.start, start_height - halfheight),
                                    (self.stop, start_height - height),
                                    (self.stop, start_height),
                                    (self.start, start_height - halfheight)
                                    ]
        else:
            if (self.strand == "+"):
                arrowstart = self.stop - halfheight
                points_as_tuples = [(self.start, start_height),
                                    (arrowstart, start_height),
                                    (self.stop, start_height - halfheight),
                                    (arrowstart, start_height - height),
                                    (self.start, start_height - height),
                                    (self.start, start_height)
                                    ]
            if (self.strand == "-"):
                arrowstart = self.start + halfheight
                points_as_tuples = [(self.start, start_height - halfheight),
                                    (arrowstart, start_height - height),
                                    (self.stop, start_height - height),
                                    (self.stop, start_height),
                                    (arrowstart, start_height),
                                    (self.start, start_height - halfheight)
                                    ]
        pg = builder.createPolygon \
            (points=builder.convertTupleArrayToPoints(points_as_tuples),
                                   strokewidth=1 ,stroke='black', fill=self.color)
        return pg


def create_xhtml_template(html_parts, page_indx, page_sizes):
    """
    Open the xhtml file and load some pre written information into the file

    :param html_parts: a list of prewriten code to put in the xhtml file at the
    beginning
    :param page_sizes: a list of integers with the size of all pages
    :param page_indx: the index of the current page in the page_sizes list
    :return: an opened file object
    """

    html_outfile = open("displaypage{}.xhtml".format(page_indx + 1),"w")
    html_outfile.write('{}  displaycblastresults("{}","all"){}'.format(html_parts[0], page_indx + 1, html_parts[1]))
    inner_string = "var list=[{},'all'];".format(",".join([str(nr + page_indx * 50 + 1) for nr in range(page_sizes[page_indx])]))
    html_outfile.write(inner_string + html_parts[2])
    html_outfile.write(inner_string + html_parts[3])
    html_outfile.write('<a class="bigtext"><br/><br/>&nbsp;Results pages: ')
    for page_nr in [nr + 1 for nr in range(len(page_sizes))]:
        html_outfile.write('<a href="displaypage{}.xhtml" class="bigtext">{}</a>'.format(page_nr, page_nr))
        if page_nr != len(page_sizes):
            html_outfile.write(", ")
        else:
            html_outfile.write("</a>")
    html_outfile.write(html_parts[4])
    return html_outfile

def write_xhtml_output(html_outfile, clusters, query_cluster, page_indx, page_size, user_options, svg_images, gene_color_dict):
    """
    Write ClusterBlast divs with pictures and description pop-up tags

    :param html_outfile: an open file object where the xhtml output is written
    into
    :param clusters: a list of Cluster objects retrieved by MultiGeneBlast
    :param query_cluster: a Cluster object that contains the query
    :param page_indx: the index of the current page that is being drawn
    :param page_size: the size of the current page (max = HITS_PER_PAGE)
    :param user_options: an Option object containing user defined options
    :param svg_images: A dictionary with classes containing information for
    :param gene_color_dict: a dictionary linking protein names to rgb colors
    creating all the svg images needed for the visual output.
    """
    screenwidth = SCREENWIDTH
    html_outfile.write('<div id="clusterblastview" class="clusterdescr">\n\n')
    #Add menu bar 3
    html_outfile.write('<div id="bartext3" style="color:#FFFFFF; font-size:1em; position:absolute; z-index:2; top:3px; left:20px;"><b>MultiGeneBlast hits</b></div>')
    html_outfile.write('<div id="descrbar3" style="position:absolute; z-index:1; top:0px;"><img src="images/bar.png" height="25" width="{}"/></div>'.format(int(screenwidth * 0.75)))
    html_outfile.write('<div class="help" id="help3" style="position:absolute; z-index:1; top:2px; left:{}px;"><a href="http://multigeneblast.sourceforge.net/usage.html"'
                      ' target="_blank"><img border="0" src="images/help.png"/></a></div>'.format(int(screenwidth * 0.75 -30)))
    page_nr = page_indx + 1
    html_outfile.write('<div id="qcluster{}">\n<br/><br/>\n<div align="left">\n<form name="clusterform{}">\n<select name="selection{}" onchange="javascript:navigate(this);">\n'.format(page_nr, page_nr, page_nr))
    html_outfile.write('<option value="">Select gene cluster alignment</option>\n')
    for index, cluster in enumerate(clusters):
        cluster_summary = "Cluster{}: {}".format(cluster.no, ",".join(list(cluster.proteins.keys())[:3]))
        if len(cluster_summary) > 50:
            cluster_summary = "{}...".format(cluster_summary[:47])
        html_outfile.write('<option value="javascript:displaycblastresults({},{})">{}</option>\n'.format(page_nr, index + 1 + page_indx * HITS_PER_PAGE, cluster_summary))
    html_outfile.write('</select>\n</form>\n\n</div>')
    html_outfile.write('<div style="position:absolute; top:33px; left:{}px;"><img src="images/button.gif" name="button{}" onclick="javascript:displaybutton({})'
                       ';"/></div>'.format(screenwidth * 0.625, page_nr, page_nr))
    for index, cluster in enumerate(clusters):
        html_outfile.write('<div id="hitcluster{}_{}">\n'.format(page_nr, index + 1 + page_indx * HITS_PER_PAGE))

        #Load svg and embed it into XHTML
        svg_image = svg_images["clusterblast_{}".format(str(cluster.no))]
        names = list(svg_image.cluster_svgs.keys())
        svg_lines = svg_image.XML.split("\n")
        html_outfile.write('\n{}id="svg{}_{}" >\n'.format(svg_lines[0][:-1], page_nr, index + 1 + page_indx * HITS_PER_PAGE))
        name_index = -1
        protein_index = 0
        for svg_line in svg_lines[1:]:
            if svg_line.startswith("</g>"):
                name_index += 1
                protein_index = 0
                html_outfile.write(svg_line + "\n")
            elif svg_line.startswith("<polygon"):
                if names[name_index] == "query":
                    html_outfile.write('<g id="{}{}_{}_{}"  >\n'.format("q", page_nr, page_indx * HITS_PER_PAGE + index + 1, protein_index))
                else:
                    html_outfile.write('<g id="{}{}_{}_{}"  >\n'.format("h", page_nr, page_indx * HITS_PER_PAGE + index + 1, protein_index))
                html_outfile.write(svg_line + "\n")
                html_outfile.write('</g>\n')
                protein_index += 1
            else:
                html_outfile.write(svg_line + "\n")

        #Insert gene cluster descriptions
        contig_desc = cluster.contig_description
        if len(contig_desc) > 90:
            contig_desc = "{}...".format(contig_desc[:87])
        if is_valid_accession(cluster.contig):
            cluster_desc = 'Cluster{}: <a href="http://www.ncbi.nlm.nih.gov/nuccore/{}" target="_blank"> {}</a> {}&nbsp;&nbsp;&nbsp;&nbsp;Total' \
                           ' score: {:.1f}&nbsp;&nbsp;&nbsp;&nbsp; Cumulative Blast bit score {:.2f}'.format(cluster.no, cluster.contig, cluster.contig,contig_desc,
                                                                                cluster.score,cluster.segmented_score("accumulated blast score")[0] * 1_000_000)
        else:
            cluster_desc = 'Cluster{}: {} {}&nbsp;&nbsp;&nbsp;&nbsp;Total score: {:.1f}&nbsp;&nbsp;&nbsp;&nbsp;' \
                           ' Cumulative Blast bit score: {:.2f}'.format(cluster.no, cluster.contig, contig_desc, cluster.scorecluster,
                                cluster.segmented_score("accumulated blast score")[0] * 1_000_000)
        query_desc = query_cluster.contig_description
        if len(query_desc) < 90:
            query_desc = "Query: {}".format(query_desc)
        else:
            query_desc = "Query: {}...".format(query_desc[:87])
        html_outfile.write('<div id="descriptionquery" style="text-align:left; position:absolute; top:60px;'
                           ' left:10px; font-size:10px; font-style:italic">{}</div>\n'.format(query_desc))
        html_outfile.write('<div id="description{}" style="text-align:left; position:absolute; top:115px;'
                           ' left:10px; font-size:10px; font-style:italic">{}</div>\n'.format(page_nr, cluster_desc))

        #Insert NCBI links
        html_outfile.write('<div id="pub_pics" style="position:absolute; top:175px; left:0px; font-size:10px"> Hit cluster cross-links: \n')
        html_outfile.write('&nbsp;&nbsp;<a href="http://www.ncbi.nlm.nih.gov/nuccore/{}"'
                           ' target="_blank"><img align="absmiddle" border="0" src="images/genbank.gif"/></a>\n'.format(cluster.contig))
        html_outfile.write('</div>\n\n')

        #Create gene pop-ups
        query_starts = [arrow.start for arrow in svg_image.cluster_svgs["query"].gene_arrows.values()]
        query_ends = [arrow.stop for arrow in svg_image.cluster_svgs["query"].gene_arrows.values()]
        for pindex, protein in enumerate(query_cluster.proteins.values()):
            html_outfile.write('<div id="q{}_{}_{}_div" class="hidden popup" style="position:absolute; top:100px; left:{}px;">\n'.format(page_nr, page_indx * HITS_PER_PAGE + index + 1 ,pindex, int(query_starts[pindex] * 0.875)))
            html_outfile.write("{}\n".format(protein.annotation.replace("_", " ")))
            if protein.protein_id != "" and is_valid_accession(protein.protein_id):
                html_outfile.write('<br/>Accession: <a href="http://www.ncbi.nlm.nih.gov/protein/{}" target="_blank">{}</a>\n'.format(protein.protein_id, protein.protein_id))
            html_outfile.write("<br/>Location:{}-{}\n".format(protein.start, protein.stop))
            html_outfile.write("</div>\n\n")
            html_outfile.write('<div id="q{}_{}_{}_divtext" class="hidden genenames" style="position:absolute; top:75px; left:'
                               '{}px;">\n'.format(page_nr, page_indx * HITS_PER_PAGE + index + 1 ,pindex, ((query_starts[pindex] + query_ends[pindex]) / 2) * 0.9375))
            if protein.locus_tag != "":
                html_outfile.write(protein.locus_tag)
            else:
                html_outfile.write(protein.protein_id)
            html_outfile.write("</div>\n\n")

        #all the protein names that are actually visually displayed
        cluster_protein_names = svg_image.cluster_svgs[cluster.no].gene_arrows.keys()
        prot_starts = [arrow.start for arrow in svg_image.cluster_svgs[cluster.no].gene_arrows.values()]
        prot_ends = [arrow.stop for arrow in svg_image.cluster_svgs[cluster.no].gene_arrows.values()]
        for pindex, protein_name in enumerate(cluster_protein_names):
            protein = cluster.get_protein(protein_name)
            html_outfile.write('<div id="h{}_{}_{}_div" class="hidden popup" style="position:absolute; top:151px; left:{}px;">\n'.format(page_nr, page_indx * HITS_PER_PAGE + index + 1 ,pindex, int(prot_starts[pindex] * 0.875)))
            html_outfile.write("{}\n".format(protein.annotation.replace("_"," ")))
            link = "http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins&amp;PROGRAM=blastp&amp;BLAST_PROGRAMS=blastp&amp;QUERY={}&amp;LINK_LOC=protein&amp;PAGE_TYPE=BlastSearch".format(protein.protein_id)
            if user_options.dbtype == "nucl":
                html_outfile.write('<br/>Accession: <a href="http://www.ncbi.nlm.nih.gov/nuccore/{}" target="_blank">{}</a>\n'.format(protein.protein_id, protein.protein_id))
            elif protein.locus_tag != "" and is_valid_accession(protein.locus_tag):
                html_outfile.write('<br/>Accession: <a href="http://www.ncbi.nlm.nih.gov/protein/{}" target="_blank">{}</a>\n'.format(protein.locus_tag, protein.locus_tag))
            html_outfile.write("<br/>Location: {}-{}\n".format(protein.start, protein.stop))
            if protein_name in cluster.blast_hit_proteins:
                blast_list = sorted(list(cluster.blast_hit_proteins[protein_name]), key=lambda x: (x.evalue, -x.bit_score))
                for blast_result in blast_list:
                    html_outfile.write("<br/><br/><b>BlastP hit with {}</b>\n<br/>Percentage identity: {} %\n".format(blast_result.query_protein.protein_id, blast_result.percent_identity))
                    html_outfile.write("<br/>BlastP bit score: {}\n<br/>Sequence coverage: {:.1f} %\n".format(blast_result.bit_score, blast_result.percent_coverage))
                    html_outfile.write("<br/>E-value: {}\n<br/>".format(blast_result.evalue))
            if is_valid_accession(protein.protein_id) and user_options.dbtype != "nucl":
                html_outfile.write('<br/><a href="{}" target="_blank"> NCBI BlastP on this gene </a>\n'.format(link))
            #muscle allignment
            if protein.name in gene_color_dict:
                color = gene_color_dict[protein.name].replace("#","")
                path = "{}{}muscle_info{}group{}_muscle.fasta".format(user_options.outdir, os.sep, os.sep, color)
                if os.path.exists(path):
                    html_outfile.write("<br/><a href=\"{}\" target=\"_blank\"> Muscle alignment of this gene with homologs </a>\n".format(path))
            html_outfile.write("</div>\n\n")
            html_outfile.write('<div id="h{}_{}_{}_divtext" class="hidden genenames" style="position:absolute; top:126px; left:'
                               '{}px;">\n'.format(page_nr, page_indx * HITS_PER_PAGE + index + 1 ,pindex,((prot_starts[pindex] + prot_ends[pindex]) / 2) * 0.9375))
            if protein.locus_tag != "":
                html_outfile.write(protein.locus_tag)
            else:
                html_outfile.write(protein.protein_id)
            html_outfile.write("</div>\n\n")
        html_outfile.write('</div>\n')
    #writing the full page with all the proteins
    all_svg_image = svg_images["clusterblast_page{}_all".format(page_indx + 1)]
    names = list(all_svg_image.cluster_svgs.keys())
    html_outfile.write('<div id="hitcluster{}_all" style="display:none">\n'.format(page_nr))
    all_svg_lines = all_svg_image.XML.split("\n")
    html_outfile.write('\n{}id="svg{}_all" >\n'.format(all_svg_lines[0][:-1], page_nr))
    name_index = -1
    protein_index = 0
    for svg_line in all_svg_lines[1:]:
        if svg_line.startswith("</g>"):
            name_index += 1
            protein_index = 0
            html_outfile.write(svg_line + "\n")
        elif svg_line.startswith("<polygon"):
            if names[name_index] == "query":
                html_outfile.write('<g id="all_{}_0_{}"  >\n'.format(page_nr, protein_index))
            else:
                html_outfile.write('<g id="all_{}_{}_{}"  >\n'.format(page_nr, page_indx * HITS_PER_PAGE + name_index, protein_index))
            html_outfile.write(svg_line + "\n")
            html_outfile.write('</g>\n')
            protein_index += 1
        else:
            html_outfile.write(svg_line + "\n")

    query_desc = query_cluster.contig_description
    if len(query_desc) < 90:
        query_desc = "Query: {}".format(query_desc)
    else:
        query_desc = "Query: {}...".format(query_desc[:87])
    html_outfile.write('<div id="descriptionquery" style="text-align:left; position:absolute; top:60px; left:10px; font-size:10px; font-style:italic">{}</div>\n'.format(query_desc))

    #make sure these are added first so the highlights dont go under the text
    for index, cluster in enumerate(clusters):
        # Insert gene cluster descriptions
        contig_desc = cluster.contig_description
        if len(contig_desc) > 90:
            contig_desc = "{}...".format(contig_desc[:87])
        if is_valid_accession(cluster.contig):
            cluster_desc = 'Cluster{}: <a href="http://www.ncbi.nlm.nih.gov/nuccore/{}" target="_blank"> {}</a> {}&nbsp;&nbsp;&nbsp;&nbsp;Total' \
                           ' score: {:.1f}&nbsp;&nbsp;&nbsp;&nbsp; Cumulative Blast bit score {:.2f}'.format(cluster.no, cluster.contig, cluster.contig,contig_desc,
                                                                                cluster.score,cluster.segmented_score("accumulated blast score")[0] * 1_000_000)
        else:
            cluster_desc = 'Cluster{}: {} {}&nbsp;&nbsp;&nbsp;&nbsp;Total score: {:.1f}&nbsp;&nbsp;&nbsp;&nbsp;' \
                           ' Cumulative Blast bit score: {:.2f}'.format(cluster.no, cluster.contig, contig_desc, cluster.scorecluster,
                                cluster.segmented_score("accumulated blast score")[0] * 1_000_000)
        html_outfile.write('<div id="description{}" style="text-align:left; position:absolute; top:{}px; left:10px; font-size:10px;'
                           ' font-style:italic">{}</div>\n'.format(page_nr, int(63 + (51.7 * (index + 1))), cluster_desc))
        #add query for the first index
    index = 0
    for index, cluster in enumerate(clusters):
        if index == 0:
            query_starts = [arrow.start for arrow in all_svg_image.cluster_svgs["query"].gene_arrows.values()]
            query_ends = [arrow.stop for arrow in all_svg_image.cluster_svgs["query"].gene_arrows.values()]
            for pindex, protein in enumerate(query_cluster.proteins.values()):
                html_outfile.write('<div id="all_{}_0_{}_div" class="hidden popup" style="position:absolute; top:100px; left:{}px;">\n'.format(page_nr ,pindex, int(query_starts[pindex] * 0.875)))
                html_outfile.write("{}\n".format(protein.annotation.replace("_", " ")))
                if protein.protein_id != "" and is_valid_accession(protein.protein_id):html_outfile.write('<br/>Accession: <a href="http://www.ncbi.nlm.nih.gov/protein/{}"'
                                                                                                     ' target="_blank">{}</a>\n'.format(protein.protein_id, protein.protein_id))
                html_outfile.write("<br/>Location:{}-{}\n".format(protein.start,protein.stop))
                html_outfile.write("</div>\n\n")
                html_outfile.write('<div id="all_{}_0_{}_divtext" class="hidden genenames" style="position:absolute; top:75px; left:''{}px;'
                                   '">\n'.format(page_nr, pindex, ((query_starts[pindex] + query_ends[pindex]) / 2) * 0.9375))
                if protein.locus_tag != "":
                    html_outfile.write(protein.locus_tag)
                else:
                    html_outfile.write(protein.protein_id)
                html_outfile.write("</div>\n\n")
        cluster_protein_names = all_svg_image.cluster_svgs[cluster.no].gene_arrows.keys()
        prot_starts = [arrow.start for arrow in all_svg_image.cluster_svgs[cluster.no].gene_arrows.values()]
        prot_ends = [arrow.stop for arrow in all_svg_image.cluster_svgs[cluster.no].gene_arrows.values()]

        for pindex, protein_name in enumerate(cluster_protein_names):
            protein = cluster.get_protein(protein_name)
            html_outfile.write('<div id="all_{}_{}_{}_div" class="hidden popup" style="position:absolute; top:{}px; left:{}px;">\n'.format(page_nr, index + 1 + page_indx * HITS_PER_PAGE ,pindex, int(100 + (51.7 * (index + 1))), int(prot_starts[pindex] * 0.875)))
            html_outfile.write("{}\n".format(protein.annotation.replace("_", " ")))
            link = "http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins&amp;PROGRAM=blastp&amp;BLAST_PROGRAMS=blastp&amp;QUERY={}&amp;LINK_LOC=protein&amp;PAGE_TYPE=BlastSearch".format(protein.protein_id)
            if user_options.dbtype == "nucl":
                html_outfile.write('<br/>Accession: <a href="http://www.ncbi.nlm.nih.gov/nuccore/{}" target="_blank">{}</a>\n'.format(protein.protein_id, protein.protein_id))
            elif protein.locus_tag != "" and is_valid_accession(protein.locus_tag):
                html_outfile.write('<br/>Accession: <a href="http://www.ncbi.nlm.nih.gov/protein/{}" target="_blank">{}</a>\n'.format(protein.locus_tag, protein.locus_tag))
            html_outfile.write("<br/>Location: {}-{}\n".format(protein.start, protein.stop))
            if protein_name in cluster.blast_hit_proteins:
                blast_list = sorted(list(cluster.blast_hit_proteins[protein_name]), key=lambda x: (x.evalue, -x.bit_score))
                for blast_result in blast_list:
                    html_outfile.write("<br/><br/><b>BlastP hit with {}</b>\n<br/>Percentage identity: {} %\n".format(blast_result.query_protein.protein_id, blast_result.percent_identity))
                    html_outfile.write("<br/>BlastP bit score: {}\n<br/>Sequence coverage: {:.1f} %\n".format(blast_result.bit_score,blast_result.percent_coverage))
                    html_outfile.write("<br/>E-value: {}\n<br/>".format(blast_result.evalue))
            if is_valid_accession(protein.protein_id) and user_options.dbtype != "nucl":
                html_outfile.write('<br/><a href="{}" target="_blank"> NCBI BlastP on this gene </a>\n'.format(link))
            if protein.name in gene_color_dict:
                color = gene_color_dict[protein.name].replace("#","")
                path = "{}{}muscle_info{}group{}_muscle.fasta".format(user_options.outdir, os.sep, os.sep, color)
                if os.path.exists(path):
                    html_outfile.write("<br/><a href=\"{}\" target=\"_blank\"> Muscle alignment of this gene with homologs </a>\n".format(path))
            html_outfile.write("</div>\n\n")
            html_outfile.write('<div id="all_{}_{}_{}_divtext" class="hidden genenames" style="position:absolute; top:{}px; left:{}px;">\n'.format(page_nr, index + 1 + page_indx * HITS_PER_PAGE ,pindex, int(75 + (51.7 * (index + 1))), ((prot_starts[pindex]+prot_ends[pindex])/2)*0.9375))
            if protein.locus_tag != "":
                html_outfile.write(protein.locus_tag)
            else:
                html_outfile.write(protein.protein_id)
            html_outfile.write("</div>\n\n")
    html_outfile.write('</div>\n')
    html_outfile.write('</div>\n\n')
    html_outfile.write('</div>\n')
    html_outfile.write('<div id="creditsbar{}" class="banner" style="position:absolute; width:{}px; align:\'left\'; height:75; top:2750px; left:0px; color:#000066; z-index:-1;">'.format(index, int(0.98 * screenwidth)))
    html_outfile.write('<div style="float:center; font-size:0.9em;">\n<div style="position:absolute; top:0px; left:30px;">\n<img src="images/ruglogo.gif" border="0"/>&nbsp;&nbsp;&nbsp;&nbsp;\n<img src="images/gbblogo.gif" border="0"/>&nbsp;&nbsp;&nbsp;&nbsp;\n</div>\n<div style="position:absolute; top:10px; left:340px;">\nDetecting sequence homology at the gene cluster level with MultiGeneBlast.\n<br/>Marnix H. Medema, Rainer Breitling &amp; Eriko Takano (2013)\n<br/><i>Molecular Biology and Evolution</i> , 30: 1218-1223.\n</div>\n</div>\n</div>')

def create_xhtml_file(clusters, query_cluster, user_options, svg_images, page_sizes, page_indx, gene_color_dict):
    """
    Create an XHTML file for a set of clusters.

    :param clusters: a list of Cluster objects
    :param query_cluster: a Cluster objects with the user defined query genes
    :param user_options: an Option object with user defined options
    :param svg_images: a Dictionary of ClusterCollectionSvg objects, that contian
    information to create the svg images
    :param page_sizes: a list of integers with the size of all pages
    :param page_indx: the index of the current page in the page_sizes list
    :param gene_color_dict: a dictionary linking protein names to rgb colors
    """

    #Create HTML file with gene cluster info in hidden div tags for all pages
    try:
        with open(MGBPATH + os.sep + "visual_copys" + os.sep + "empty.xhtml", "r") as htmlfile:
            html = htmlfile.read()
            html = html.replace("\r", "\n")
            html_parts = html.split("<SPLIT HERE>")
    except:
        logging.critical("empty.xhml is missing from {}. Cannot create "
                         "the full html page. Exiting...".format(MGBPATH + os.sep + "visual_copys"))
        raise MultiGeneBlastException(
            "empty.xhml is missing from {}. Cannot create "
            "the full html page.".format(MGBPATH + os.sep + "visual_copys"))
    #note this returns an open file stream. Take care addign code in this loop
    html_outfile = create_xhtml_template(html_parts, page_indx, page_sizes)
    write_xhtml_output(html_outfile, clusters, query_cluster, page_indx, page_sizes[page_indx], user_options, svg_images, gene_color_dict)

    html_outfile.write(html_parts[-1])
    html_outfile.close()
