Usage of format_embl.py
- Prepare a tab-delimited TXT file ‘annotationtable.txt’ containing a table with the information on each exon / gene structured in the following columns: 1) name of contig FASTA file, 2) gene locus tag (must be unique), 3) 5’ exon or gene start, 4) 3’ exon or gene end, 5) gene annotation.
- Copy ‘annotationtable.txt’ and the FASTA files of the contigs into the folder in which you saved the EMBL formatting script.
- In the command line, type ‘python format_embl.py’
- If the script finishes successfully, input EMBL files for each contig will have been created which you can use as input for antiSMASH.
- If you have many contigs resulting in many EMBL files, you can use the batch-processing script available from the download page on the antiSMASH web site.

Example files for Aspergillus fumigatus have been provided.

Copyright: Marnix Medema, University of Groningen.