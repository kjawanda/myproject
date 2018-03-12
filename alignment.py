from Bio import SwissProt
from Bio import SeqIO
from Bio import AlignIO
from Bio import Phylo
from Bio import pairwise2
import gzip
import io
import subprocess
import sys
from Bio.pairwise2 import format_alignment
from Bio.SubsMat.MatrixInfo import blosum62
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align.Applications import MuscleCommandline

def format_conversion(file_name, current_format):
    """This function converts the file into the FASTA format. """
    name = os.path.splitext(file_name)[0]
    if os.path.splitext(file_name)[1] == ".fasta":
            print("This file is already in fasta format!")
    else:
        name = name + ".fasta"
        count = SeqIO.convert(file_name, current_format, name, "fasta")
        print("converted %i records" % count)

#format_conversion("cor6_6.gb", "genbank")
#format_conversion("opuntia.fasta", "fasta")

#Now that the files are in fasta format, it is important to decide whether pairwise alignment or multiple alignment needs to be used. The following functions are for pairwise alignment.
#Depending on whether the user wants global or local alignment, and the format of their sequences. If the sequences are in strings, then pairwise2_global or pairwise2_local is used.

def pairwise2_global(sequence1, sequence2):
    """This function aligns two sequences (string format) using global pairwise alignment. The algorithm used is Needleman-Wunsch."""
    alignments = pairwise2.align.globalms(sequence1, sequence2, 2, -1, -0.5, -0.1)
    for a in alignments:
        print(format_alignment(*a))

#pairwise2_global("ACT", "ACGT")

def pairwise2_global_file(sequence_file1, sequence_file2):
    """This function aligns two sequences (sequences in file) using global pairwise alignment. The algorithm used is Needleman-Wunsch."""
    seq1 = SeqIO.read(sequence_file1, "fasta")
    seq2 = SeqIO.read(sequence_file2, "fasta")
    alignments = pairwise2.align.globalds(seq1.seq, seq2.seq, blosum62, -10, -0.5)
    for a in alignments:
        print(format_alignment(*a))

#pairwise2_global_file("alpha.faa", "beta.faa")

def pairwise2_local(sequence1, sequence2):
    """This function aligns two sequences (string format) using local pairwise alignment. The algorithm used is Smith-Waterman."""
    alignments = pairwise2.align.localms(sequence1, sequence2, 2, -1, -0.5, -0.1)
    for a in alignments:
        print(format_alignment(*a))

#pairwise2_global("ACT", "ACGT")

def pairwise2_local_file(sequence_file1, sequence_file2):
    """This function aligns two sequences (sequences in file) using local pairwise alignment. The algorithm used is Smith-Waterman."""
    seq1 = SeqIO.read(sequence_file1, "fasta")
    seq2 = SeqIO.read(sequence_file2, "fasta")
    alignments = pairwise2.align.localds(seq1.seq, seq2.seq, blosum62, -10, -0.5)
    for a in alignments:
        print(format_alignment(*a))

#pairwise2_global_file("alpha.faa", "beta.faa")

#The following are the functions for multiple sequence alignment:

def clustalw(file):
    """ This is a function that aligns sequences within a file using clustalw.
    Then the sequences are represented in a phylogenetic tree.""" 
    cline = ClustalwCommandline("clustalw", infile=file)
    print(cline)
    stdout, stderr = cline()
    alignment_name = os.path.splitext(file)[0]
    #"{}.aln".format(alignment_name)
    alignment_name = alignment_name + ".aln"
    align = AlignIO.read(alignment_name,"clustal")
    print(align)
    print("The sequences can be represented as a phylogenetic tree: ")
    tree_name = os.path.splitext(file)[0]
    tree_name = tree_name + ".dnd"
    tree = Phylo.read(tree_name, "newick")
    Phylo.draw_ascii(tree) 

#clustalw("opuntia.fasta")

#For muscle, there are two types of functions, one for when the input is small and one where the input is large.

user_answer = input("Is the number of sequences for MSA less than 10?").lower().strip()
if user_answer == "yes":
    print("Please use the 'muscle_smallinput(file)' function where file is your file name")
    
if user_answer == "no":
    print("lets use muscle_large function:")
    print("Please use the 'muscle_smallinput(file)' function where file is your file name")

def muscle_smallinput(file):
    
    muscle_cline = MuscleCommandline(input=file)
    stdout, stderr = muscle_cline()
    muscle_align = AlignIO.read(io.StringIO(stdout), "fasta")
    print(muscle_align)

def muscle_largeinput(file):
    muscle_cline = MuscleCommandline(input=file)
    child = subprocess.Popen(str(muscle_cline),
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             universal_newlines=True,
                             shell=(sys.platform!="linux2"))
    muscle_align = AlignIO.read(child.stdout, "fasta")
    print(muscle_align)

#muscle_largeinput("opuntia.fasta")




