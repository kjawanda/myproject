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
from Bio.Align.Applications import TCoffeeCommandline
from Bio.Emboss.Applications import NeedleCommandline
from Bio.Emboss.Applications import WaterCommandline
from Bio.Emboss.Applications import NeedleallCommandline

handle = gzip.open("uniprot_sprot.xml.gz")
records = SeqIO.parse(handle, "uniprot-xml")

"""in order to do this, we first need to download clustal to the computer. we do
this by going into the downloads directory in the cline and typing
'sudo apt-get install clustalw'
i will be using opuntia.fasta, the downloaded file is found on http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc81
we first gonna try this and then move onto my coursework"""

def clustalw(file):
    cline = ClustalwCommandline("clustalw", infile=file)
    print(cline)
    stdout, stderr = cline()
    alignment_name = os.path.splitext(file)[0]
    #"{}.aln".format(alignment_name)
    alignment_name = alignment_name + ".aln"
    align = AlignIO.read(alignment_name,"clustal")
    print(align)
    tree_name = os.path.splitext(file)[0]
    tree_name = tree_name + ".dnd"
    tree = Phylo.read(tree_name, "newick")
    Phylo.draw_ascii(tree) 

clustalw("opuntia.fasta")

def muscle_smallinput(file):
    muscle_cline = MuscleCommandline(input=file)
    stdout, stderr = muscle_cline()
    muscle_align = AlignIO.read(io.StringIO(stdout), "fasta")
    print(muscle_align)

muscle_smallinput("opuntia.fasta")

#The approach used above is simple but if you are dealing with very large output text the fact that all of stdout and
#stderr is loaded into memory as a string can be a potential drawback.
#Using the 'subprocess' module we can work directly with handles instead: 

def muscle_largeinput(file):
    muscle_cline = MuscleCommandline(input=file)
    child = subprocess.Popen(str(muscle_cline),
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             universal_newlines=True,
                             shell=(sys.platform!="linux2"))
    muscle_align = AlignIO.read(child.stdout, "fasta")
    print(muscle_align)

muscle_largeinput("opuntia.fasta")


#this doesnt work at all!!!!
"""
def tcoffee():
    print("T-Coffee alignment:")
    tcoffee_cline = TCoffeeCommandline(infile="opuntia.fasta",
                                   output="clustalw",
                                   outfile="tcoffee_aligned.aln")

    tcoffee_aligned = AlignIO.read("tcoffee_aligned.aln", "clustal")
    print(align)





child = subprocess.Popen(str(tcoffee_cline),
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE,
                         universal_newlines=True,
                         shell=(sys.platform!="linux2"))
tcoffee_align = AlignIO.read(child.stdout, "fasta")
print(tcoffee_align)
"""
#The input for sequence1 and sequence 2 in pairwise2_global/pairwise2_local must be a sequence in a string.
#If the sequences are composed in files, use function pairwise2_global_file/pairwise2_local_file 

def pairwise2_global(sequence1, sequence2):
    #X = sequence1
    #Y = sequence2
    alignments = pairwise2.align.globalms(sequence1, sequence2, 2, -1, -0.5, -0.1)
    for a in alignments:
        print(format_alignment(*a))

pairwise2_global("ACT", "ACGT")

def pairwise2_global_file(sequence_file1, sequence_file2):
    seq1 = SeqIO.read(sequence_file1, "fasta")
    seq2 = SeqIO.read(sequence_file2, "fasta")
    alignments = pairwise2.align.globalds(seq1.seq, seq2.seq, blosum62, -10, -0.5)
    for a in alignments:
        print(format_alignment(*a))

pairwise2_global_file("alpha.faa", "beta.faa")

def pairwise2_local(sequence1, sequence2):
    #X = sequence1
    #Y = sequence2
    alignments = pairwise2.align.localms(sequence1, sequence2, 2, -1, -0.5, -0.1)
    for a in alignments:
        print(format_alignment(*a))

pairwise2_global("ACT", "ACGT")

def pairwise2_local_file(sequence_file1, sequence_file2):
    seq1 = SeqIO.read(sequence_file1, "fasta")
    seq2 = SeqIO.read(sequence_file2, "fasta")
    alignments = pairwise2.align.localds(seq1.seq, seq2.seq, blosum62, -10, -0.5)
    for a in alignments:
        print(format_alignment(*a))

pairwise2_global_file("alpha.faa", "beta.faa")

#needleall - this doesnt work either
"""
print("pairwise where many sequences are used, with needle")
many_needle_cline = NeedleallCommandline(asequence = "alpha.faa",bsequence = "gamma.faa",gapopen=10, gapextend=0.5, outfile="many_needle.txt")
print(many_needle_cline)

stdout, stderr = many_needle_cline()
print(stdout + stderr)

many_needle_align = AlignIO.read("many_needle.txt", "emboss")
print(many_needle_align)
"""
#https://www.biostars.org/p/237209/
#http://emboss.sourceforge.net/apps/release/6.2/emboss/apps/needleall.html
#doesnt do many-many pairwise, only one pairwise alignment
#http://emboss.sourceforge.net/apps/release/6.6/emboss/apps/needleall.html





