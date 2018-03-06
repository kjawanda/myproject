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
"""
print("Clustal alignment (includes alignment and phylogenetic tree)")
cline = ClustalwCommandline("clustalw", infile="opuntia.fasta")
print(cline)
stdout, stderr = cline()
align = AlignIO.read("opuntia.aln","clustal")
print(align)
tree = Phylo.read("opuntia.dnd", "newick")
Phylo.draw_ascii(tree)

print("MUSCLE Alignment:")
print("Muscle with a small input:")
muscle_cline = MuscleCommandline(input="opuntia.fasta")
stdout, stderr = muscle_cline()
muscle_align = AlignIO.read(io.StringIO(stdout), "fasta")
print(muscle_align)

The approach used above is simple but if you are dealing with very large output text the fact that all of stdout and
stderr is loaded into memory as a string can be a potential drawback.
Using the 'subprocess' module we can work directly with handles instead: 

print("Muscle with a large input:")

muscle_cline = MuscleCommandline(input="opuntia.fasta")
child = subprocess.Popen(str(muscle_cline),
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE,
                         universal_newlines=True,
                         shell=(sys.platform!="linux2"))
muscle_align = AlignIO.read(child.stdout, "fasta")
print(muscle_align)

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


print("pairwise global alignment:")
X = "ACGGGT"
Y = "ACG"
alignments = pairwise2.align.globalms(X, Y, 2, -1, -0.5, -0.1)
for a in alignments:
    print(format_alignment(*a))


print("pairwise local alignment:")
X = "ACGGGT"
Y = "ACG"
alignments = pairwise2.align.localms(X, Y, 2, -1, -0.5, -0.1)
for a in alignments:
    print(format_alignment(*a))

print("pairwise alignments using EMBOSS")
print("pairwise with needle")
needle_cline = NeedleCommandline(asequence = "alpha.faa", bsequence = "beta.faa",
                                 gapopen=10, gapextend=0.5, outfile="needle.txt")
print(needle_cline)
stdout, stderr = needle_cline()
print(stdout + stderr)

needle_align = AlignIO.read("needle.txt", "emboss")
print(needle_align)


print("pairwise with water")      
water_cline = WaterCommandline(asequence = "alpha.faa", bsequence = "beta.faa",
                                 gapopen=10, gapextend=0.5, outfile="water.txt")
print(water_cline)
stdout, stderr = water_cline()
print(stdout + stderr)

water_align = AlignIO.read("water.txt", "emboss")
print(water_align)

print("this is the name of the script:", sys.argv[1])


"""
print("pairwise where many sequences are used, with needle")
many_needle_cline = NeedleallCommandline(asequence = "alpha.faa",bsequence = "gamma.faa",gapopen=10, gapextend=0.5, outfile="many_needle.txt")
print(many_needle_cline)

stdout, stderr = many_needle_cline()
print(stdout + stderr)

many_needle_align = AlignIO.read("many_needle.txt", "emboss")
print(many_needle_align)

#https://www.biostars.org/p/237209/
#http://emboss.sourceforge.net/apps/release/6.2/emboss/apps/needleall.html
#doesnt do many-many pairwise, only one pairwise alignment
#http://emboss.sourceforge.net/apps/release/6.6/emboss/apps/needleall.html

"""



