from Bio import SwissProt
from Bio import SeqIO
from Bio import AlignIO
from Bio import Phylo
import io
import subprocess
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align.Applications import MuscleCommandline
import gzip

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
"""
print("MUSCLE Alignment:")
print("Muscle with a small input:")
muscle_cline = MuscleCommandline(input="opuntia.fasta")
stdout, stderr = muscle_cline()
muscle_align = AlignIO.read(io.StringIO(stdout), "fasta")
print(muscle_align)

"""The approach used above is simple but if you are dealing with very large output text the fact that all of stdout and
stderr is loaded into memory as a string can be a potential drawback.
Using the 'subprocess' module we can work directly with handles instead: """

print("Muscle with a large input:")

muscle_cline = MuscleCommandline(input="opuntia.fasta")
child = subprocess.Popen(str(muscle_cline),
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE,
                         universal_newlines=True,
                         shell=(sys.platform!="win32"))
muscle_align = AlignIO.read(child.stdout, "fasta")
print(muscle_align)
#name sys is not defined. ask jen


