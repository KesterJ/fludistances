
import Bio
from Bio.Align import Applications
from Bio import SeqIO
from Bio import AlignIO
from io import StringIO
import copy

#Align them
def align_genes(in_file):
	mafft_cline = Applications.MafftCommandline('/usr/local/bin/mafft',input=in_file)
	stdout, stderr = mafft_cline()
	align = AlignIO.read(StringIO(stdout), "fasta")
	return align


def get_third_pos(alignlist):
	thirdlist = []
	for align in alignlist:
		thirdpos = ''.join( [align[r] for r in range(len(align)) if r%3==2] )
		tempseq = copy.deepcopy(align)
		tempseq.seq = thirdpos
		thirdlist.append(tempseq)
	return thirdlist


def get_differing_pos(seqlist):
	diffpos = []
	for nucl in range(len(thirdlist[0])):
		templist = ''.join( [thirdlist[r][nucl] for r in range(len(thirdlist))])
		if templist.count (templist[0]) != len (templist):
			diffpos.append(str(nucl)+templist)
	return diffpos


def hamming_distance(s1, s2):
    """Return the Hamming distance between equal-length sequences"""
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(bool(ord(ch1) - ord(ch2)) for ch1, ch2 in zip(s1, s2))


def get_distances(seqlist):
	distdict = {}
	for i in range(len(seqlist)):
		for j in range(i+1, len(seqlist)):
			#As the isolate is stored within the description of a SeqRecord along with a lot of other data,
			#the next line is to extract specifically that part.
			#Split this off into a function to make code more understandable?
			key = '%s %s' %seqlist[i].description.split('|')[1] %seqlist[j].description.split('|')[1]
			distdict[key] = hamming_distance(seqlist[i].seq, seqlist[j].seq)
	return distdict


def main():
	filename = 'arbitrary-cds-data.fasta'
	#genes = readseqs(filename)
	align = align_genes (filename)
	print(align[0].id)
	print(align[0].description)
	print(align[0].name)
	print(align[0].dbxrefs)
	print(align[0].features)
	print(align[0].annotations)
	#Something here to check they look like coding seqs - just check start with ATG?
	"""thirdlist = get_third_pos(align)
	diffpos = getd_iffering_pos (thirdlist)

	for seq in diffpos:
		print (seq)
	"""
main()