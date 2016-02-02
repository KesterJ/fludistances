
import Bio
from Bio.Align import Applications
from Bio import SeqIO
from Bio import AlignIO
from io import StringIO
import copy
import matplotlib.pyplot as plt
import csv


#CONSTANTS
HAFILE = 'Human H3N2 data 2010-2015/HA-cds-data.fasta'
NAFILE = 'Human H3N2 data 2010-2015/NA-cds-data.fasta'

#Align them
def align_genes(in_file):
	"""
	Takes a file name (in_file) and extracts the sequences in that file into a Biopython Align object.
	"""
	mafft_cline = Applications.MafftCommandline('/usr/local/bin/mafft',input=in_file)
	stdout, stderr = mafft_cline()
	align = AlignIO.read(StringIO(stdout), "fasta")
	return align


def get_third_pos(alignlist):
	"""
	Takes a list of alignments, and returns a list of sequence objects dervied from the same sequences,
	but with all positions apart from third codons removed. Other info remains the same.
	This doesn't check for reading frames etc - just assumes that the start of the alignment is the start
	of the reading frame and takes every third one from there.
	"""
	thirdlist = []
	for align in alignlist:
		thirdpos = ''.join( [align[r] for r in range(len(align)) if r%3==2] )
		tempseq = copy.deepcopy(align)
		tempseq.seq = thirdpos
		thirdlist.append(tempseq)
	return thirdlist

#Prob don't need this unless want to calculate e.g. some stats on which positions are most variable
def get_differing_pos(seqlist):
	"""
	Takes a list of sequences, and returns a list of strings containing only those positions that vary
	between the sequences, and appends the position in the seq to the front of each.

	Example:
	seqlist: [SeqRecord('actg'), SeqRecord('acgg'), SeqRecord('attg')]
	returns: ['1cct', '2tgt']
	"""
	diffpos = []
	for nucl in range(len(seqlist[0])):
		templist = ''.join( [seqlist[r][nucl] for r in range(len(seqlist))])
		if templist.count (templist[0]) != len (templist):
			diffpos.append(str(nucl)+templist)
	return diffpos


def hamming_distance(s1, s2):
    """Return the Hamming distance between equal-length sequences"""
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(bool(ord(ch1) - ord(ch2)) for ch1, ch2 in zip(s1, s2))


def get_distances(seqlist):
	"""
	Takes a list of sequences, and returns a dictionary of Hamming distances between each pair - the
	dictionary keys are the concatenation of the names of the two isolates.
	"""
	distdict = {}
	for i in range(len(seqlist)):
		print('For sequence %s...' %i)
		for j in range(i+1, len(seqlist)):
			#As the isolate name is stored within the description of a SeqRecord along with a lot of other data,
			#the next line is to extract specifically that part.
			#Split this off into a function to make code more understandable?
			key = seqlist[i].description.split('|')[1] + seqlist[j].description.split('|')[1]
			distdict[key] = hamming_distance(seqlist[i].seq, seqlist[j].seq)
	return distdict


def file_to_distances(in_file):
	"""
	This should just be a wrapper for all the other functions, that calls them in order.
	"""
	print('Reading files...')
	alignlist = align_genes(in_file)
	#Something here to check they look like coding seqs - just check start with ATG?
	print('Extracting third positions')
	thirdposlist = get_third_pos(alignlist)
	#The line below (calling get_differing_pos) isn't strictly necessary -
	#could just compare whole sequence, which may be more appropriate. But effect should just scale
	#graph differently - not affect final result.
	#diffposlist = get_differing_pos(thirdposlist)
	print('Finding Hamming distances...')
	distdict = get_distances(thirdposlist)
	return distdict


def prep_for_plot(dict1, dict2):
	"""
	Takes two dictionaries with the same key list, and turns them into two lists ordered in the same
	way. The point is to make two lists where the points are correctly paired for matplotlib's scatter
	function.
	"""
	#Add sth to deal with differing keys
	list1 = []
	list2 = []
	for key in dict1:
		if key in dict2:
			list1.append(dict1[key])
			list2.append(dict2[key])
	return (list1, list2)


def dict_to_csv(dict1, filename):
	"""
	Quick wrapper for writing dictionary to a csv file, using csv module. 
	"""
	with open(filename, 'w') as in_hndl:
		writer = csv.DictWriter(in_hndl, dict1.keys())
		writer.writeheader()
		writer.writerow(dict1)


def main():
	file1 = HAFILE
	file2 = NAFILE
	dict1 = file_to_distances(file1)
	dict2 = file_to_distances(file2)
	dict_to_csv(dict1, 'HA-dict.csv')
	dict_to_csv(dict2, 'NA-dict.csv')
	print('Reformatting dictionaries...')
	plotlists = prep_for_plot(dict1, dict2)
	print('Plotting...')
	plt.scatter(plotlists[0], plotlists[1], marker='+')
	plt.show()
	

if __name__ == '__main__':
	try:
	    main (sys.argv[1:])
	except Exception as e:
	    print (e)
	    raise