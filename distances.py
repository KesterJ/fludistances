

#IMPORTS

import Bio
from Bio.Align import Applications
from Bio import SeqIO
from Bio import AlignIO
from io import StringIO
import copy
import csv
import sys


#CONSTS
DNA = 'ACTGactg-'


#CODE

def align_file(in_file, format):
	"""
	Takes a file name (in_file) and extracts the sequences in that file into a Biopython Align object.
	"""
	mafft_cline = Applications.MafftCommandline('/usr/local/bin/mafft',input=in_file)
	stdout, stderr = mafft_cline()
	align = AlignIO.read(StringIO(stdout), format)
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


def check_ACTG(seqlist):
	onlyACTGlist = []
	for seq in seqlist:
		if all(i in DNA for i in seq):
			onlyACTGlist.append(seq)
		else:
			print(seq)
	return onlyACTGlist


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
    """Return the Hamming distance (the number of positions at which they differ) between equal-length
    strings or sequence objects"""
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
			#???: Split this off into a function to make code more understandable?
			key = seqlist[i].description.split('|')[1] + seqlist[j].description.split('|')[1]
			distdict[key] = hamming_distance(seqlist[i].seq, seqlist[j].seq)
	return distdict


def file_to_distances(in_file, format):
	"""
	This should just be a wrapper for all the other functions, that calls them in order.
	"""
	print('Reading files...')
	alignlist = align_file(in_file, format)
	#TODO: KJ - Add something here to check they look like coding seqs - just check start with ATG?
	print('Extracting third positions')
	thirdposlist = get_third_pos(alignlist)
	print('Length thirdposlist = %s' %len(thirdposlist))
	cleanlist = check_ACTG(thirdposlist)
	print('Length cleanlist = %s' %len(cleanlist))
	print('Finding Hamming distances...')
	distdict = get_distances(cleanlist)
	return distdict



def dict_to_csv(dict1, filename):
	"""
	Quick wrapper for writing dictionary to a csv file, using csv module. 
	"""
	with open(filename, 'w') as in_hndl:
		writer = csv.DictWriter(in_hndl, dict1.keys())
		writer.writeheader()
		writer.writerow(dict1)



def parse_clargs (clargs):
	import argparse
	aparser = argparse.ArgumentParser()

	aparser.add_argument ("infile",
		help='File containing full gene sequences')

	aparser.add_argument("outfile",
		help='Name of the file you want to output to')

	aparser.add_argument("-f", "--format",
		help='Format of the input sequence file',
		type=str,
		default='fasta'
	)

	args = aparser.parse_args()

	return args



def main(clargs):
	args = parse_clargs(clargs)
	dict1 = file_to_distances(args.infile, args.format)
	dict_to_csv(dict1, args.outfile)


if __name__ == '__main__':
	try:
	    main (sys.argv[1:])
	except Exception as e:
	    print (e)
	    raise