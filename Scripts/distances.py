

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
	Takes a list of alignments, and returns a list of sequence objects derived from the same sequences,
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



def check_coding(seqlist):
	"""
	Takes a list of sequence objects, and checks that the first three non '-' characters of each are
	'ATG' (can be in caps or lower case). Returns a list containing only those sequences which fulfil
	this criterion.
	"""
	codinglist = []
	notcodinglist = []
	for seq in seqlist:
		x = 0
		#While loop to set x and make sure we start from the first nucleotide if there are gaps
		#at the start of the alignment
		while seq[x] == '-' and x <= len(seq)-3:
			x += 1
		if seq[x+0] in ('a', 'A') and seq[x+1] in ('t', 'T') and seq[x+2] in ('g', 'G'):
			codinglist.append(seq)
		else:
			notcodinglist.append(seq)
	if notcodinglist:
		seqlist_to_file(notcodinglist, 'notcoding%s'%args.outfile)
	return codinglist



def check_ACTG(seqlist):
	"""
	Takes a list of sequence objects, and checks that the only characters present are 'A' 'T' 'C' 'G'
	and '-' (letters can be lower case or upper case). (i.e. rejects sequences containing ambiguous
	sites). Other chars are replaced with gaps (this isn't biologically perfect, but works for our
	purposes here as the Hamming distance function essentially ignores gaps. Returns a list containing
	only those sequences which fulfil this criterion.
	"""
	onlyACTGlist = copy.deepcopy(seqlist)
	for seq in onlyACTGlist:
		for base in seq.seq:
			if base not in DNA:
				base = '-'
				#Alternatively, can just drop the sequence out - but retaining it with reduced info seems better.
	return onlyACTGlist



def hamming_distance(s1, s2, ignoregaps):
    """Return the Hamming distance (the number of positions at which they differ) between equal-length
    strings or sequence objects. For simplicity, sites which have a gap in either sequence are currently
    ignored and do not contribute to the distance."""
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    if ignoregaps == True:
    	return sum(bool(ord(ch1) - ord(ch2)) for ch1, ch2 in zip(s1, s2) if ch1 != '-' and ch2 != '-')
    else:
    	return sum(bool(ord(ch1) - ord(ch2)) for ch1, ch2 in zip(s1, s2))


def get_distancedict(seqlist, ignoregaps):
	"""
	Takes a list of sequences, and returns a data frame of Hamming distances between each pair - only the
	upper triangular matrix is filled, and the rest return NaN.
	"""
	distdict = {}
	for i in range(len(seqlist)):
		print('For sequence %s...' %i)
		for j in range(i+1, len(seqlist)):
			#As the isolate name is stored within the description of a SeqRecord along with a lot of other data,
			#the next line is to extract specifically that part.
			#???: Split this off into a function to make code more understandable?
			key = seqlist[i].description.split('|')[1] + seqlist[j].description.split('|')[1]
			distdict[key] = hamming_distance(seqlist[i].seq, seqlist[j].seq, ignoregaps)
	return distdict



def list_to_file(outlist, savepath):
	"""
	Writes a list to a file, one item per line.
	"""
	with open(savepath, 'w') as out_hndl:
		for item in outlist:
  			out_hndl.write('%s\n' % item)


def seqlist_to_file(seqlist, savepath):
	"""
	Writes a list of sequences to a file from sequence object, one item per line.
	"""
	with open(savepath, 'w') as out_hndl:
  		for item in seqlist:
  			out_hndl.write('%s\n' %item.description)	
  			out_hndl.write('%s\n' %str(item.seq))


def dict_to_csv(dict1, filename):
	"""
	Quick wrapper for writing dictionary to a csv file, using csv module. 
	"""
	with open(filename, 'w') as in_hndl:
		writer = csv.DictWriter(in_hndl, dict1.keys())
		writer.writeheader()
		writer.writerow(dict1)


def file_to_distances(in_file, out_file, format, ignoregaps):
	"""
	This should just be a wrapper for all the other functions, that calls them in order.
	"""
	print('Reading files...')
	alignlist = align_file(in_file, format)

	print('Not coding seqs:')
	codinglist = check_coding(alignlist)
	print('Length coding list = %s' %len(codinglist))

	print('Extracting third positions')
	thirdposlist = get_third_pos(codinglist)
	print('Length thirdposlist = %s' %len(thirdposlist))

	cleanlist = check_ACTG(thirdposlist)
	print('Length cleanlist = %s' %len(cleanlist))

	print('Finding Hamming distances...')
	distdict = get_distancedict(cleanlist, ignoregaps)
	
	dict_to_csv(distdict, out_file)





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
		default='fasta')

	aparser.add_argument("-i", "--ignoregaps",
		help='Should gap nucleotides be ignored when calculating Hamming distance',
		type=bool,
		default=True)

	args = aparser.parse_args()

	return args



def main(clargs):
	args = parse_clargs(clargs)
	distdict = file_to_distances(args.infile, args.outfile, args.format, args.ignoregaps)


if __name__ == '__main__':
	try:
	    main (sys.argv[1:])
	except Exception as e:
	    print (e)
	    raise