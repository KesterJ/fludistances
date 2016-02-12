

#IMPORTS

import Bio
from Bio.Align import Applications
from Bio import SeqIO
from Bio import AlignIO
from io import StringIO
import copy
import csv
import sys
import pandas as pd

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
	sites). Returns a list containing only those sequences which fulfil
	this criterion.
	"""
	onlyACTGlist = []
	notACTGlist = []
	for seq in seqlist:
		if all(i in DNA for i in seq):
			onlyACTGlist.append(seq)
		else:
			notACTGlist.append(seq)
	seqlist_to_file(notACTGlist, 'notACTG%s'%args.outfile)
	return onlyACTGlist


#Rewrite this for pandas, and to just return a number of differing positions instead of a list
def count_differing_pos(seqlist):
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
	return templist


def hamming_distance(s1, s2):
    """Return the Hamming distance (the number of positions at which they differ) between equal-length
    strings or sequence objects. For simplicity, sites which have a gap in either sequence are currently
    ignored and do not contribute to the distance."""
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    if args.ignoregaps:
    	return sum(bool(ord(ch1) - ord(ch2)) for ch1, ch2 in zip(s1, s2) if ch1 != '-' and ch2 != '-')
    else:
    	return sum(bool(ord(ch1) - ord(ch2)) for ch1, ch2 in zip(s1, s2))

#Older data structure since replaced with data frame - code retained for now
def get_distancedict(seqlist):
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
			distdict[key] = hamming_distance(seqlist[i].seq, seqlist[j].seq)
	return distdict


def get_distanceframe(seqlist):
	"""
	Takes a list of sequences, and returns a data frame of Hamming distances between each pair - only the
	upper triangular matrix is filled, and the rest return NaN.
	"""
	descriptions = [r.description for r in seqlist]
	distframe = pd.DataFrame(index=descriptions, columns=descriptions)
	for i in range(len(seqlist)):
		print('For sequence %s...' %i)
		for j in range(i+1, len(seqlist)):
			distframe.loc[descriptions[i], descriptions[j]] = hamming_distance(seqlist[i].seq, seqlist[j].seq)
	return distframe



def drop_multi_indices(df):
	"""
	Drops all indices which are repeated after removing the GenBank ID (assume this is due to two
	separate groups sequencing the same isolate). This isn't an ideal solution, but it'll do for now.
	NB: This should ideally check if sequences are equal for duplicates - if they are, keep one; if not, drop both (as it's hard
	to know why they're different and so they're unreliable). 
	"""
	#TODO: KJ - add sth to keep one of each pair duplicates if their seqs are identical.
	duplicates = set([x for x in df.index.tolist() if df.index.tolist().count(x) > 1])
	df = df.drop(duplicates)
	df = df.drop(duplicates, 1)
	#Outputs the list of dropped isolates - in case useful.
	if duplicates:
		list_to_file(duplicates, 'duplicates-%s'%args.outfile)
	return df



def prep_frame(df):
	"""
	Cuts out anything after the ID and description from keys of data frame index and columns. This
	makes it much easier to match indices for frames containing info about different proteins later.
	"""
	newindex = []
	for i in df.index:
		newindex.append(i.split('|')[1])
	df.index = newindex
	df.columns = df.index
	df = drop_multi_indices(df)
	return df


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


def file_to_distances(in_file, format):
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
	distframe = get_distanceframe(cleanlist)
	return distframe


#Older data structure since replaced with data frame - code retained for now
def dict_to_csv(dict1, filename):
	"""
	Quick wrapper for writing dictionary to a csv file, using csv module. 
	"""
	with open(filename, 'w') as in_hndl:
		writer = csv.DictWriter(in_hndl, dict1.keys())
		writer.writeheader()
		writer.writerow(dict1)


def frame_to_csv(frame1, filename):
	with open(filename, 'w') as in_hndl:
		frame1.to_csv(in_hndl)


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
	global args
	args = parse_clargs(clargs)
	frame1 = file_to_distances(args.infile, args.format)
	prepframe1 = prep_frame(frame1)
	nodupframe1 = drop_multi_indices(prepframe1)
	frame_to_csv(nodupframe1, args.outfile)


if __name__ == '__main__':
	try:
	    main (sys.argv[1:])
	except Exception as e:
	    print (e)
	    raise