
import sys
import matplotlib.pyplot as plt
import random
import pandas as pd



def get_subindex(frame1, phrase, includephrase):
	"""
	Takes a square data frame where index == columns, and a phrase, and finds all keys in the index
	that contain that phrase - if includephrase is True - or those that do not include that phrase -
	if includephrase is False. Used to find e.g. all data from a certain year or place.
	NB: This specifically works for the data from FluDB, because all the info is stored in a single
	string, which is used as the key here.
	"""
	indexlist = []
	if includephrase:
		for key in frame1.index:
			if phrase in key:
				indexlist.append(key)
	else:
		for key in frame1.index:
			if phrase not in key:
				indexlist.append(key)
	return indexlist



def index_overlap(frame1, frame2):
	"""
	Takes two data frames, each of which is square and where index == columns, and returns a list of the
	common items in the indices of the two. This will NOT necessarily be sorted in the same order as the
	input indices.
	"""
	overlap = frame1.index.intersection(frame2.index)
	overlap = overlap.tolist()
	return overlap



def get_subframe(frame1, indexlist):
	"""
	Takes a dataframe where indices and columns are the same, and a pre-prepared and filtered index list,
	which must be a subset of the indices of that dataframe. Will return a new frame containing only
	those values where at least one of its index and column are in the list.
	"""
	import time
	newframe = pd.DataFrame(index = indexlist, columns = indexlist)
	newframe.update(frame1)
	return newframe



def plot_frames(frame1, frame2, color):
	"""
	Plots two lists of data (first element of first list against first element of second list,
	second element against second element, etc)
	"""
	overlap = index_overlap(frame1, frame2)
	f1 = get_subframe(frame1, overlap)
	f2 = get_subframe(frame2, overlap)
	print('Plotting...')
	plt.scatter(f1, f2, marker='.', color=color)



def parse_clargs (clargs):
	import argparse
	aparser = argparse.ArgumentParser()

	aparser.add_argument ("infile1",
		help='Distances for gene 1')

	aparser.add_argument ("infile2",
		help='Distances for gene 2')

	aparser.add_argument("outfile",
		help='Name of the file you want to output to')

	aparser.add_argument("-f", "--filterword",
		help='Phrase to plot in alternate colour',
		type=str,
		default=None)

	args = aparser.parse_args()

	return args	


def main(clargs):
	args = parse_clargs(clargs)
	frame1 = pd.read_csv(args.infile1, index_col = 0)
	frame2 = pd.read_csv(args.infile2, index_col = 0)
	if args.filterword:
		#Only need to find subframes for one of our two, as the plot_frames function already finds the
		#overlap betwen the two frames being plotted - this overlap will naturally include the subframe
		#of the other too.
		subindex = get_subindex(frame1, args.filterword, True)
		antisubindex = get_subindex(frame1, args.filterword, False)
		subframe1 = get_subframe(frame1, subindex)
		antisubframe1 = get_subframe(frame1, antisubindex)
		plot_frames(antisubframe1, frame2, 'b')
		plot_frames(subframe1, frame2, 'r')
		plt.savefig(args.outfile, format = 'png')
	else:	
		plot_frames(frame1, frame2, 'b')
		plt.savefig(args.outfile, format = 'png')


if __name__ == '__main__':
	try:
	    main (sys.argv[1:])
	except Exception as e:
	    print (e)
	    raise