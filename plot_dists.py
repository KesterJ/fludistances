
import sys
import csv
import matplotlib.pyplot as plt
import random


def get_subdict(indict, phrase):
	subdict = {}
	for key in indict:
		if phrase in key:
			subdict[key] = indict[key]
	return subdict



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



def csv_to_dict(filename):
	"""
	Simple wrapper for csv DictReader functionality.
	"""
	with open(filename, 'r') as in_hndl:
		indict = [i for i in csv.DictReader(in_hndl)]
	return indict[0]



def random_sample(list1, list2, factor):
	random.seed()
	length = min(len(list1), len(list2))
	samplelist = random.sample(range(0,length), int(length/factor))
	return samplelist



def get_subset(datalist, randlist):
	#This needs to be rewritten to take a subset from the original dicts, not the lists generated
	#from them ready for plotting.
	sublist = []
	for rand in randlist:
		sublist.append(datalist[rand])
	return sublist



def plot_dicts(dict1, dict2):
	print('Reformatting dictionaries...')
	plotlists = prep_for_plot(dict1, dict2)
	list1 = plotlists[0]
	list2 = plotlists[1]
	#print('Taking subset...')
	#randsamples = random_sample(list1, list2, 100)
	#sublist1 = get_subset(list1, randsamples)
	#sublist2 = get_subset(list2, randsamples)
	print('Plotting...')
	plt.scatter(list1, list2, marker='+')
	plt.show()
	#Add sth here to save plot to outfile if one if provided in command line


def get_multi_plots(dict1, dict2, phraselist):
	for i in range(len(phraselist)-1):
		for j in range(i+1, len(phraselist)):
			subdict1 = get_subdict(dict1, phraselist[i])
			subdict2 = get_subdict(dict2, phraselist[j])
			plot_dicts(subdict1, subdict2)
			savepath = 'HA-NA-%s-%s.png' %(phraselist[i], phraselist[j])
			plt.savefig(savepath)
			plt.clf()


def parse_clargs (clargs):
	import argparse
	aparser = argparse.ArgumentParser()

	aparser.add_argument ("infile1",
		help='Dictionary of distances for gene 1')

	aparser.add_argument ("infile2",
		help='Dictionary of distances for gene 2')

	aparser.add_argument("-o" "--outfile",
		help='Name of the file you want to output to')

	args = aparser.parse_args()

	return args	


def main(clargs):
	args = parse_clargs(clargs)
	dict1 = csv_to_dict(args.infile1)
	dict2 = csv_to_dict(args.infile2)
	phraselist = ['2010', '2011', '2012', '2013', '2014', '2015']
	get_multi_plots(dict1, dict2, phraselist)


if __name__ == '__main__':
	try:
	    main (sys.argv[1:])
	except Exception as e:
	    print (e)
	    raise