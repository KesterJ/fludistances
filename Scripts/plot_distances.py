"""
This takes csv files of gene distances (produced by distances.py) and plots them against each other - usually with
a selected feature highlighted.
"""


import sys
import csv
import matplotlib.pyplot as plt
import random

#CONSTS
#TODO: These are currently taken manually from relevant datasets; ideally want a way to get them automatically.
YEARLIST = ['2010', '2011', '2012', '2013', '2014', '2015']
PLACELIST = ['Brisbane', 'Perth', 'Sydney', 'Victoria', 'Santa Cruz', 'Santiago', 
	'HuNan', 'Nanjing', 'Suzhou', 'Bogota', 'Santo Domingo', 'San Salvador', 'Pavia',
	'Hokkaido', 'Kyoto', 'Nagasaki', 'Niigata', 'Mexico', 'Netherlands', 'Managua',
	'Nicaragua', 'Peru', 'Singapore', 'Thailand', 'Alabama', 'Alaska', 'Arizona',
	'Arkansas', 'California', 'Santa Clara', 'Colorado', 'Connecticut', 'Delaware',
	'District Of Columbia', 'Florida', 'Georgia', 'Hawaii', 'Idaho', 'Chicago', 'Illinois',
	'Indiana', 'Iowa', 'Kansas', 'Kentucky', 'Louisiana', 'Maine', 'Bethesda', 'Maryland',
	'Boston', 'Massachusetts', 'Michigan', 'Minnesota', 'Mississippi', 'Missouri', 'Montana',
	'Nebraska', 'Nevada', 'New Hampshire', 'New Jersey', 'New Mexico', 'New York',
	'North Carolina', 'North Dakota', 'Ohio', 'Oklahoma', 'Oregon', 'Pennsylvania',
	'Rhode Island', 'South Carolina', 'South Dakota', 'Tennessee', 'Houston', 'Texas',
	'Utah', 'Vermont', 'Washington', 'Virginia', 'Wisconsin', 'Wyoming']

#CODE

def get_subdict(indict, phrase):
	"""
	Takes a phrase and a dictionary, and returns a subdictionary consisting of all the elements
	that contain that phrase in their key.
	"""
	subdict = {}
	for key in indict:
		if phrase in key:
			subdict[key] = indict[key]
	return subdict



def get_anti_subdict(indict, phrase):
	"""
	Takes a phrase and a dictionary, and returns a subdictionary consisting of all the elements
	that DO NOT contain that phrase in their key.
	"""
	subdict = {}
	for key in indict:
		if phrase not in key:
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




def plot_dicts(dict1, dict2, color):
	"""
	Plots the two provided dictionaries against each other, in the provided colour.
	"""
	plotlists = prep_for_plot(dict1, dict2)
	list1 = plotlists[0]
	list2 = plotlists[1]
	plt.scatter(list1, list2, marker='.', color=color)




def plot_feature(dict1, dict2, feature, savepath):
	"""
	This plots the data using the provided feature; all points where at least one of the two strains being compared
	is part of that feature (e.g a place, or year) will be displayed in green. The others will be displayed in red.
	"""
	subdict1a = get_subdict(dict1, feature)
	subdict1b = get_anti_subdict(dict1, feature)
	print("Plotting for %s" %feature)
	plot_dicts(subdict1b, dict2, 'g')
	plot_dicts(subdict1a, dict2, 'r')
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
	for place in PLACELIST:
		plot_feature(dict1, dict2, place, 'HA-NA-%s.png' %place)



if __name__ == '__main__':
	try:
	    main (sys.argv[1:])
	except Exception as e:
	    print (e)
	    raise