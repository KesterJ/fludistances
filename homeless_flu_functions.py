

#CONSTS
COLOURS = ['r', 'b', 'g', 'm', 'y', 'k']
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
	

def plot_feature(dict1, dict2, placename):
	subdict1a = get_subdict(dict1, placename)
	subdict1b = get_anti_subdict(dict1, placename)
	plot_dicts(subdict1b, dict2, 'g')
	plot_dicts(subdict1a, dict2, 'r')
	savepath = 'HA-NA-%s.png' %placename
	plt.savefig(savepath)
	plt.clf()


#This needs to be a proper rigorous function using the hypergeometric and not just a quick bodge to
#look at some rough numbers.
def find_reassortants(dict1, dict2):
	"""
	Takes two dictionaries with the same key list, compares numbers and prints those with a ratio >3.
	This finds outliers on the graph, which are the reassortants.
	"""
	#Add sth to deal with differing keys
	list1 = []
	list2 = []
	for key in dict1:
		if key in dict2:
			if int(dict1[key]) > 0 and int(dict2[key]) > 0:
				if int(dict1[key])/int(dict2[key]) > 3 or int(dict2[key])/int(dict1[key]) > 3:
					print(key)