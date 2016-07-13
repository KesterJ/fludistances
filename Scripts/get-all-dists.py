import distances
import os


for filename in os.listdir('../Input FASTA files from FluDB/'):
	if not filename == '.DS_Store':
		print('Working on %s'%filename)
		distances.file_to_distances('../Input FASTA files from FluDB/%s'%filename, 'Working files/%s-distances.csv'%filename, 'fasta', True)