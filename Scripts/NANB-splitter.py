from Bio import SeqIO

def get_recs(seqlist, word):
	newseqlist = []
	for seq in seqlist:
		if word in seq.description:
			newseqlist.append(seq)
	return newseqlist

def main():
	with open('Human Flu B data all time/NANB-cds-all-fluB.fasta', 'r') as in_hndl:
		seqlist = [r for r in SeqIO.parse(in_hndl, 'fasta')]
	print('Original length', len(seqlist))
	NAlist = get_recs(seqlist, 'Neuraminidase')
	NBlist = get_recs(seqlist, 'NB protein')
	print('NA length', len(NAlist))
	print('NB length', len(NBlist))
	with open('NA-cds-all-fluB.fasta', 'w') as out_hndl:
		SeqIO.write(NAlist, out_hndl, 'fasta')
	with open('NB-cds-all-fluB.fasta', 'w') as out_hndl:
		SeqIO.write(NBlist, out_hndl, 'fasta')
	for seq in seqlist:
		if 'Neuraminidase' not in seq.description and 'NB protein' not in seq.description:
			print(seq)

main()