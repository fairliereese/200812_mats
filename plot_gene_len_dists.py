import argparse
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

def get_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('-f', dest='infile',
		help='Input shared genes file')
	parser.add_argument('-o', dest='oprefix',
		help='Output file prefix m')
	parser.add_argument('-c', dest='celltype')
	args = parser.parse_args()
	return args

def main():

	args = get_args()
	infile = args.infile
	celltype = args.celltype
	oprefix = args.oprefix

	ab_df = pd.read_csv(infile)
	len_df = pd.read_csv('gencode_v29_gene_lens.tsv', sep='\t')
	df = ab_df.merge(len_df, how='inner', left_on='gene_ID', right_on='gid')
	lens = df.gene_len.tolist()
	lens = [i for i in lens if i < 500000]
	ax = sns.distplot(lens)
	ax.set_xlim([0,500000])
	plt.title('{} Illumina/PacBio Shared Gene Lengths'.format(celltype))
	plt.savefig(oprefix+'_shared_gene_len_dist.png')
	plt.clf()

if __name__ == '__main__':
	main()