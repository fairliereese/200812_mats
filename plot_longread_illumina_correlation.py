import pandas as pd
import argparse
import scipy.stats as stats
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

def get_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('-f', dest='ab_file',
		help='Unfiltered TALON abundance file')
	parser.add_argument('-ik1', dest='ik1', 
		help='Abundance file from Kallisto, rep1')
	parser.add_argument('-ik2', dest='ik2', 
		help='Abundance file from Kallisto, rep2')
	parser.add_argument('-dname', dest='dname',
		help='Dataset name in TALON abundance file')
	parser.add_argument('-celltype', dest='celltype',
		help='Celltype of data')
	parser.add_argument('-gt_type', dest='gt_type',
		help='Correlation type to run. Choose gene or transcript.')
	parser.add_argument('-o', dest='oprefix', 
		help='Output file prefix')
	args = parser.parse_args()
	return args

def process_longread_data(ab_file, count_col, gt_type):

	df = pd.read_csv(ab_file, sep='\t')
	usecols = ['annot_gene_id', 'annot_transcript_id',
		'gene_novelty', 'transcript_novelty', count_col]
	df = df[usecols]

	# only get entries expressed above 0 counts
	df = df.loc[df[count_col] > 0]

	# transcripts: only get known transcripts
	if gt_type == 'transcript':
		id_col = 'annot_transcript_id'
		df = df.loc[df.transcript_novelty == 'Known']

	# genes: only get known genes and aggregate on genes
	elif gt_type == 'gene':
		id_col = 'annot_gene_id'
		df = df.loc[df.gene_novelty == 'Known']
		df = df.loc[df.transcript_novelty != 'Genomic']
		df = df[['annot_gene_id', count_col]].groupby(['annot_gene_id']).sum()
		df.reset_index(inplace=True)

	# compute tpm
	total_counts = df[count_col].sum()
	df['tpm'] = (df[count_col]/total_counts)*1000000
	df.rename({'annot_gene_id':'gid', 'annot_transcript_id':'tid'},
		axis=1, inplace=True)

	if gt_type == 'transcript':
		usecols = ['gid', 'tid', 'tpm']
	elif gt_type == 'gene':
		usecols = ['gid', 'tpm']

	df = df[usecols]

	return df

def process_shortread_data(ab_file, gt_type):

	df = pd.read_csv(ab_file, sep='\t')

	df = df.loc[df.tpm > 0]

	# pull the gene and transcript ids out 
	df[['tid', 'gid', 'blah']] = df.target_id.str.split(pat='|',
		expand=True, n=2)
	df.drop('blah', axis=1, inplace=True)

	# genes: aggregate on genes
	if gt_type == 'gene':
		df = df[['gid', 'tpm']].groupby(['gid']).sum()
		df.reset_index(inplace=True)

	if gt_type == 'transcript':
		usecols = ['gid', 'tid', 'tpm']
	elif gt_type == 'gene':
		usecols = ['gid', 'tpm']

	df = df[usecols]

	return df

def merge_shortread_reps(r1, r2, gt_type):
	if gt_type == 'transcript':
		df = r1.merge(r2, on=['gid', 'tid'])
	elif gt_type == 'gene':
		df = r1.merge(r2, on='gid')

	# only allow things detected in both replicates
	df = df.loc[(df.tpm_x > 0) & (df.tpm_y) > 0]

	# average out tpm
	df['tpm'] = df['tpm_x'] + df['tpm_y']
	df['tpm'] = df['tpm']/2

	if gt_type == 'transcript':
		usecols = ['gid', 'tid', 'tpm']
	elif gt_type == 'gene':
		usecols = ['gid', 'tpm']
	df = df[usecols]

	return df

def merge_long_short(pb_df, i_df, gt_type):

	print('{} {}s in PacBio'.format(len(pb_df.index), gt_type))
	print('{} {}s in Illumina'.format(len(i_df.index), gt_type))

	if gt_type == 'transcript':
		df = pb_df.merge(i_df, on=['gid','tid'], suffixes=['_long', '_short'])
	elif gt_type == 'gene':
		df = pb_df.merge(i_df, on='gid', suffixes=['_long', '_short'])

	print('{} shared {}s'.format(len(df.index), gt_type))

	return df

def plot_corr(df, gt_type, celltype, oprefix, desc=None):

	df['log_tpm_long'] = np.log10(df.tpm_long+0.1)
	df['log_tpm_short'] = np.log10(df.tpm_short+0.1)

	if gt_type == 'transcript':
		gt_type_2 = 'Transcript'
		color = '#009E73'
	elif gt_type == 'gene':
		gt_type_2 = 'Gene'
		color = 'navy'

	ax = sns.jointplot(data=df, x='log_tpm_long', y='log_tpm_short',
		edgecolors=None, s=10, color=color, alpha=0.5)
	xlabel = '{} log10(TPM+0.1) in {} PacBio'.format(gt_type_2, celltype.upper())
	ylabel = '{} log10(TPM+0.1) in Illumina'.format(gt_type_2, celltype.upper())
	ax.set_axis_labels(xlabel, ylabel)


	fname = '{}{}_{}_expression_correlation'.format(oprefix, celltype, gt_type)

	if desc:
		fname+='_{}'.format(desc)
	fname+='.png'	

	plt.savefig(fname, bbox_inches='tight')
	plt.clf()

def plot_corr_color_len(df, gt_type, celltype, oprefix):

	if gt_type == 'transcript':
		gt_type_2 = 'Transcript'
	elif gt_type == 'gene':
		gt_type_2 = 'Gene'

	len_col = '{}_len'.format(gt_type)

	df['log_tpm_long'] = np.log10(df.tpm_long+0.1)
	df['log_tpm_short'] = np.log10(df.tpm_short+0.1)
	df['log10(len)'] = np.log10(df[len_col])

	x = df.log_tpm_long.tolist()
	y = df.log_tpm_short.tolist()
	c = df['log10(len)'].tolist()

	plt.hexbin(x, y, C=c)
	plt.xlabel('{} log10(TPM+0.1) in {} PacBio'.format(gt_type_2, celltype.upper()))
	plt.ylabel('{} log10(TPM+0.1) in Illumina'.format(gt_type_2, celltype.upper()))
	cbar = plt.colorbar()
	cbar.set_label('log10({} length)'.format(gt_type))

	fname = '{}{}_{}_expression_correlation_len.png'.format(oprefix, celltype, gt_type)
	plt.savefig(fname, bbox_inches='tight')
	plt.clf()

def compute_corr(df, gt_type, celltype, oprefix, desc=None):

	short_data = df.tpm_short.tolist()
	long_data = df.tpm_long.tolist()
	pearson, _ = stats.pearsonr(short_data, long_data)
	spearman, _ = stats.spearmanr(short_data, long_data)

	s = 'Correlations for {} {}s'.format(celltype, gt_type)
	if desc:
		s+=', {}:'.format(desc)
	else:
		s+=':'

	s += '\nPearson r: {} \nSpearman rho: {}'.format(pearson, spearman)
	print(s)

	fname = '{}{}_{}_expression_correlations'.format(oprefix, celltype, gt_type)

	if desc:
		fname+='_{}'.format(desc)
	fname+='.txt'	
	
	with open(fname, 'w') as ofile:
		ofile.write('Pearson r: {}\n'.format(pearson))
		ofile.write('Spearman rho: {}'.format(spearman))

def rm_short_genes(df):
	df = df.loc[df.gene_len > 200]
	return df
	

def add_gene_lens(df):
	len_df = pd.read_csv('gencode_v29_gene_lens.tsv', sep='\t')
	df = df.merge(len_df, on='gid')
	df.rename({'len':'gene_len'}, axis=1, inplace=True)

	return df

def add_transcript_lens(df):
	len_df = pd.read_csv('gencode_v29_transcript_lens.tsv', sep='\t')
	df = df.merge(len_df, on=['gid','tid'])
	df.rename({'len':'transcript_len'}, axis=1, inplace=True)

	return df

def get_gencode_info(gt_type):
	df = pd.read_csv('/data/users/freese/mortazavi_lab/ref/gencode.v29/gencode.v29.annotation.gtf',
			sep='\t', usecols=[0,2,8], names=['chr', 'entry_type', 'fields'],
			comment='#')
	df = df.loc[df.entry_type == gt_type]

	# remove pseudogenes
	df = df.loc[df.fields.str.contains('pseudogene') == False]

	# remove mitochondrial genes
	df = df.loc[df.chr != 'chrM']

	# get gid and tid
	df['gid'] = df.fields.str.extract(r'gene_id "([A-z0-9-.]+)',
		expand=False)
	if gt_type == 'transcript':
		df['tid'] = df.fields.str.extract(r'transcript_id "([A-z0-9-.]+)',
			expand=False)

	if gt_type == 'transcript':
		usecols = ['gid', 'tid']
	elif gt_type == 'gene':
		usecols = ['gid']
	df = df[usecols]

	return df

def merge_gencode(df, g_df, gt_type):
	if gt_type == 'transcript':
		df = df.merge(g_df, on=['tid', 'gid'])
	elif gt_type == 'gene':
		df = df.merge(g_df, on='gid')
	return df

def get_single_isoform_genes(df):
	df = df.loc[~df.gid.duplicated(keep=False)]
	df.drop('tid', axis=1, inplace=True)
	return df

def main():
	args = get_args()
	ab_file = args.ab_file
	ik1 = args.ik1
	ik2 = args.ik2 
	dname = args.dname
	celltype = args.celltype
	gt_type = args.gt_type
	oprefix = args.oprefix

	pb_df = process_longread_data(ab_file, dname, gt_type)
	i1_df = process_shortread_data(ik1, gt_type)
	i2_df = process_shortread_data(ik2, gt_type)
	i_df = merge_shortread_reps(i1_df, i2_df, gt_type)

	df = merge_long_short(pb_df, i_df, gt_type)

	df = add_gene_lens(df)
	if gt_type == 'transcript':
		df = add_transcript_lens(df)

	# remove genes that are shorter than 200 bp
	df = rm_short_genes(df)

	# remove mitochondrial genes and pseudogenes
	g_df = get_gencode_info(gt_type)
	df = merge_gencode(df, g_df, gt_type)

	print('{} {}s after filtering'.format(len(df.index), gt_type))

	# plot correlation between transcripts
	plot_corr(df, gt_type, celltype, oprefix)
	compute_corr(df, gt_type, celltype, oprefix)

	# plot correlation between transcipts, color by gene or
	# transcript length
	plot_corr_color_len(df, gt_type, celltype, oprefix)

	# make another plot but restrict it only to genes that contain
	# one isoform (should this be isoforms detected or annotated?)
	if gt_type == 'gene':
		g_df = get_gencode_info('transcript')
		g_df = get_single_isoform_genes(g_df)
		df = merge_gencode(df, g_df, gt_type)

		plot_corr(df, gt_type, celltype, oprefix, 'single_isoform')
		compute_corr(df, gt_type, celltype, oprefix, 'single_isoform')
	print()

if __name__ == '__main__': main()