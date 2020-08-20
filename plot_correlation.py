import pandas as pd
import argparse
import scipy.stats as stats
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

def get_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('-s1_type', dest='s1_type',
		help='Type of sample1 to use in labels')
	parser.add_argument('-s2_type', dest='s2_type',
		help='Type of sample2 to use in labels')
	parser.add_argument('-i1_1', dest='i1_1', 
		help='Abundance file from Kallisto sample1, rep1')
	parser.add_argument('-i1_2', dest='i1_2', 
		help='Abundance file from Kallisto sample1, rep2')
	parser.add_argument('-i2_1', dest='i2_1', 
		help='Abundance file from Kallisto sample2, rep1')
	parser.add_argument('-i2_2', dest='i2_2', 
		help='Abundance file from Kallisto sample2, rep2')
	parser.add_argument('-celltype', dest='celltype',
		help='Celltype of data')
	parser.add_argument('-gt_type', dest='gt_type',
		help='Correlation type to run. Choose gene or transcript.')
	parser.add_argument('-o', dest='oprefix', 
		help='Output file prefix')
	args = parser.parse_args()
	return args

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

def merge_samples(df1, df2, gt_type):

	print('{} {}s in Sample 1'.format(len(df1.index), gt_type))
	print('{} {}s in Sample 2'.format(len(df2.index), gt_type))

	if gt_type == 'transcript':
		df = df1.merge(df2, on=['gid','tid'], suffixes=['_1', '_2'])
	elif gt_type == 'gene':
		df = df1.merge(df2, on='gid', suffixes=['_1', '_2'])

	print('{} shared {}s'.format(len(df.index), gt_type))

	return df

def plot_corr(df, gt_type, celltype, oprefix, s1, s2, desc=None):

	df['log_tpm_1'] = np.log10(df.tpm_1+0.1)
	df['log_tpm_2'] = np.log10(df.tpm_2+0.1)

	if gt_type == 'transcript':
		gt_type_2 = 'Transcript'
		color = '#009E73'
	elif gt_type == 'gene':
		gt_type_2 = 'Gene'
		color = 'navy'

	ax = sns.jointplot(data=df, x='log_tpm_1', y='log_tpm_2',
		edgecolors=None, s=10, color=color, alpha=0.5)
	xlabel = '{} log10(TPM+0.1) in {} {}'.format(gt_type_2, celltype.upper(), s1)
	ylabel = '{} log10(TPM+0.1) in {} {}'.format(gt_type_2, celltype.upper(), s2)
	ax.set_axis_labels(xlabel, ylabel)

	fname = '{}{}_{}_expression_correlation'.format(oprefix, celltype, gt_type)

	if desc:
		fname+='_{}'.format(desc)
	fname+='.png'	

	plt.savefig(fname, bbox_inches='tight')
	plt.clf()

def plot_corr_color_len(df, gt_type, celltype, oprefix, s1, s2):

	if gt_type == 'transcript':
		gt_type_2 = 'Transcript'
	elif gt_type == 'gene':
		gt_type_2 = 'Gene'

	len_col = '{}_len'.format(gt_type)

	df['log_tpm_1'] = np.log10(df.tpm_1+0.1)
	df['log_tpm_2'] = np.log10(df.tpm_2+0.1)
	df['log10(len)'] = np.log10(df[len_col])

	x = df.log_tpm_1.tolist()
	y = df.log_tpm_2.tolist()
	c = df['log10(len)'].tolist()

	plt.hexbin(x, y, C=c)
	plt.xlabel('{} log10(TPM+0.1) in {} {}'.format(gt_type_2, celltype.upper(), s1))
	plt.ylabel('{} log10(TPM+0.1) in {} {}'.format(gt_type_2, celltype.upper(), s2))
	cbar = plt.colorbar()
	cbar.set_label('log10({} length)'.format(gt_type))

	fname = '{}{}_{}_expression_correlation_len.png'.format(oprefix, celltype, gt_type)
	plt.savefig(fname, bbox_inches='tight')
	plt.clf()

def compute_corr(df, gt_type, celltype, oprefix, desc=None):

	short_data = df.tpm_2.tolist()
	long_data = df.tpm_1.tolist()
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
	p_1 = args.i1_1
	p_2 = args.i1_2
	r_1 = args.i2_1
	r_2 = args.i2_2
	s1 = args.s1_type
	s2 = args.s2_type 
	celltype = args.celltype
	gt_type = args.gt_type
	oprefix = args.oprefix

	p1_df = process_shortread_data(p_1, gt_type)
	p2_df = process_shortread_data(p_2, gt_type)
	p_df = merge_shortread_reps(p1_df, p2_df, gt_type)

	r1_df = process_shortread_data(r_1, gt_type)
	r2_df = process_shortread_data(r_2, gt_type)
	r_df = merge_shortread_reps(r1_df, r2_df, gt_type)

	df = merge_samples(p_df, r_df, gt_type)

	df = add_gene_lens(df)
	if gt_type == 'transcript':
		df = add_transcript_lens(df)

	# remove genes that are shorter than 200 bp
	df = rm_short_genes(df)

	# remove mitochondrial genes and pseudogenes
	g_df = get_gencode_info(gt_type)
	df = merge_gencode(df, g_df, gt_type)

	print('{} shared {}s after filtering'.format(len(df.index), gt_type))

	# plot correlation between transcripts
	plot_corr(df, gt_type, celltype, oprefix, s1, s2)
	compute_corr(df, gt_type, celltype, oprefix)

	# plot correlation between transcipts, color by gene or
	# transcript length
	plot_corr_color_len(df, gt_type, celltype, oprefix, s1, s2)

	# make another plot but restrict it only to genes that contain
	# one isoform (should this be isoforms detected or annotated?)
	if gt_type == 'gene':
		g_df = get_gencode_info('transcript')
		g_df = get_single_isoform_genes(g_df)
		df = merge_gencode(df, g_df, gt_type)

		plot_corr(df, gt_type, celltype, oprefix,
			s1, s2, desc='single_isoform')
		compute_corr(df, gt_type, celltype, oprefix,desc='single_isoform')
	print()

if __name__ == '__main__': main()