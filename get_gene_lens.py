# get_gene_lens

import pandas as pd

gtf = '/data/users/freese/mortazavi_lab/ref/gencode.v29/gencode.v29.annotation.gtf'

df = pd.read_csv(gtf, sep='\t', usecols=[2,3,4,8],names=['entry_type', 'start', 'stop', 'fields'],comment='#')

df = df.loc[df.entry_type == 'gene']
df['gene_len'] = df.start-df.stop
df['gene_len'] = df.gene_len.abs()

df['gid'] = df.fields.str.extract(r'gene_id "([A-z]+[0-9.]+)";',expand=False)
# df['tid'] = df.fields.str.extract(r'transcript_id "([A-z]+[0-9.]+)";',expand=False)
# df['gname'] = df.fields.str.extract(r'gene_name "([A-z]+[0-9.]+)";',expand=False)



df = df[['gid', 'gene_len']]

df.to_csv('gencode_v29_gene_lens.tsv', sep='\t')