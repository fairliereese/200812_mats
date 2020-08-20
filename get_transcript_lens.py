import pandas as pd

infile = '/data/users/freese/mortazavi_lab/ref/gencode.v29/gencode.v29.transcripts.fa'

dictuel = dict()

f = open(infile, 'r')
first = True
for line in f:
	if line.startswith('>'):
		if not first:
			entry = {'gid': gid,
					 'len': curr_len}
			dictuel[tid] = entry
		line = line[1:-1]
		beep = line.split('|')
		gid = beep[1]
		tid = beep[0]
		curr_len = 0
		first = False
	else:
		line = line[:-1]
		curr_len+= len(line)

# add last transcript
entry = {'gid': gid,
		 'len': curr_len}
dictuel[tid] = entry

df = pd.DataFrame.from_dict(dictuel, orient='index')
df.to_csv('gencode_v29_transcript_lens.tsv', sep='\t')

f.close()
