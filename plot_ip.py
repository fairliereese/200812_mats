import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


df = pd.read_csv('ip_numbers.tsv', sep='\t')
ax = sns.barplot(x='dataset', y='percent_ip', data=df)
ax.set_xlabel('Sample')
ax.set_ylabel('Percentage of internally-primed reads')

plt.savefig('figures/ip_bar.png')
