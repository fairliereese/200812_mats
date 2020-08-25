import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


df = pd.read_csv('sirv_numbers.tsv', sep='\t')
ax = sns.barplot(x='dataset', y='percent_sirvs', data=df)
ax.set_xlabel('Sample')
ax.set_ylabel('Percentage of SIRV reads')

plt.savefig('figures/sirv_bar.png')