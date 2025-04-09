#importing libraries
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import gzip
#using wget to get the human chromosome from https://www.gencodegenes.org/human/
#function to look through human genome data to get gene position and chromosome
def load_gencode(gtf_file):
    gencode_data = []
    
    # Open the compressed GTF file
    with gzip.open(gtf_file, "rt", encoding="utf-8") as file:
        for line in file:
            if line.startswith("#"):  # Skip header lines
                continue
            columns = line.strip().split("\t")
            if columns[2] == "gene":  # Only process gene entries
                chrom = columns[0]  # Chromosome
                start = int(columns[3])  # Gene start position
                
                # Extract attributes (gene_name)
                attributes = {kv.split()[0]: kv.split()[1].strip('"') for kv in columns[8].split(";") if kv}
                if "gene_name" in attributes:
                    gencode_data.append([attributes["gene_name"], chrom, start])
    
    return pd.DataFrame(gencode_data, columns=["gene_name", "chromosome", "position"])

#importing results from twas
results = pd.read_csv('MetaXcan/software/TWAS_result.csv')
#creating -log p column
results['-logp'] = np.log(results['pvalue'])
#loads in gene data
gencode_df = load_gencode("gencode.v47.basic.annotation.gtf.gz")  # Update with actual file name
#merges the two bases
merged_df = results.merge(gencode_df, on="gene_name", how="left")
#renames some of the chromosome positions from x and y to their numbered form
merged_df['chromosome'] = merged_df['chromosome'].replace({'chr1': 1, 'chr2': 2, 'chr3': 3, 'chr4': 4, 'chr5': 5, 'chr6': 6, 'chr7': 7, 'chr8': 8, 'chr9': 9, 'chr10': 10, 'chr11': 11, 'chr12': 12, 'chr13': 13, 'chr14': 14, 'chr15': 15, 'chr16': 16, 'chr17': 17, 'chr18': 18, 'chr19': 19, 'chr20': 20, 'chr21': 21, 'chr22': 22, 'chrX': 23, 'chrY': 24})
#sorts based off chromosome and position
merged_df = merged_df.sort_values(['chromosome', 'position'])
# dropping na
merged_df.dropna(subset= ['chromosome'], inplace = True)
# Compute cumulative position for Manhattan plot
merged_df['cumulative_pos'] = 0 
running_pos = 0
cumulative_pos ={}

for chrom in sorted(merged_df['chromosome'].unique()):
    chrom_subset = merged_df[merged_df['chromosome'] == chrom]
    cumulative_pos[chrom] = running_pos
    merged_df.loc[merged_df['chromosome'] == chrom, 'cumulative_pos'] = chrom_subset['position'] + running_pos
    running_pos += chrom_subset['position'].max() 

merged_df['SNP_number'] = merged_df.index
#for bug fixing
merged_df.to_csv('mergeddata.csv', index = False)
# Set the color palette
colors = sns.color_palette("husl", n_colors=len(merged_df['chromosome'].unique()))

#define figure size
sns.set(rc={"figure.figsize":(10, 5)}) #width=8, height=4

# Create the plot
g = sns.scatterplot(
    x=merged_df['cumulative_pos'], 
    y=-(merged_df['-logp']), 
    hue='chromosome', 
    palette='Set1', 
    data=merged_df, 
    legend=None
)

# Add labels and title
g.set_xlabel('chromosome')
g.set_ylabel('-log10(p-value)')
g.set_title('Manhattan Plot of TWAS Results')

# Compute median positions for chromosome labels
xticks = merged_df.groupby('chromosome')['cumulative_pos'].median()
#print(xticks)
# Add chromosome labels
g.set_xticks(xticks)
g.set_xticklabels(xticks.index.astype(int)) 
# Add a significance threshold line
g.axhline(y=-np.log10(5e-8), color='red', linestyle='--')
#save as png
fig = g.get_figure()
fig.savefig("manhattan_plot.png", dpi=300)
# Show the plot