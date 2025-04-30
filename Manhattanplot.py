import argparse
import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import gzip
import os

def check_arg(args=None):
    parser = argparse.ArgumentParser(description="Multi-file Manhattan Plot Script.")
    parser.add_argument("-i", "--input", nargs='+', help="Input CSV files", required=True)
    parser.add_argument("-o", "--output", help="Output PNG filename", required=True)
    return parser.parse_args(args)

arguments = check_arg(sys.argv[1:])
input_files = arguments.input
output_file = arguments.output

def load_gencode(gtf_file):
    gencode_data = []
    
    #opens the compressed file
    with gzip.open(gtf_file, "rt", encoding="utf-8") as file:
        for line in file:
            if line.startswith("#"):  #skips header lines
                continue
            columns = line.strip().split("\t")
            if columns[2] == "gene":  #only processes gene entries
                chrom = columns[0]  #chromosome
                start = int(columns[3])  #gene start position
                
                attributes = {kv.split()[0]: kv.split()[1].strip('"') for kv in columns[8].split(";") if kv}
                if "gene_name" in attributes:
                    gencode_data.append([attributes["gene_name"], chrom, start])
    
    return pd.DataFrame(gencode_data, columns=["gene_name", "chromosome", "position"])

#loads and combines data from multiple input files
all_results = []
for file in input_files:
    results = pd.read_csv(file)
    results['-logp'] = -np.log10(results['pvalue'])  # Create -log10(p-value) column
    
    # Extract file name and cell type
    name = os.path.basename(file) 
    cell = name.split(".")[0]
    results['cell_type'] = cell.removesuffix("_output")  # Add a column for cell type
    results['file_type'] = os.path.basename(file)  # Add a column for file type
    
    all_results.append(results)
    
#combines all results into a single DataFrame
combined_results = pd.concat(all_results, ignore_index=True)

#loads the gene data
gencode_df = load_gencode("gencode.v47.basic.annotation.gtf.gz")  # Update with actual file name

#merges the combined results with gene data
merged_df = combined_results.merge(gencode_df, on="gene_name", how="left")

#renames chromosome positions from X/Y to numeric form
merged_df['chromosome'] = merged_df['chromosome'].replace({
    'chr1': 1, 'chr2': 2, 'chr3': 3, 'chr4': 4, 'chr5': 5, 'chr6': 6, 'chr7': 7, 'chr8': 8, 
    'chr9': 9, 'chr10': 10, 'chr11': 11, 'chr12': 12, 'chr13': 13, 'chr14': 14, 'chr15': 15, 
    'chr16': 16, 'chr17': 17, 'chr18': 18, 'chr19': 19, 'chr20': 20, 'chr21': 21, 'chr22': 22, 
    'chrX': 23, 'chrY': 24
})

#sort by chromosome and position
merged_df = merged_df.sort_values(['chromosome', 'position'])

#drop rows with missing chromosome data
merged_df.dropna(subset=['chromosome'], inplace=True)

#finds cumulative position for Manhattan plot with equal spacing between chromosomes
merged_df['cumulative_pos'] = 0
chrom_spacing = 150_000_000  #equal spacing between chromosomes
chrom_offsets = {chrom: i * chrom_spacing for i, chrom in enumerate(sorted(merged_df['chromosome'].unique()))}

#assigns x-axis positions based on chromosome offsets
merged_df['cumulative_pos'] = merged_df['position'] + merged_df['chromosome'].map(chrom_offsets)

#save the merged data to a CSV file
merged_df.to_csv('mergeddata.csv', index=False)

#color palette
colors = sns.color_palette("bright", n_colors=len(merged_df['file_type'].unique()))

sns.set(rc={"figure.figsize": (12, 6)})  # width=12, height=6 for better spacing

g = sns.scatterplot(
    x=merged_df['cumulative_pos'], 
    y=merged_df['-logp'], 
    hue=merged_df['file_type'],  #colors by file type
    palette=colors, 
    data=merged_df, 
    legend='full',
    edgecolor="white",  #adds a white border around points
    linewidth=0.5       #sets border thickness
)

#move the legend (file type box) away from the graph
g.legend(loc='upper left', bbox_to_anchor=(1.05, 1), title="File Type")

g.set_xlabel('Chromosome')
g.set_ylabel('-log10(p-value)')
g.set_title('Manhattan Plot of TWAS Results')

#finds the median positions for chromosome labels with equal spacing
xticks = [chrom_offsets[chrom] + chrom_spacing / 2 for chrom in sorted(chrom_offsets.keys())]
g.set_xticks(xticks)
g.set_xticklabels([int(chrom) for chrom in sorted(chrom_offsets.keys())])  #convert chromosome labels to integers

#adds a significance threshold line
g.axhline(y=-np.log10(5e-8), color='red', linestyle='--')

fig = g.get_figure()
fig.savefig(output_file, dpi=300, bbox_inches='tight')  # Save with tight layout to avoid cutting off the legend

plt.show()

'''
so to run the script, you would use the command line like this:
python testManhat.py -i file1.csv file2.csv file3.csv -o output_plot.png

example:
python /home/project4/testManhat.py 
-i /home/project4/MetaXcan/software/CD4-positive_alpha-beta_cytotoxic_T_cell_output.csv 
/home/project4/MetaXcan/software/CD4-positive_alpha-beta_T_cell_output.csv 
/home/project4/MetaXcan/software/CD8-positive_alpha-beta_T_cell_output.csv 
/home/project4/MetaXcan/software/CD14-low_CD16-positive_monocyte_output.csv 
/home/project4/MetaXcan/software/CD14-positive_monocyte_output.csv 
/home/project4/MetaXcan/software/CD16-negative_CD56-bright_natural_killer_cell_human_output.csv 
/home/project4/MetaXcan/software/central_memory_CD4-positive_alpha-beta_T_cell_output.csv 
/home/project4/MetaXcan/software/central_memory_CD8-positive_alpha-beta_T_cell_output.csv 
/home/project4/MetaXcan/software/conventional_dendritic_cell_output.csv 
/home/project4/MetaXcan/software/dendritic_cell_output.csv 
-o /home/project4/manhattan_plot.png

-sara
'''
