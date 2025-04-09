'''
this was the work i had put in before but because it did not focus on twas and instead gwas. I will be fixing this file to focus on the twas data using matplotlib and seaborn and ensuring that the file is able to work with multiple files.
the output file has been saved as output.csv
'''
import sys, os
import argparse
import difflib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt  # Corrected import from matplotlib as plt

#function to parse command line arguments
def check_arg(args=None):
     parser = argparse.ArgumentParser(
     description="ADD TITLE OF SCRIPT HERE (shows on help -h)")
     parser.add_argument("-i", "--input",
     help="input file",
     required=True)
     parser.add_argument("-o", "--output",
     help="output file",
     required=True)
     return parser.parse_args(args)

# retrieve command line arguments
arguments = check_arg(sys.argv[1:])
infile = arguments.input
outfile = arguments.output

#known column names
known_columns = ["snp_column", "effect_allele_column", "non_effect_allele_column", "reference_allele_column", "other_allele_column", "alleles_column", "beta_column", "pvalue_column", "se_column", "freq_column", "position_column", "chromosome_column", 'locus_column', "No match found" ]

#GWAS headers to match
gwas_headers = ["CHR", "BP", "RS", "SNP_hg38", "A1", "A2", "FREQ", "BETA", "SE", "P", 
 "locus", "how", "ID", "REF", "ALT", "#CHROM", "POS", 
 "Pvalue", "AF_Allele2", "alleles", "#CHROM", "Pvalue"]

explicit_mappings = {
    "P": "pvalue", #it wasnt working with pvalue_column because the file header is pvalue..
    "#CHROM": "chromosome_column",
    "CHR": "chromosome_column",
    "FREQ": "freq_column",
    "AF_Allele2": "freq_column",
    'BETA': 'beta_column', 
    'A1':'effect_allele_column',
    'A2': 'non_effect_allele_column',
    'REF': 'reference_allele',
    'ALT': 'other_allele',
    'SE': 'se_column',
    'Pvalue': 'pvalue',
    'BP': 'position_column',
    'POS': 'position_column',
    'SNP': 'snp_column',
    'RS': 'snp_column',
    'how': 'No match found',
}

#mapping headers to closest known column names
mapped_columns = {}
for header in gwas_headers:
    if header in explicit_mappings:
        mapped_columns[header] = explicit_mappings[header]
    else:
        matches = difflib.get_close_matches(header.lower(), [col.lower() for col in known_columns], n=1, cutoff=0.2)
        if matches:
            mapped_columns[header] = matches[0]
        else:
            mapped_columns[header] = "No match found"

# Print the mapping
print(mapped_columns)

#Reading in the results and creating a manhatten plot of the results section -Brendon

'''
https://python-graph-gallery.com/manhattan-plot-with-matplotlib/
https://stackoverflow.com/questions/37463184/how-to-create-a-manhattan-plot-with-matplotlib-in-python

This is the reference link I used as a template to make the plot. -sara
'''

#directory containing GWAS files
data_dir = "/home/project4/MetaXcan/software/data/GWAS/"
chromosome_files = [f"chr{i}.assoc.dosage.gz" for i in range(1, 23)]

#read and concatenate all chromosome files
dataframes = []
for file in chromosome_files:
    file_path = os.path.join(data_dir, file)
    if os.path.exists(file_path):
        df = pd.read_csv(file_path, delim_whitespace=True)
        df["chromosome_column"] = int(file.split(".")[0][3:])  # Extract chromosome number as integer
        dataframes.append(df)
    else:
        print(f"Warning: {file_path} not found. Skipping.")

#combining all chromosome data into one DataFrame
if dataframes:
    gwas_data = pd.concat(dataframes, ignore_index=True)
else:
    print("Error: No GWAS files were loaded.")
    sys.exit(1)

gwas_data = gwas_data.rename(columns=mapped_columns)

gwas_data['minuslog10pvalue'] = -np.log10(gwas_data['pvalue'])
gwas_data['chromosome_column'] = gwas_data['chromosome_column'].astype(str)
gwas_data['chromosome_column'] = pd.Categorical(gwas_data['chromosome_column'], categories=[str(i) for i in range(1, 23)], ordered=True)

gwas_data = gwas_data.sort_values('chromosome_column')
gwas_data['ind'] = range(len(gwas_data))
df_grouped = gwas_data.groupby('chromosome_column', observed=False)

fig = plt.figure(figsize=(14, 8))
ax = fig.add_subplot(111)
colors = ['darkred', 'darkgreen', 'darkblue', 'gold']
x_labels = []
x_labels_pos = []

for num, (name, group) in enumerate(df_grouped): #https://stackoverflow.com/questions/37463184/how-to-create-a-manhattan-plot-with-matplotlib-in-python
    group.plot(kind='scatter', x='ind', y='minuslog10pvalue', color=colors[num % len(colors)], ax=ax)
    x_labels.append(name)
    x_labels_pos.append((group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0]) / 2))

ax.set_xticks(x_labels_pos)
ax.set_xticklabels(x_labels)
ax.set_xlim([0, len(gwas_data)])
ax.set_ylim([0, 3.5])
ax.set_xlabel('Chromosome')
plt.savefig("manhattan_plot.png", dpi=300)
print("Plot saved as manhattan_plot.png")
plt.show()

gwas_data.to_csv(outfile, index=False)
print(f"{outfile}")

