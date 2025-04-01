import sys, os 
import argparse
import difflib
<<<<<<< HEAD
import gzip

# function to parse command line arguments
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


# Open and read the content of a .gz file
with gzip.open(infile, 'rt') as file:
    head = file.readline()
    header = " ".join(head.split())
    headers_list = header.split(" ")


# known column names
known_columns = ["snp_column", "effect_allele_column", "non_effect_allele_column", "reference_allele_column", "other_allele_column", "alleles_column", "beta_column", "pvalue_column", "se_column", "freq_column", "position_column", "chromosome_column", 'locus_column', "No match found" ]

# GWAS headers to match
gwas_headers = headers_list 
=======

#function to parse command line arguments
# def check_arg(args=None):
#     parser = argparse.ArgumentParser(
#     description="ADD TITLE OF SCRIPT HERE (shows on help -h)")
#     parser.add_argument("-i", "--input",
#     help="input file",
#     required=True)
#     parser.add_argument("-o", "--output",
#     help="output file",
#     required=True)
#     return parser.parse_args(args)

# #retrieve command line arguments
# arguments = check_arg(sys.argv[1:])
# infile = arguments.input
# outfile = arguments.output


#known column names
known_columns = ["snp_column", "effect_allele_column", "non_effect_allele_column", "reference_allele_column", "other_allele_column", "alleles_column", "beta_column", "pvalue_column", "se_column", "freq_column", "position_column", "chromosome_column", 'locus_column', "No match found" ]

#GWAS headers to match
gwas_headers = ["CHR", "BP", "RS", "SNP_hg38", "A1", "A2", "FREQ", "BETA", "SE", "P", 
 "locus", "how", "ID", "REF", "ALT", "#CHROM", "POS", 
 "Pvalue", "AF_Allele2", "alleles", "#CHROM", "Pvalue"]
>>>>>>> 0868a33121473d109c776ee161f6aee6d33400e8

explicit_mappings = {
    "P": "pvalue_column",
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
    'Pvalue': 'pvalue_column',
    'BP': 'position_column',
    'POS': 'position_column',
    'SNP': 'snp_column',
    'RS': 'snp_column',
    'how': 'No match found',
    }

<<<<<<< HEAD
# mapping headers to closest known column names
=======
#mapping headers to closest known column names
>>>>>>> 0868a33121473d109c776ee161f6aee6d33400e8
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

<<<<<<< HEAD
# print the mapping
print(mapped_columns)

# commands that use the dicionary of mapped columns  to assign header name:
# example: 
# snp_column SNP 
# effect_allele_column A1
# non_effect_allele_column A2
# then use os to call SPrediXcan.py with the necessary info/info from the GWAS
# basic info will include: SNP, effect allele, non effect allele, beta, pval. 
# and hope it works :) 

=======
# Print the mapping
print(mapped_columns)



#Reading in the results and creating a manhatten plot of the results section -Brendon
>>>>>>> 0868a33121473d109c776ee161f6aee6d33400e8
