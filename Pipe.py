import sys, os 
import argparse
import difflib
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

# mapping headers to closest known column names
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

# print the mapping
print(mapped_columns)

def run_bash(command): 
    print(f"Running command: {command}")
    os.system(command)

# using the mapping to create the command that calls SPrediXcan;
command = ""

# Iterate through all known columns and find the corresponding headers
for known_column in known_columns:
    matching_header = [key for key, value in mapped_columns.items() if value == known_column]
    if matching_header and known_column != "No match found":  # exclude "No match found"
        command += f"--{known_column} {matching_header[0]} \\\n"

# Print or execute the command
print(f"Generated command:\n{command}")
run_bash(command)
gwas_path = os.path.dirname(infile) #get path from infile
"--gwas_folder {gwas_path} \\n" #use path

os.chdir('MetaXcan/software') #moves directory so the main command cand find it's stuff

# main command that calls SPrediXcan.py and gives path info for models and data - taken from github
main_c = (
    "python /home/project4/scPrediXcan_Group_Project/MetaXcan/software/SPrediXcan.py \\\n"
    "--model_db_path data/DGN-WB_0.5.db \\\n"
    "--covariance data/covariance.DGN-WB_0.5.txt.gz \\\n"
    "--gwas_folder data/GWAS \\\n"   #removed MetaXcan/software/ b/c of addition of os.chdir above
    "--gwas_file_pattern \".*gz\" \\\n"
)

# make the column info portion of the command using mapping from above
info_c = ""
for known_column in known_columns:
    matching_header = [key for key, value in mapped_columns.items() if value == known_column]
    if matching_header and matching_header[0] != "No match found":
        info_c += f"--{known_column} {matching_header[0]} \\\n"

# final output location
output_c = f"--output_file {outfile}"

# Combine all parts into the full command
command = main_c + info_c + output_c

# Print or execute the full command
print(f"Generated full command:\n{command}")
run_bash(command)

#updated 4/2/2025 4:03pm By Jude