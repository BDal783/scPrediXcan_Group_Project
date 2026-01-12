import sys, os 
import argparse
import difflib
import gzip
import urllib.request

# function to parse command line arguments
def check_arg(args=None):
    parser = argparse.ArgumentParser(
    description="SPrediXcan Cell Data Processing Script.")
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

# link for the basic gene annotation GTF file
gtf_url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.basic.annotation.gtf.gz"

output_path = "./gencode_bga.gtf.gz"

# check if the file exists
if not os.path.exists(output_path):
    print(f"File not found locally. Downloading from {gtf_url}...")
    try:
        # download the file
        urllib.request.urlretrieve(gtf_url, output_path)
        print(f"Download complete: {output_path}")
    except Exception as e:
        print(f"An error occurred while downloading the file: {e}")
else:
    print(f"File already exists: {output_path}")

# Open and read the contents of a .gz file
with gzip.open(infile, 'rt') as file:
    head = file.readline()
    header = " ".join(head.split())
    headers_list = header.split(" ")

# known column names
known_columns = ["snp_column", "effect_allele_column", "non_effect_allele_column", "beta_column", "pvalue_column", "se_column", "freq_column", "position_column", "chromosome_column", "gwas_h2", "gwas_N", "No match found" ] #deleted locus column

# GWAS headers to match
gwas_headers = headers_list 

explicit_mappings = { #list of known column mappings
    "P": "pvalue_column",
    "#CHROM": "chromosome_column",
    "CHR": "chromosome_column",
    "FREQ": "freq_column",
    "AF_Allele2": "freq_column",
    'BETA': 'beta_column', 
    'A1':'effect_allele_column',
    'A2': 'non_effect_allele_column',
    'REF': 'non_effect_allele_column', 
    'ALT': 'effect_allele_column', 
    'SE': 'se_column',
    'Pvalue': 'pvalue_column',
    'BP': 'position_column',
    'POS': 'position_column',
    'SNP': 'snp_column',
    'SNP_hg38': 'snp_column',
    'ID': 'snp_column',
    'how': "No match found",
    'h2': 'gwas_h2',
    'heritability': 'gwas_h2',
    'sample_size': 'gwas_N',
    'N': 'gwas_N',
    'alleles': "No match found" #metaxcan doesn't recognize 'alleles column'
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

# snp_column formating to match .db files format: chr_pos_ref_alt_db
def format_snp_row(row, headers_mapping, genome_build="b38"): 
    # get mapped column names  
    snp_column_header = next((key for key, value in headers_mapping.items() if value == "snp_column"), None)
    ref_column_header = next((key for key, value in headers_mapping.items() if value == "non_effect_allele_column"), None)
    alt_column_header = next((key for key, value in headers_mapping.items() if value == "effect_allele_column"), None)
    chr_column_header = next((key for key, value in headers_mapping.items() if value == "chromosome_column"), None)
    pos_column_header = next((key for key, value in headers_mapping.items() if value == "position_column"), None)

    # can't work if snp header doesn't exists
    if not snp_column_header or snp_column_header not in row:
        print("Error: SNP column not found in row")
        return "unknown_format"

    snp_value = row[snp_column_header]

    # if the snp col contains "_" it is likely the same format as the .db files.  
    if "_" in snp_value:
        if snp_value.endswith(genome_build):  # already includes the genome build  
            formatted_snp = snp_value  # leave alone  
        else:  # add the genome build if missing  
            formatted_snp = f"{snp_value}_{genome_build}"
    elif ":" in snp_value:  
        parts = snp_value.split(":")  # split the SNP string by colon  
        if len(parts) == 2 and ref_column_header and alt_column_header:  # Format: chr:pos  
            ref = row.get(ref_column_header, "")
            alt = row.get(alt_column_header, "")
            formatted_snp = f"{parts[0]}_{parts[1]}_{ref}_{alt}_{genome_build}"
        elif len(parts) == 4:  # Format: chr:pos:ref:alt  
            formatted_snp = f"{parts[0]}_{parts[1]}_{parts[2]}_{parts[3]}_{genome_build}"
        elif len(parts) == 5:  # could be Format: chr:pos:ref:alt:db  
            formatted_snp = f"{parts[0]}_{parts[1]}_{parts[2]}_{parts[3]}_{parts[4]}"
    elif "rs" in snp_value and chr_column_header and pos_column_header and ref_column_header and alt_column_header:  # RSID format  
        chromosome = row.get(chr_column_header, "") # builds a proper snp format for .db compatibilty
        position = row.get(pos_column_header, "")
        ref = row.get(ref_column_header, "")
        alt = row.get(alt_column_header, "")
        formatted_snp = f"chr{chromosome}_{position}_{ref}_{alt}_{genome_build}"
    else:
        formatted_snp = "unknown_format"  # unexpected formats  
        print(f"Warning: Unrecognized SNP format: {snp_value}")
    return formatted_snp

# opens and reads input GWAS file, formats snp column values, and writes the altered data to a new gz file
def process_gwas_file(infile, outfile, headers_mapping, genome_build="b38"):
    
    gwas_alt_folder = "./GWAS_alt" # GWAS_alt folder path

    # create the folder only if it doesn't exist
    if not os.path.exists(gwas_alt_folder):
        os.makedirs(gwas_alt_folder)

    # new output file path based on '-o' argument
    alt_outfile = os.path.join(gwas_alt_folder, f"{outfile}.gz")

    # open the input GWAS file for processing and outfile for writing
    with gzip.open(infile, "rt") as gwas_file, gzip.open(alt_outfile, "wt") as out_file:
        # read and modify headers
        headers = gwas_file.readline().strip().split()
        headers.append("new_snp_column")  # add formatted SNP column
        out_file.write("\t".join(headers) + "\n")  # write updated headers

        # process each row and format SNP value 
        for line in gwas_file:
            row = dict(zip(headers, line.strip().split()))
            row["new_snp_column"] = format_snp_row(row, headers_mapping, genome_build)  # format SNP
            out_file.write("\t".join([str(row.get(h, "")) for h in headers]) + "\n")  # write updated row

    print(f"New modified GWAS file saved to: {alt_outfile}")

process_gwas_file(infile, outfile, mapped_columns)

# make the column info portion of the command using mapping from above, not including any "no match found" values
info_c = ""
for known_column in known_columns:
    matching_header = [key for key, value in mapped_columns.items() if value == known_column and value != "No match found"]
    if matching_header:
        if known_column == "snp_column": #skips over the orginal snp column so the formatted version can be used
            continue
        info_c += f"--{known_column} \\{matching_header[0]} \\\n" #added \\ so it takes # as a str not a command
info_c += f"--snp_column new_snp_column \\\n" #adds the formated snp column to the command.

# runs os.sys command and prints for debugging purposes
def run_bash(command): 
    print(f"Running command: {command}")
    os.system(command)

#get path from infile
gwas_path = os.path.dirname(os.path.abspath(infile)) 

#full list of cell data files
file_list = [
    "CD14-low_CD16-positive_monocyte_covariances.txt.gz",
    "CD14-low_CD16-positive_monocyte.db",
    "CD14-positive_monocyte_covariances.txt.gz",
    "CD14-positive_monocyte.db",
    "CD16-negative_CD56-bright_natural_killer_cell_human_covariances.txt.gz",
    "CD16-negative_CD56-bright_natural_killer_cell_human.db",
    "CD4-positive_alpha-beta_cytotoxic_T_cell_covariances.txt.gz",
    "CD4-positive_alpha-beta_cytotoxic_T_cell.db",
    "CD4-positive_alpha-beta_T_cell_covariances.txt.gz",
    "CD4-positive_alpha-beta_T_cell.db",
    "CD8-positive_alpha-beta_T_cell_covariances.txt.gz",
    "CD8-positive_alpha-beta_T_cell.db",
    "central_memory_CD4-positive_alpha-beta_T_cell_covariances.txt.gz",
    "central_memory_CD4-positive_alpha-beta_T_cell.db",
    "central_memory_CD8-positive_alpha-beta_T_cell_covariances.txt.gz",
    "central_memory_CD8-positive_alpha-beta_T_cell.db",
    "conventional_dendritic_cell_covariances.txt.gz",
    "conventional_dendritic_cell.db",
    "dendritic_cell_covariances.txt.gz",
    "dendritic_cell.db",
    "double_negative_thymocyte_covariances.txt.gz",
    "double_negative_thymocyte.db",
    "effector_memory_CD4-positive_alpha-beta_T_cell_covariances.txt.gz",
    "effector_memory_CD4-positive_alpha-beta_T_cell.db",
    "effector_memory_CD8-positive_alpha-beta_T_cell_covariances.txt.gz",
    "effector_memory_CD8-positive_alpha-beta_T_cell.db",
    "erythrocyte_covariances.txt.gz",
    "erythrocyte.db",
    "gamma-delta_T_cell_covariances.txt.gz",
    "gamma-delta_T_cell.db",
    "hematopoietic_precursor_cell_covariances.txt.gz",
    "hematopoietic_precursor_cell.db",
    "innate_lymphoid_cell_covariances.txt.gz",
    "innate_lymphoid_cell.db",
    "memory_B_cell_covariances.txt.gz",
    "memory_B_cell.db",
    "mucosal_invariant_T_cell_covariances.txt.gz",
    "mucosal_invariant_T_cell.db",
    "naive_B_cell_covariances.txt.gz",
    "naive_B_cell.db",
    "naive_thymus-derived_CD4-positive_alpha-beta_T_cell_covariances.txt.gz",
    "naive_thymus-derived_CD4-positive_alpha-beta_T_cell.db",
    "naive_thymus-derived_CD8-positive_alpha-beta_T_cell_covariances.txt.gz",
    "naive_thymus-derived_CD8-positive_alpha-beta_T_cell.db",
    "natural_killer_cell_covariances.txt.gz",
    "natural_killer_cell.db",
    "peripheral_blood_mononuclear_cell_covariances.txt.gz",
    "peripheral_blood_mononuclear_cell.db",
    "plasmablast_covariances.txt.gz",
    "plasmablast.db",
    "plasmacytoid_dendritic_cell_covariances.txt.gz",
    "plasmacytoid_dendritic_cell.db",
    "platelet_covariances.txt.gz",
    "platelet.db",
    "regulatory_T_cell_covariances.txt.gz",
    "regulatory_T_cell.db",
    "transitional_stage_B_cell_covariances.txt.gz",
    "transitional_stage_B_cell.db"
]

#seperates files based on type indicators
covariance_files = [f for f in file_list if f.endswith("_covariances.txt.gz")]
db_files = [f for f in file_list if f.endswith(".db")]

cell_type_tracker = [] #tracks all cell types used
base_path = os.path.abspath('./CellData') #finds the base path for CellData in the users system, used to give metaxcan the full path
gwas_path = os.path.abspath("./GWAS_alt")
os.chdir('./MetaXcan/software/') # change directory so metaxcan can be called 

# loop through each covariance file and database file
for cov_file in covariance_files:
    # extract the cell type from covariance file name
    cell_type = cov_file.replace("_covariances.txt.gz", "")
    cell_type_tracker.append(cell_type) #for future use calling Manhattenplot.py
    # find corresponding database file, this way nothing needs to be in order = less user error
    db_file = next((db for db in db_files if db.startswith(cell_type)), None)
    if db_file:
        # output file name based on cell type
        output_file = f"--output_file {cell_type}_output.csv"
        main_command = (
            f"python SPrediXcan.py \\\n"
            f"--model_db_path {base_path}/{db_file} \\\n"
            f"--covariance {base_path}/{cov_file} \\\n"
            f"--gwas_file {gwas_path}/{outfile}.gz \\\n"
            "--overwrite \\\n"
            "--chromosome all \\\n"
            "--keep_non_rsid \\\n"
        )
        command = main_command + info_c + output_file #construct full command
        #run the command
        run_bash(command)
    else:
        print(f"No matching database file found for {cell_type}.") #debugging


per_graph = 10 # 10 cell types per plot
main_folder = os.path.dirname(gwas_path) 
os.chdir(main_folder) #moving back to the project folder
manhat_name = os.path.basename(infile).replace(".txt.gz", "")
manhat_folder = "./ManhattanPlots" # folder to store manhattan plots
if not os.path.exists(manhat_folder): #make folder if it doesnt exist
    os.makedirs(manhat_folder)

for i in range(0, len(cell_type_tracker), per_graph):
    plot = cell_type_tracker[i:i+per_graph] #10 cell types to plot 
    inputs = " ".join([os.path.abspath(f"./MetaXcan/software/{cell}_output.csv") for cell in plot]) # get csv files using output naming convention
    outputs = os.path.join(manhat_folder, f"{manhat_name}_{i//per_graph +1}.png") #names output using the input GWAS for organization
    man_command = f"python Manhattanplot.py -i {inputs} -o {outputs}"
    run_bash(man_command)

# last edit: 4/30/2025 2:00pm by Hannah
