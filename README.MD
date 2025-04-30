# Comp 383 scPrediXcan Group Project

## Quick Background
The purpose of this repo was to take the scPrediXcan code created by hakyimlab, be able to run it and improve upon it. We are doing this by improving the way that it reads in GWAS files. This by using such tools like difflib and it being able to read in multiple GWAS files. We also improved upon the file by adding to the TWAS Results the gene's chromosome and SNP, this so that it automatically creates a Manhatten plot of the TWAS results.
## Installation
>[!CAUTION]
>This project uses packages only available for Linux and MacOS; attempting to create the necessary environment on Windows will fail.<br />

>[!WARNING]
>This project uses git and miniconda, which must be installed and set up first.<br />
>Miniconda here: https://www.anaconda.com/download/success<br />
>Git here: https://git-scm.com/downloads<br />

To begin installation, download a copy of this repo. This will download the project files to the current directory. 

```bash
git clone https://github.com/BDal783/scPrediXcan_Group_Project.git
```

After finishing the clone download, ensure conda is activated and you're cd'd into the project folder.

```bash
conda activate
cd ./scPrediXcan_Group_Project
```

Now, the conda environment needed to run scPrediXcan must be activated.

```bash
cd MetaXcan/software
conda env create -f conda_env.yaml
conda activate imlabtools2
```

You can run this code to test the scPrediXcan install: 

```bash
./SPrediXcan.py --model_db_path 'data/DGN-WB_0.5.db' --covariance 'data/covariance.DGN-WB_0.5.txt.gz' --gwas_folder data/GWAS --gwas_file_pattern ".*gz" --snp_column SNP --effect_allele_column A1 --non_effect_allele_column A2 --beta_column BETA --pvalue_column P --output_file 'TWAS_result.csv'
```
It will create a TWAS_result.csv in MetaXcan/software, if it runs into issues take debugging to the first acknowledgement link (scPrediXcan Github)


## Dependencies<br />
Creating the conda environment will download all necessary dependencies

## Running the Python Wrapper<br />
List your cell data files in the file list of Pipe.py starting at line 121 (MAKE SURE THIS REMAINS 121 IN FINAL VERSION)<br />
Run Pipe.py with a single GWAS sumstats file<br />

```bash
python Pipe3.py -i /path/to/gwas_sumstats_file/gwas.txt.gz -o altered_gwas_name.txt
```


-
## Acknowledgments
https://github.com/hakyimlab/scPrediXcan?tab=readme-ov-file<br />
https://python-graph-gallery.com/manhattan-plot-with-matplotlib/<br />
https://stackoverflow.com/questions/37463184/how-to-create-a-manhattan-plot-with-matplotlib-in-python<br />

