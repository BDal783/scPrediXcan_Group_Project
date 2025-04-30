# Comp 383 scPrediXcan Group Project

## Quick Background
The purpose of this repo is to take the scPrediXcan code created by hakyimlab and increase its user accessibility by automating the process of GWAS header mapping and Manhattan plot visualization of TWAS results. We are creating a pipeline file that reads in a GWAS sumstats file and reformats headers for processing, as there is no standard GWAS header format.
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
Move all cell model data files to the CellData folder.
List your cell model data files in the variable file_list of Pipe.py starting at line 190, make sure to remove any cell data files that aren't in use (MAKE SURE THIS REMAINS 121 IN FINAL VERSION)<br />
Run Pipe.py with a single GWAS sumstats file<br />

```bash
python Pipe3.py -i /path/to/gwas_sumstats_file/gwas.txt.gz -o altered_gwas_file.txt
```
Pipe.py reformats the gwas_sumstats into -o altered_gwas_file.txt for processing<br />
You can find the cell type outputs at Metascan/software named [cell_type].csv<br />
The pipeline pulls these outputs and creates Manhattan plots 10 cell types at a time in .pngs located in the ManhattanPlot folder<br />


## Acknowledgments
https://github.com/hakyimlab/scPrediXcan?tab=readme-ov-file<br />
https://python-graph-gallery.com/manhattan-plot-with-matplotlib/<br />
https://stackoverflow.com/questions/37463184/how-to-create-a-manhattan-plot-with-matplotlib-in-python<br />

