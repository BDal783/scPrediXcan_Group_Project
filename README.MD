Comp 383 scPrediXcan Group Project

1. Quick Background
-The purpose of this repo was to take the scPrediXcan code created by hakyimlab, be able to run it and improve upon it. We are doing this by improving the way that it reads in GWAS files. This by using such tools like difflib and it being able to read in multiple GWAS files. We also improved upon the file by adding to the TWAS Results the gene's chromosome and SNP, this so that it automatically creates a Manhatten plot of the TWAS results.
2. Installation
To begin installation one need to first creat a copy of this repo by using git clone.

```bash
git clone https://github.com/BDal783/scPrediXcan_Group_Project.git
```

After finishing that install one needs to set up the conda enviornment to run the scPrediXcan code. If one dosent have conda installed make sure to cd to your personal folder and install miniconda. (I'll include a current version example but make sure to download the most up to date version)

```bash
 mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh
~/miniconda3/bin/conda init bash
~/miniconda3/bin/conda init zsh
```

Once conda is installed make sure to cd into the miniconda folder and to activate it using "conda activate".

Now with conda installed you need to cd into the project folder and create a conda enviornment. The code can be seen below.

```bash
cd scPrediXcan_Group_Project/MetaXcan/software
conda env create -f conda_env.yaml
conda activate imlabtools
```
```bash
run this test code: ./SPrediXcan.py --model_db_path 'data/DGN-WB_0.5.db' --covariance 'data/covariance.DGN-WB_0.5.txt.gz' --gwas_folder data/GWAS --gwas_file_pattern ".*gz" --snp_column SNP --effect_allele_column A1 --non_effect_allele_column A2 --beta_column BETA --pvalue_column P --output_file 'TWAS_result.csv'
```

3. Dependencies
-
4. Running the Python Wrapper
-
5. Acknowledgments

