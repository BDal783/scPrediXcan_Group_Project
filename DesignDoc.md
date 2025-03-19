# Overview
Genome-wide association studies (GWAS) is a method of study used to identify associations between specific genes and diseases and/or traits through analysing the frequency of single nucleotide polymorphisms (SNPs) in the genomes of a target population (E.x. diabetes patients) (MedlinePlus). Comparing the frequency of SNPs in the target population against a general population highlights commonly occurring differences in the genomes of those with the target disease/trait. Transcriptome-wide association studies (TWAS) take GWAS a step further by utilizing gene datasets to showcase gene expression levels in relation to the trait/disease in focus (Mai et al 2023). The advantage of TWAS is that it detects functioning genes and allows analysis of the mechanisms behind the trait/disease being studied and can even lead to determining disease prediction/risk. For TWAS to function, access to an expression quantitative trait loci (eQTL) is necessary for the creation of a predictive model and to act as a reference sheet to help identify genes of interest. <br />
One issue that scientists ran into when trying to access single cell data was the lack of full data. To overcome this, a deep learning algorithm and scRNA-seq data was used to create a gene expression prediction model as a method to predict epigenetic features at cell-type level, ctPred (Zhous et al 2024). This method isn’t restricted to the tissue-level predictions allowed by analysis of bulk RNA-seq data, as it uses cell specific expression and gene behavior to predict epigenetic traits such as cellular behavior under different environmental factors. scPrediXcan is a framework that utilizes ctPred and single-cell data (GWAS) to perform a TWAS. scPrediXcan allows ctPred to work with the scale of GWAS information by reducing the computational demand that ctPred would otherwise be unable to handle. scPrediXcan works in three steps: utilizing scPred to predict gene expression, converting the information into a linear model based on SNPs, and testing associations between the genes and the disease/trait by analyzing the data alongside GWAS summary statistics. <br />
Sources: <br />
Alquicira-Hernandez J, Sathe A, Ji HP, Nguyen Q, Powell JE. scPred: accurate supervised method for cell-type classification from single-cell RNA-seq data. Genome Biol. 2019 Dec 12;20(1):264. doi: 10.1186/s13059-019-1862-5. PMID: 31829268; PMCID: PMC6907144.<br />
Mai, J., Lu, M., Gao, Q. et al. Transcriptome-wide association studies: recent advances in methods, applications and available databases. Commun Biol 6, 899 (2023). https://doi.org/10.1038/s42003-023-05279-y<br />
MedlinePlus. Bethesda (MD): National Library of Medicine (US); What are genome-wide association studies?. Available from: https://medlineplus.gov/genetics/understandi ng/genomicresearch/gwastudies/<br />
Wainberg M, Sinnott-Armstrong N, Mancuso N, Barbeira AN, Knowles DA, Golan D, Ermel R, Ruusalepp A, Quertermous T, Hao K, Björkegren JLM, Im HK, Pasaniuc B, Rivas MA, Kundaje A. Opportunities and challenges for transcriptome-wide association studies. Nat Genet. 2019 Apr;51(4):592-599.<br />
Zhous Y, Adeluwa T, Zhu L, Salazar-Magaña S, Sumner S, Kim H, Gona S, Nyasimi F, Kulkarni R, Powell J, Madduri R, Liu B, Chen M, Im HK. scPrediXcan integrates advances in deep learning and single-cell data into a powerful cell-type-specific transcriptome-wide association study framework. bioRxiv [Preprint]. 2024 Nov 14:2024.11.11.623049.<br />
<br />
# Context
The basis of this model will focus on three central steps. This is training a model, creating a linear learning model of ctPred for SNP data usage, and trait association. However, through the model developed by the Hakyim lab that has already completed the creation and training of the model. This leaves us with developing a Python script that will perform the trait associations. We will start by using GWAS data from a Polycystic ovary syndrome (PCOS) study, but we plan on having this program function properly with any GWAS and produce an accurate TWAS result. <br />

# Goals
1. Build a Python pipeline that runs Step 3: Trait association of the scPrediXcan software with previously trained cell type specific prediction models (steps 1-2, linked from their repo) and GWAS summary statistics.
2. Create a flexible pipeline that works for any GWAS summary statistics.
3. Automatically generate a Manhattan plot for the scPrediXcan cell types. 
4. Adjust code so newest versions of numpy and pandas can be used.
![](https://www.biorxiv.org/content/biorxiv/early/2024/11/14/2024.11.11.623049/F1.large.jpg)
# Non-Goals
1. Label expression value problem, where cells have an expression value of 0 and are not ranked in the same order each run. In a one cell type/one model, there will not be ranking disparities.

# Proposed Solution
1. Setting up the environment → Downloading models from scPrediXcan GitHub repository, GWAS summary statistics, sample data that is required for the initial runs, installing any dependencies and importing all libraries
2. Preparing the Data → Checking the dataset structure and converting or unzipping the files, including argparse for the code to run on different types of data
3. Running the pipeline → Running S-PrediXcan and scPrediXcan to find the associations between genes and traits
4. Result Interpretation → Creating gene association tables, Manhattan plots and identifying significant gene-trait relationships
5. Documentation → Writing all steps into the READ.me file. Documenting the installation steps such as dependencies, providing explanation of where training model and sample data was found from, then noting step-by-step instructions on how to run the pipeline.

# Milestones
| Week | Date | Assignment | Sara | Hannah | Brendon | Jude |
|------|------|------------|------|--------|---------|------|
|Week 1|3/12|Work on Design Doc|Implement Plan Slides. <br />Proposed solution.|Milestones and general slide format. Goals, non-goals.|Introduction- Context.|Introduction- Overview.|
|Week 2|3/19|Design Doc + Presentation|Apply feedback to ‘Proposed solution’. Ensure code functionality.|Apply feedback to ‘Goals /NonGoals’ and ‘Milestones’. Edit base code to use arg parse.|Apply feedback to ‘Intro- Context’. Edit base code to exclude np.object.|Apply feedback to ‘Intro- Overview’. Ensure code functionality|
|Week 3|3/26|Repo Check #1|Begin implementing automatic Manhattan Plot production.|Modify code to function with specific GWAS headers present in sample data (or common).|Begin implementing automatic Manhattan Plot production.|Begin implementing automatic Manhattan Plot production.|
|Week 4|4/2|5-min progress presentation|Debug automatic Manhattan Plot code using smallest GWAS sample data|Debug automatic Manhattan Plot code using smallest GWAS sample data|Debug automatic Manhattan Plot code using smallest GWAS sample data|Debug automatic Manhattan Plot code using smallest GWAS sample data|
|Week 5|4/9|Repo Check #2|Check functionality of entire pipeline with other GWAS sample data. Debug if needed. |Check functionality of entire pipeline with other GWAS sample data. Debug if needed. |Check functionality of entire pipeline with other GWAS sample data. Debug if needed. |Check functionality of entire pipeline with other GWAS sample data. Debug if needed.|
|Week 6|4/16|5-min progress presentation and Rough Draft|Begin final presentation.|Begin final presentation.|Begin final presentation.|Begin final presentation.|
|Week 7|4/23|Final Presentation and Repo Check #3|Repo polishing|Code polishing and optimization.|Code polishing and optimization.|Repo polishing|
|Week 8|4/30|Final Project Code and Final Application Note|||||
