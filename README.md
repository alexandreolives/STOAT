# STOAT

<p align="center">
    <a href="https://www.python.org/downloads/release/python-3100/"><img src="https://img.shields.io/badge/Python-3.10-blue.svg"></a>
    <a href="https://github.com/brentp/cyvcf2/releases/tag/v0.31.1"><img src="https://img.shields.io/badge/cyvcf2-0.31.1-green.svg"></a>
    <a href="https://github.com/vgteam/libbdsg/releases/tag/v0.3"><img src="https://img.shields.io/badge/bdsg-0.3-green.svg"></a>
    <a href="https://github.com/statsmodels/statsmodels/releases/tag/v0.14.4"><img src="https://img.shields.io/badge/statsmodels-0.14.4-green.svg"></a>
    <a href="https://github.com/ShujiaHuang/qmplot/releases/tag/v0.3.1"><img src="https://img.shields.io/badge/qmplot-0.3.3-green.svg"></a>
</p>

## Project Overview
STOAT is a specialized tool developed for conducting Genome-Wide Association Studies (GWAS) with a unique focus on snarl structures within pangenome graphs. Unlike traditional GWAS tools that analyze linear genome variants, STOAT processes VCF files to extract and analyze snarl regions—complex structural variations that capture nested and overlapping variant patterns within a pangenome. This approach allows for a more nuanced understanding of genetic variations in diverse populations and complex traits.

STOAT supports both binary and quantitative phenotypes:

- For binary phenotypes (e.g., case vs control studies), it utilizes chi-squared tests and Fisher’s exact test to evaluate associations between phenotype groups and snarl variants, providing robust statistical validation even in cases of sparse data.

- For quantitative phenotypes (e.g., traits measured on a continuous scale), the tool employs linear regression models to assess the association between snarl structures and phenotype values, allowing for continuous trait mapping with greater precision.

## Installation

````bash
git clone https://github.com/Plogeur/STOAT.git
cd STOAT
pip install -r requirements.txt

# install bdsg version > 3.0.0 
# DO NOT use pip install bdsg cause version 3.0.0 instead do :
git clone --recursive https://github.com/vgteam/libbdsg.git
cd libbdsg
pip install .

# see more installation information on bdsg github
````

## Dependencies
- Python version 3.10
- cyvcf2
- numpy
- pandas
- bdsg
- ...


## Input format file

Required files :
- pg_file : Pangenome graph file, formats accepted: .pg or .xg.
- dist_file : Distance file generated with vg dist, format: .dist.
- vcf_ref : Reference node positions in the pangenome graph generated with vg deconstruct, formats: .vcf.
- vcf : Merged VCF file, created using bcftools merge, formats: .vcf or .vcf.gz.
- pheno : Three-column file with FID (family/sample name), IID (sample name), and PHENO (integer/float). Format: .txt or .tsv (tab-separated).

Optional file : 
- list_paths : Two-column file containing snarl names and the list of paths through the snarl's netgraph, separated by tabs. Format: .txt or .tsv.

## Usage

Use `stoat.py` if you want to launch the full tool at once, starting from snarl path identification (identifying the multiple paths that can be taken by a sample based on the pangenome graph) and ending with the results plots (Manhattan plot and QQ plot).

- Run full tool :
```bash
# binary trait
python3 stoat.py -p <path_to_pg_file.pg> -d <path_to_dist_file.dist> -v <path_to_vcf_file.vcf.gz> -r <path_to_vcf_reference_file.vcf.gz> -b <path_to_group_file.txt> -o output.tsv

# quantative trait
python3 stoat.py -p <path_to_pg_file.pg> -d <path_to_dist_file.dist> -v <path_to_vcf_file.vcf.gz> -r <path_to_vcf_reference_file.vcf.gz> -q <path_to_pheno_file.txt> -o output.tsv
```

Explanation of all options:
```bash 
-p           Path to the input pangenome `.pg` file (mutually exclusive, required if `-l` not provided).
-d           Path to the input distance index `.dist` file (mutually exclusive, required if `-l` not provided).
-t           Specifies the children threshold (optional).
-v           Path to the merged VCF file (`.vcf` or `.vcf.gz`) (required).
-r           Path to the VCF file referencing all snarl positions (`.vcf` or `.vcf.gz`) (mutually exclusive, required if `-p and -d` not provided).
-l, --listpath  
             Path to the list of paths file (optional).
-b, --binary  
             Path to the binary group file (`.txt` or `.tsv`) (mutually exclusive, required if `-q` not provided).
-q, --quantitative  
             Path to the quantitative phenotype file (`.txt` or `.tsv`) (mutually exclusive, required if `-b` not provided).
-c, --covariate  (working progress)
             Path to the covariate file (`.txt` or `.tsv`), LMM analysis (optional).
-g, --gaf    
             Prepare binary GWAS output for a GAF file and create a GAF file with the top 10 significant paths (optional).
-o, --output  
             Specifies the base path for the output directory (optional).
```

Alternatively, you can specify the script you want to launch, depending on your desired task:

- **`list_snarl_paths.py`**: Identifies snarl paths. (By default, it excludes snarls with more than 50 children or 2 second computation by snarl)  
  **Input**: Graph information files (`.dist` and `.pg`)  
  **Output**: A two-column string file (`.tsv`) containing the snarl reference and the corresponding list of paths.

- **`children_distribution.py`**: An optional script to help decide the thresholds for `list_snarl_paths.py` based on the children distribution in the netgraph (pangenome graph).  
  **Input**: Graph information files (`.dist` and `.pg`)  
  **Output**: None

- **`snarl_analyser.py`**: The main script. It parses and analyzes all input files to return p-value associations.  
  **Input**:  
  - Snarl list paths (`.tsv` or `.txt`)  
  - Binary or quantitative phenotype (`.tsv` or `.txt`)  
  - Merged VCF (`.vcf.gz` or `.vcf`)  
  **Output**: P-value associations.

- **`gaf_creator.py`**: Creates a GAF file from a list of binary p-value associations (output of `snarl_analyser.py`).  
  **Input**: Binary p-value association file.  
  **Output**: GAF file.

- **`p_value_analysis.py`**: Creates a Manhattan plot and QQ plot (output of `snarl_analyser.py`).  
  **Input**: Binary or Quantitative p-value association file.  
  **Output**: 2 PNG file plots.

- Example of specific script of tool : 

```bash
# decompose pangenome
python3 list_snarl_paths.py -p <path_to_pg_file.pg> -d <path_to_dist_file.dist> -o <output.tsv>

# binary trait with list_path already computed (can also do with quantitative trait) and gaf creation 
python3 stoat.py -l <list_paths_snarl.txt> -v <path_to_vcf_file.vcf.gz> -r <path_to_vcf_reference_file.vcf.gz> -b <path_to_group_file.txt> --gaf -o output.tsv

# binary trait
python3 snarl_analyser.py <path_to_vcf_file.vcf.gz> <list_paths_snarl.txt> <path_to_vcf_reference_file.vcf.gz> -b <path_to_group_file.txt> -o output.txt

# quantitative trait 
python3 snarl_analyser.py <path_to_vcf_file.vcf.gz> <list_paths_snarl.txt> <path_to_vcf_reference_file.vcf.gz> -q <path_to_pheno_file.txt> -o output.txt
```

## Output

| Column Name       | Description                                                                                   |
|-------------------|-----------------------------------------------------------------------------------------------|
| **CHR**           | Chromosome name where the variation occurs.                                                   |
| **POS**           | Position of the snarl within the chromosome.                                                  |
| **SNARL**         | Identifier for the variant, snarl name/id                                                     |
| **TYPE**          | Type of genetic variation (e.g., SNP, INS, DEL).                                              |
| **REF**           | Sequence of the reference allele.                                                             |
| **ALT**           | Sequence of the alternate allele.                                                             |
| **P_FISHER**      | P-value calculated using Fisher's exact test (binary analysis).                               |
| **P_CHI2**        | P-value calculated using the Chi-squared test (binary analysis).                              |
| **TOTAL_SUM**     | Total number of samples that pass for this snarl. (binary analysis).                          |
| **MIN_ROW_INDEX** | Minimum group of samples that pass through one path of the snarl. (binary analysis).          |
| **NUM_COLUM**     | Number of paths in the snarl. (binary analysis).                                              |
| **INTER_GROUP**   | Sum of the minimum samples that pass through each path. (binary analysis).                    |
| **AVERAGE**       | Average number of total samples passing through this snarl, divided by the number of paths. 
(binary analysis).  |
| **P**             | P-value calculated using linear regression (quantitative analysis).                           |
| **RSQUARED**      | R-squared value indicating the proportion of variance explained by the model (quantitative analysis).          |
| **SE**            | Standard error of the estimated coefficients in the linear regression(quantitative analysis). |
| **BETA**          | Beta coefficients representing the estimated effect sizes of the predictors(quantitative analysis).
|

### Example of Output:

Below is an example of the output for a binary phenotype analysis:

```bash
CHR POS SNARL           TYPE  REF ALT   P_FISHER  P_CHI2  TOTAL_SUM  MIN_ROW_INDEX NUM_COLUM INTER_GROUP AVERAGE
1   12  5262721_5262719 SNP   A   T     0.4635    0.5182  286        2             137       46          143.0
1   15  5262719_5262717 INS   A   ATT   0.8062    0.8747  286        2             141       34          143.0
1   18  5262717_5262714 DEL   AA  T     0.2120    0.2363  286        2             134       32          143.0
```

Below is an example of the output for a quantitative phenotype analysis (-q option) :

```bash
CHR	POS	SNARL	          TYPE	  REF	ALT	RSQUARED	  BETA	      SE	        P
1   12  5262721_5262719	SNP     A   T   8.3697e-01	1.3878e+01	6.5108e+00	4.0376e-01
1   15  5262719_5262717	INS     A   ATT 4.4237e-01	1.3238e+01	6.5345e+00	4.6574e-01
1   18  5262717_5262714	DEL     AA  T   6.3237e-01	1.6458e+01	6.6453e+00	4.7484e-01
1   19  5262717_5262714	COMPLEX C   NA  4.2342e-01	2.3242e+01	5.3251e+00	1.3245e-01
```

## Visualization

### Manhattan and QQ plots 

STOAT will generated a manhattan and a QQ plot for binary and quantitatif analysis.

### SequenceTube

Use `gaf_creator.py` to geneate a GAF file and [SequenceTube](https://github.com/vgteam/sequenceTubeMap) tool to visualize your gwas binary region results.

```bash 
python3 gaf_creator.py -s <binary_gwas_stoat_output.tsv> -l <decomposition_paths.tsv> -p <pg_file_path>
```

*Output Plots and sequenceTube*
