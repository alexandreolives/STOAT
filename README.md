# STOAT

<p align="center">
    <a href="https://github.com/STOAT" title="STOAT GitHub">
        <img src="https://img.shields.io/badge/STOAT-GWAS-blue.svg" />
    </a>
    <a href="https://www.python.org/downloads/release/python-310/" title="Python 3.10">
        <img src="https://img.shields.io/badge/Python-3.10-blue.svg" />
    </a>
    <a href="https://pypi.org/project/cyvcf2/0.31.1/" title="cyvcf2 version 0.31.1">
        <img src="https://img.shields.io/badge/cyvcf2-0.31.1-blue.svg" />
    </a>
    <a href="https://pypi.org/project/bdsg/0.3/" title="bdsg version 0.3">
        <img src="https://img.shields.io/badge/bdsg-0.3-purple.svg" />
    </a>
    <a href="https://pypi.org/project/statsmodels/0.14.4/" title="statsmodels version 0.14.4">
        <img src="https://img.shields.io/badge/statsmodels-0.14.4-yellow.svg" />
    </a>
    <a href="https://pypi.org/project/qmplot/0.3.3/" title="qmplot version 0.3.3">
        <img src="https://img.shields.io/badge/qmplot-0.3.3-red.svg" />
    </a>
</p>

## Project Overview
STOAT is a specialized tool developed for conducting Genome-Wide Association Studies (GWAS) with a unique focus on snarl structures within pangenome graphs. Unlike traditional GWAS tools that analyze linear genome variants, STOAT processes VCF files to extract and analyze snarl regions—complex structural variations that capture nested and overlapping variant patterns within a pangenome. This approach allows for a more nuanced understanding of genetic variations in diverse populations and complex traits.

STOAT supports both binary and quantitative phenotypes:

- For binary phenotypes (e.g., case vs. control studies), it utilizes chi-squared tests and Fisher’s exact test to evaluate associations between phenotype groups and snarl variants, providing robust statistical validation even in cases of sparse data.

- For quantitative phenotypes (e.g., traits measured on a continuous scale), the tool employs linear regression models to assess the association between snarl structures and phenotype values, allowing for continuous trait mapping with greater precision.

Conventional GWAS tools typically rely on a single-reference genome, often overlooking structural variants and variations present in underrepresented populations. By using a pangenome graph, STOAT enables analysis across multiple genomes simultaneously, capturing a broader spectrum of structural variations, including insertions, deletions, and nested variants. The focus on snarls—graph features that encapsulate complex variant structures—provides a powerful means to map associations that would be difficult to detect in linear-based analyses.

## Installation

````bash
git clone https://github.com/Plogeur/STOAT.git
cd STOAT
pip install -r requirements.txt
````

## Dependencies
- Python version 3.10
- cyvcf2
- numpy
- pandas
- ...

## Usage

Use `STOAT.py` if you want to launch the full tool at once, starting from snarl path identification (identifying the multiple paths that can be taken by a sample based on the pangenome graph) and ending with the results plots (Manhattan plot and QQ plot).

- Run full tool :
```bash
# binary trait
python3 STOAT.py -p <path_to_pg_file.pg> -d <path_to_dist_file.dist> -v <path_to_vcf_file.vcf.gz> -r <path_to_vcf_reference_file.vcf.gz> -b <path_to_group_file.txt> -o output.tsv

# quantative trait
python3 STOAT.py -p <path_to_pg_file.pg> -d <path_to_dist_file.dist> -v <path_to_vcf_file.vcf.gz> -r <path_to_vcf_reference_file.vcf.gz> -q <path_to_pheno_file.txt> -o output.tsv
```

Alternatively, you can specify the script you want to launch, depending on your desired task:

- **`list_snarl_paths.py`**: Identifies snarl paths. (By default, it excludes snarls with more than 50 children or 10,000 nodes.)  
  **Input**: Graph information files (`.dist` and `.pg`)  
  **Output**: A two-column string file (`.tsv`) containing the snarl reference and the corresponding list of paths.

- **`children_distribution.py`**: An optional script to help decide the thresholds for `list_snarl_paths.py` based on the children distribution in the netgraph (pangenome graph).  
  **Input**: Graph information files (`.dist` and `.pg`)  
  **Output**: None

- **`snarl_vcf_parser.py`**: The main script. It parses and analyzes all input files to return p-value associations.  
  **Input**:  
  - Snarl list paths (`.tsv` or `.txt`)  
  - Binary or quantitative phenotype (`.tsv` or `.txt`)  
  - Merged VCF (`.vcf.gz` or `.vcf`)  
  **Output**: P-value associations.

- **`gaf_creator.py`**: Creates a GAF file from a list of binary p-value associations (output of `snarl_vcf_parser.py`).  
  **Input**: Binary p-value association file.  
  **Output**: GAF file.

- **`p_value_analysis.py`**: Creates a Manhattan plot and QQ plot (output of `snarl_vcf_parser.py`).  
  **Input**: Binary or Quantitative p-value association file.  
  **Output**: 2 PNG file plots.

- Example of specific script of tool : 
```bash
# decompose pangenome
python3 list_snarl_paths.py -p <path_to_pg_file.pg> -d <path_to_dist_file.dist> -o <output.tsv>

# binary trait
python3 snarl_vcf_parser.py <path_to_vcf_file.vcf.gz> <path_to_snarl_file.txt> <path_to_vcf_reference_file.vcf.gz> -b <path_to_group_file.txt> -o output.txt

# quantitative trait 
python3 snarl_vcf_parser.py <path_to_vcf_file.vcf.gz> <path_to_snarl_file.txt> <path_to_vcf_reference_file.vcf.gz> -q <path_to_pheno_file.txt> -o output.txt
```

## Input format file

Required files :
- pg_file : Pangenome graph file, formats accepted: .pg or .xg.
- dist_file : Distance file generated with vg dist, format: .dist.
- vcf_ref : Reference node positions in the pangenome graph generated with vg deconstruct, formats: .vcf or .vcf.gz.
- vcf : Merged VCF file, created using bcftools merge, formats: .vcf or .vcf.gz.
- group/pheno : 
    - Binary phenotype: Two-column file with SAMPLE and GROUP labels. Format: .txt or .tsv (tab-separated).

    - Quantitative phenotype: Three-column file with FID (family/sample name), IID (sample name), and PHENO (integer/float). Format: .txt or .tsv (tab-separated).

Optional file : 
- snarl : Two-column file containing snarl names and the list of paths through the snarl's netgraph, separated by tabs. Format: .txt or .tsv.

## Output


## Output

| Column Name       | Description                                                                                   |
|-------------------|-----------------------------------------------------------------------------------------------|
| **CHR**           | Chromosome number where the variation occurs.                                                 |
| **POS**           | Position of the variation within the chromosome.                                              |
| **SNARL**         | Identifier for the variant region, specifying start and end positions.                        |
| **TYPE**          | Type of genetic variation (e.g., SNP, INS, DEL).                                              |
| **REF**           | Reference allele.                                                                             |
| **ALT**           | Alternate allele.                                                                             |
| **P_FISHER**      | P-value calculated using Fisher's exact test. (binary analysis)                               |
| **P_CHI2**        | P-value calculated using the Chi-squared test. (binary analysis)                              |
| **P_value**       | P-value calculated using the linear regression. (quantitatif analysis)                        |
| **TOTAL_SUM**     | Total sum of all values in the DataFrame.                                                     |
| **MIN_ROW_INDEX** | Minimum sum among the row-wise sums of the DataFrame.                                         |
| **NUM_COLUM**     | Number of columns in the DataFrame.                                                           |
| **INTER_GROUP**   | Sum of the minimum values for each column.                                                    |
| **AVERAGE**       | Average of the total sum divided by the number of columns.                                    |

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
CHR POS Snarl           TYPE    REF ALT	P_value
1   12  5262721_5262719	SNP     A   T   0.9099
1   15  5262719_5262717	INS     A   ATT 0.9093
1   18  5262717_5262714	DEL     AA  T   0.9008
```

### Visualization

### Manhattan and QQ plots 

STOAT will generated a manhattan and a QQ plot for binary and quantitatif analysis.

### SequenceTube

Use `src/gaf_creator.py` to geneate a GAF file and [SequenceTube](https://github.com/vgteam/sequenceTubeMap) tool to visualize your gwas binary region results.

```bash 
python3 src/gaf_creator.py -s <binary_gwas_stoat_output.tsv> -l <decomposition_paths.tsv> -p <pg_file_path>
```

*Output Plots and sequenceTube*

