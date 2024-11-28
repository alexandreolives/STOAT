# snarl_project
### Project Overview
Snarl_project is a specialized tool developed for conducting Genome-Wide Association Studies (GWAS) with a unique focus on snarl structures within pangenome graphs. Unlike traditional GWAS tools that analyze linear genome variants, Snarl_project processes VCF files to extract and analyze snarl regions—complex structural variations that capture nested and overlapping variant patterns within a pangenome. This approach allows for a more nuanced understanding of genetic variations in diverse populations and complex traits.

Snarl_project supports both binary and quantitative phenotypes:

- For binary phenotypes (e.g., case vs. control studies), it utilizes chi-squared tests and Fisher’s exact test to evaluate associations between phenotype groups and snarl variants, providing robust statistical validation even in cases of sparse data.

- For quantitative phenotypes (e.g., traits measured on a continuous scale), the tool employs linear regression models to assess the association between snarl structures and phenotype values, allowing for continuous trait mapping with greater precision.

Conventional GWAS tools typically rely on a single-reference genome, often overlooking structural variants and variations present in underrepresented populations. By using a pangenome graph, Snarl_project enables analysis across multiple genomes simultaneously, capturing a broader spectrum of structural variations, including insertions, deletions, and nested variants. The focus on snarls—graph features that encapsulate complex variant structures—provides a powerful means to map associations that would be difficult to detect in linear-based analyses.

### Installation

````bash
git clone https://github.com/Plogeur/snarl_project.git
cd snarl_project
pip install -r requirements.txt
````
### snarl_project

### Dependencies
- Python version 3.10
- cyvcf2
- numpy
- pandas
- ...
more information available on the [README](https://github.com/Plogeur/snarl_project/blob/main/README.md)

### Usage

Use `snarl_project.py` if you want to launch the full tool at once, starting from snarl path identification (identifying the multiple paths that can be taken by a sample based on the pangenome graph) and ending with the results plots (Manhattan plot and QQ plot).

- Run full tool :
```bash
# binary trait
python3 snarl_project.py -p <path_to_pg_file.pg> -d <path_to_dist_file.dist> -v <path_to_vcf_file.vcf.gz> -r <path_to_vcf_reference_file.vcf.gz> -b <path_to_group_file.txt> -o output.tsv
python3 snarl_project.py -p <path_to_pg_file.pg> -d <path_to_dist_file.dist> -v <path_to_vcf_file.vcf.gz> -r <path_to_vcf_reference_file.vcf.gz> -q <path_to_pheno_file.txt> -o output.tsv
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

### Input format file

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

### Output
Example of binary phenotype analysis (-b option) :

```bash
CHR POS Snarl   TYPE    REF ALT P_value_Fisher  P_value_Chi2	Table_sum	Inter_group	Average
1   12  5262721_5262719	SNP A   T   0.46359729745943223	0.518292726549784	286	2	137	143.0
1   15  5262719_5262717	INS A   ATT 0.8062220214636773	0.8747410243373839	286	2	141	143.0
1   18  5262717_5262714	DEL AA  T   0.2120778233741457	0.2363840346684607	286	2	134	143.0
```

Example of quantitative phenotype analysis (-q option) :

```bash
CHR POS Snarl   TYPE    REF ALT	P_value
1   12  5262721_5262719	SNP A   T   0.9099354411737626
1   15  5262719_5262717	INS A   ATT 0.9099354411737626
1   18  5262717_5262714	DEL AA  T   0.9008729687523177
```

*Plots and sequenceTube output*
coming soon

