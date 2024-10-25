# snarl_project

## Project Overview
Snarl_project is a specialized tool developed for conducting Genome-Wide Association Studies (GWAS) with a unique focus on snarl structures within pangenome graphs. Unlike traditional GWAS tools that analyze linear genome variants, Snarl_project processes VCF files to extract and analyze snarl regions—complex structural variations that capture nested and overlapping variant patterns within a pangenome. This approach allows for a more nuanced understanding of genetic variations in diverse populations and complex traits.

Snarl_project supports both binary and quantitative phenotypes:

- For binary phenotypes (e.g., case vs. control studies), it utilizes chi-squared tests and Fisher’s exact test to evaluate associations between phenotype groups and snarl variants, providing robust statistical validation even in cases of sparse data.

- For quantitative phenotypes (e.g., traits measured on a continuous scale), the tool employs linear regression models to assess the association between snarl structures and phenotype values, allowing for continuous trait mapping with greater precision.

Conventional GWAS tools typically rely on a single-reference genome, often overlooking structural variants and variations present in underrepresented populations. By using a pangenome graph, Snarl_project enables analysis across multiple genomes simultaneously, capturing a broader spectrum of structural variations, including insertions, deletions, and nested variants. The focus on snarls—graph features that encapsulate complex variant structures—provides a powerful means to map associations that would be difficult to detect in linear-based analyses.

## Installation

````bash
git clone https://github.com/Plogeur/snarl_project.git
cd snarl_project
pip install -r requirements.txt
````

## Dependencies
- Python version 3.10
- cyvcf2
- numpy
- pandas
- ...

## Usage
Run the script from the command line, providing the paths to your VCF file, group file, and snarl file :

- Run full tool :
```bash
# binary trait
python3 snarl_project.py -p <path_to_pg_file.pg> -d <path_to_dist_file.dist> -v <path_to_vcf_file.vcf.gz> -r <path_to_vcf_reference_file.vcf.gz> -b <path_to_group_file.txt> -o output.tsv
```

```bash
# quantitative trait 
python3 snarl_project.py -p <path_to_pg_file.pg> -d <path_to_dist_file.dist> -v <path_to_vcf_file.vcf.gz> -r <path_to_vcf_reference_file.vcf.gz> -q <path_to_pheno_file.txt> -o output.tsv
```

- Specific script of tool : 
```bash
# decompose pangenome
python3 list_snarl_paths.py -p <path_to_pg_file.pg> -d <path_to_dist_file.dist> -o <output.tsv>

# binary trait
python3 snarl_vcf_parser.py <path_to_vcf_file.vcf.gz> <path_to_snarl_file.txt> <path_to_vcf_reference_file.vcf.gz> -b <path_to_group_file.txt> -o output.txt

# quantitative trait 
python3 snarl_vcf_parser.py <path_to_vcf_file.vcf.gz> <path_to_snarl_file.txt> <path_to_vcf_reference_file.vcf.gz> -q <path_to_pheno_file.txt> -o output.txt
```

### Input format file
Obligatory file :
- pg_file : Pangenome graph file, formats accepted: .pg or .xg.
- dist_file : Distance file generated with vg dist, format: .dist.
- vcf_ref : Reference node positions in the pangenome graph generated with vg deconstruct, formats: .vcf or .vcf.gz.
- vcf : Merged VCF file, created using bcftools merge, formats: .vcf or .vcf.gz.
- group/pheno : 
    - Binary phenotype: Two-column file with SAMPLE and GROUP labels. Format: .txt or .tsv (tab-separated).

    - Quantitative phenotype: Three-column file with FID (family/sample name), IID (sample name), and PHENO (integer/float). Format: .txt or .tsv (tab-separated).

Optional file : 
- snarl : Two-column file containing snarl names and their associated snarl list, separated by tabs. Format: .txt or .tsv.

### Output
Example of binary phenotype analysis (-b option) :
chrom, pos, snarl, type_var, ref, alt, fisher_p_value, chi2_p_value, total_sum, numb_colum, inter_group, average

```bash
CHR POS Snarl   TYPE    REF ALT P_value_Fisher  P_value_Chi2	Table_sum	Inter_group	Average
1   12  5262721_5262719	SNP A   T   0.46359729745943223	0.518292726549784	286	2	137	143.0
1   15  5262719_5262717	INS A   ATT 0.8062220214636773	0.8747410243373839	286	2	141	143.0
1   18  5262717_5262714	DEL AA  T   0.2120778233741457	0.2363840346684607	286	2	134	143.0
```

Example of quantitative phenotype analysis (-q option) :
```bash
CHR POS Snarl   TYPE    REF ALT	P_value
1   12  >5262719>5262721	SNP A   T   0.9099354411737626
1   15  >5262719>5262720	INS A   ATT 0.9099354411737626
1   18  >5262717>5262719	DEL AA  T   0.9008729687523177
```
