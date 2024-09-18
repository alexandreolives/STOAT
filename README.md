# snarl_project

## Project Overview
Snarl_project is tool designed to parse and analyze VCF (Variant Call Format) files, particularly focusing on snarl structures. The tool processes VCF files and extracts snarl information, in order to generated snarl's statistique on binary or quantitative phenotype.

## Installation

````bash
git clone https://github.com/yourusername/snarl-vcf-parser.git
cd snarl-vcf-parser
pip install -r requirements.txt
````

## Dependencies
- Python version 3.10
- cyvcf2
- numpy
- pandas

## Usage
Run the script from the command line, providing the paths to your VCF file, group file, and snarl file:

```bash
python3 snarl_vcf_parser.py <path_to_vcf_file.vcf.gz> <path_to_snarl_file.txt> -b <path_to_group_file.txt> -o output.tx
```

### Input format file
- vcf_path : Merged VCFs using `bcftools merge`, input format : .vcf or .vcf.gz.
- snarl : File that containt 2 column : snarl name and list of snarl (separeted by an tabulation), input format : .txt or .tsv.
- group : File that containt 2 column : Sample name for the groupe 1 and sample name for the groupe 2 (separeted by an tabulation), input format : .txt or .tsv.

### Run test
```bash
python3 src/snarl_vcf_parser.py test/small_vcf.vcf test/list_snarl_short.txt -b test/group.txt
```
