# esm-variants & ESM Variants Command-Line Tools
This repository contains resources, tools, and command-line tools developed for the paper, ["Genome-wide prediction of disease variants with a deep protein language model"](https://www.biorxiv.org/content/10.1101/2022.08.25.505311v1) by Nadav Brandes, Grant Goldman, Charlotte H. Wang, Chun Jimmie Ye, and Vasilis Ntranos. A complete catalog of missense variant effect predictions is accessible [here](https://huggingface.co/spaces/ntranoslab/esm_variants).
## Repository Contents
### Notebooks/Scripts
- `esm_variants_utils.ipynb`
- `esm_variants_utils.py`
### Command-Line Tools
- `esm_score_missense_mutations.py`
- `esm_score_multi_residue_mutations.py`
### Benchmarks
- `ClinVar_gnomAD_benchmark_with_predictions.csv`
- `ClinVar_indel_benchmark_with_predictions.csv`
- `ClinVar_stop_gains_benchmark_with_predictions.csv.gz`
- `dms_assays.zip`
### Results
- `Table_of_results.xlsx`
These files contain all benchmark data, VEP predictions used for performance evaluation, and results, except HGMD variants (see below).
## Data Availability
Most data used in this work is already within the public domain. Exceptions and other data sources are detailed in the paper and below:
- The HGMD dataset is a private resource owned by the Institute of Medical Genetics in Cardiff University. [Request access here](https://www.hgmd.cf.ac.uk/ac/index.php).
- ClinVar labels of missense variants in all protein isoforms and of indels and stop-gains are available [here](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/).
- Full details on how the datasets and benchmarks were processed are available in the [Supplementary Methods](https://www.biorxiv.org/content/10.1101/2022.08.25.505311v1).
- Predicted effect scores for most VEP methods were downloaded from [dbNSFP](http://database.liulab.science/dbNSFP).
## Related Resources
- ESM1b: This study leverages and expands the use of ESM1b, a protein language model developed by MetaAI. The code and pre-trained parameters for ESM1b were taken from the model’s [official GitHub repository](https://github.com/facebookresearch/esm).
- Web portal: We created a [web portal](https://huggingface.co/spaces/ntranoslab/esm_variants) allowing researchers to query, visualize, and download missense variant effect predictions for all protein isoforms in the human genome.
## Installation
The following dependencies are required:
```
pip3 install tqdm numpy pandas biopython torch fair-esm
```
Clone the repository:
```
git clone https://github.com/ntranoslab/esm-variants.git
cd esm-variants
```
## Usage
### Scoring Missense Mutations
```
python3 esm_score_missense_mutations.py --input-fasta-file /path/to/input.fasta --output-csv-file /path/to/output.csv
```
### Scoring Multi-Residue Mutations
```
python3 esm_score_multi_residue_mutations.py --input-csv-file /path/to/input.csv --output-csv-file /path/to/output.csv
```
The input CSV file for multi-residue mutations should have three fields:
- `wt_seq`: the wild type (original) protein sequence
- `mut_seq`: the mutated protein sequence
- `start_pos`: the starting position (1-indexed) of the mutation relative to the wild type sequence
## Examples
### Example for Scoring Missense Mutations
Assuming an example FASTA file named `example.fasta`:
```
>seq1
FISHWISHFQRCHIPSTHATARECRISP
>seq2
RAGEAGAINSTTHEMACHINE
```
You can calculate ESM scores for all possible missense mutations in these sequences:
```
python3 esm_score_missense_mutations.py --input-fasta-file example.fasta --output-csv-file esm_scores.csv
```
This will create a CSV file (`esm_scores.csv`) that starts like this:
```
seq_id,mut_name,esm_score
seq1,F1K,-3.2310808
seq1,F1R,-2.872289
seq1,F1H,-3.4361703
...
```
Each row represents a possible missense mutation and its ESM score.
### Example for Scoring Multi-Residue Mutations
Assuming the following `example.csv`:
```
wt_seq,mut_seq,start_pos
FISHWISHFQRCHIPSTHATARECRISP,FISHWISHFQRCHEESETHATARECRISP,14
MARGTYNMGKHFDA,MGTYNMGKHFDA,2
```
You can calculate ESM (PLLR) scores for the specified multi-residue mutations:
```
python3 esm_score_multi_residue_mutations.py --input-csv-file example.csv --output-csv-file esm_multi_residue_scores.csv
```
This will create a CSV file (`esm_multi_residue_scores.csv`) that starts like this:
```
wt_seq,mut_seq,start_pos,esm_score
FISHWISHFQRCHIPSTHATARECRISP,FISHWISHFQRCHEESETHATARECRISP,14,-1.0078125
MARGTYNMGKHFDA,MGTYNMGKHFDA,2,1.0056415
```
Each row represents a multi-residue mutation and its ESM (PLLR) score.
