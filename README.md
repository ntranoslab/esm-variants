# esm-variants

## Genome-wide prediction of disease variants with a deep protein language model

This repository contains resources and tools developed for the paper, ["Genome-wide prediction of disease variants with a deep protein language model"](https://www.biorxiv.org/content/10.1101/2022.08.25.505311v1) by Nadav Brandes, Grant Goldman, Charlotte H. Wang, Chun Jimmie Ye, and Vasilis Ntranos.

A complete catalog of missense variant effect predictions is accessible [here](https://huggingface.co/spaces/ntranoslab/esm_variants).

## Repository Contents

### Notebooks/Scripts

- `esm_variants_utils.ipynb`
- `esm_variants_utils.py`

### Benchmarks
- `ClinVar_gnomAD_benchmark_with_predictions.csv`
- `ClinVar_indel_benchmark_with_predictions.csv`
- `ClinVar_stop_gains_benchmark_with_predictions.csv.gz`
- `dms_assays.zip`

These files contain all benchmark data and VEP predictions used for performance evaluation, except HGMD variants (see below).

### Results

- `Table_of_results.xlsx`

## Data Availability

Most data used in this work is already within the public domain. Exceptions and other data sources are detailed in the paper and below:

- The HGMD dataset is a private resource owned by the Institute of Medical Genetics in Cardiff University. [Request access here](https://www.hgmd.cf.ac.uk/ac/index.php).
- ClinVar labels of missense variants in all protein isoforms and of indels and stop-gains are available [here](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/).
- Full details on how the datasets and benchmarks were processed are available in the [Supplementary Methods](https://www.biorxiv.org/content/10.1101/2022.08.25.505311v1).
- Predicted effect scores for most VEP methods were downloaded from [dbNSFP](http://database.liulab.science/dbNSFP).

## Related Resources

- ESM1b: This study leverages and expands the use of ESM1b, a protein language model developed by MetaAI. The code and pre-trained parameters for ESM1b were taken from the model’s [official GitHub repository](https://github.com/facebookresearch/esm).
- Web portal: We created a [web portal](https://huggingface.co/spaces/ntranoslab/esm_variants) allowing researchers to query, visualize, and download missense variant effect predictions for all protein isoforms in the human genome.
