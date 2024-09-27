# Chemical space visualisations for bacterial growth screen compounds
This repository contains the scripts and output for the visualisations of the chemical space coverage for the compounds screened for their impact on bacterial growth. The results are presented in [Lindell, Kamrad &amp; Roux et al](https://www.biorxiv.org/content/10.1101/2024.09.05.610817v1).


## Set up environment
To create a suitable environment to run the script, we suggest using conda/mamba. The following commands will create a working environment with the dependencies installed.

```bash
conda create -n compound_visualisation python=3.10
conda activate compound_visualisation
mamba install -y -c conda-forge -c bioconda umap-learn pandas numpy scipy \
    matplotlib seaborn chemplot bokeh==2.4.3 \
    r-ggplot2 r-optparse r-dplyr r-devtools
R -e "devtools::install_github('VPetukhov/ggrastr')"
```

To generate the UMAP plot, run:

```bash
python scripts/chemical_space_UMAP.py -pc raw/CID-SMILES_random_25K.tsv -sc raw/compound2smile.tsv --log umap.log --output-prefix results/umap
```

To plot the Tanimoto distance between PubChem compounds and screened compounds, run the following. Note that `results/umap_min_jaccard.tsv` and `results/umap_min_pollutant_drug_jaccard.tsv` are files created by `scripts/chemical_space_UMAP.py`.
```bash
scripts/plot_tanimoto.R --min-jacc-inf results/umap_min_jaccard.tsv --min-jacc-pol-drug-inf results/umap_min_pollutant_drug_jaccard.tsv --plot-prefix results/plot 
```
