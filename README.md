# Chemical space visualisations for bacterial growth screen compounds
This repository contains the scripts and output for the visualisations of the chemical space coverage for the compounds screened for their impact on bacterial growth. The results are presented in Kamrad &amp; Lindell et al (manuscript link TBC).

To generate the UMAP plot, run:

```bash
python scripts/chemical_space_UMAP.py -pc raw/CID-SMILES_random_25K.tsv -sc raw/compound2smile.tsv --log umap.log --output-prefix results/plots/umap
```

To plot the Tanimoto distance between PubChem compounds and screened compounds, run:
```bash
(ADD COMMAND HERE)
```