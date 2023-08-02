Tom Smith
02 August 2023

# Obtaining SMILES for PubChem compounds

From [this post](thttps://chemistry.stackexchange.com/questions/122109/how-to-get-the-smiles-of-all-compounds-on-pubchem):
_'For example, if you want the unfiltered SMILES of every CID in PubChem, the URL is ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-SMILES.gz'_

Following commands were used to create a file containing the SMILES for a random set of 25K PubChem compounds. Note that this does not ensure that the randomly selected compounds do not have the same SMILES or are truely representative of the complete chemical space of PubChem compounds

```bash
wget ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-SMILES.gz;
zcat < CID-SMILES.gz| shuf -n25000 --random-source=<(yes 1986) > CID-SMILES_random_25K.tsv;
rm -f CID-SMILES.gz
```