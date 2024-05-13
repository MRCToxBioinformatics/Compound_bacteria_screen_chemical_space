'''
================================
chemical_space.py
================================

This script takes two files of SMILES, representing compounds that have been screened
(foreground) and a random subset of PubChem compounds (background).

A UMAP manifold is obtained from the background, onto which the foreground is also
plotted.

The foreground is further expected to have a column describing the screen library

'''
import sys
import argparse
import logging

import umap

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from chemplot import Plotter

import numpy as np

from scipy.spatial.distance import cdist, jaccard

# Functions to use UMAP to embed and project and then plot
def get_UMAP_vs_background(plotter_obj,
                           background_name='PubChem',
                           foreground_name=None,
                           n_neighbors = 20,
                           min_dist = 0.2,
                           random_state=1985):
   
	data = plotter_obj._Plotter__df_descriptors
	data_bg = data.loc[[x==background_name for x in plotter_obj._Plotter__target],:]

	df_type = [background_name]*data_bg.shape[0]

	if foreground_name:
		data_fg = data.loc[[x==foreground_name for x in plotter_obj._Plotter__target],:]
		df_type.extend([foreground_name]*data_fg.shape[0])
	else:
		data_fg = data.loc[[x!=background_name for x in plotter_obj._Plotter__target],:]
		df_type.extend([x for x in plotter_obj._Plotter__target if x!=background_name])

	umap_fit = umap.UMAP(n_neighbors=n_neighbors, min_dist=min_dist, random_state=random_state,
						 n_components=2, metric='jaccard')
		
	trans = umap_fit.fit(data_bg)

	bg_proj = pd.DataFrame(trans.embedding_)
	fg_proj = pd.DataFrame(trans.transform(data_fg))

	df_data = pd.concat([bg_proj, fg_proj])
	df_data['type'] = df_type

	return(df_data)

def plot_UMAP_vs_background(df_data,
                            palette=('#E2E2E2',
                                     "#F0E442",
                                     "#56B4E9",
                                     "#D55E00",
                                     "#E69F00",
                                     "#0072B2",
                                     "#000000",
                                     "#CC79A7",
                                     "#009E73"),
                            outfile_prefix = None,
                            highlight_compounds=None):
	
	size=20/2.54
	hue = 'type'

	palette = palette
	sns.set_style("white")
	sns.set_context("notebook", font_scale=size*0.15)
	fig, ax = plt.subplots(figsize=(size,size/2))

	plot = sns.scatterplot(x=0, y=1, hue=hue,
                               hue_order=['PubChem',
                                          'Pharmaceutical drugs',
                                          'Pesticide',
                                          'Pesticide metabolite',
                                          'Pesticide-related',
                                          'Industrial chemical',
                                          'Mycotoxin'],
                               palette=palette, data=df_data,
                               s=size*0.5,
                               rasterized=True)

	if highlight_compounds:
		highlight_df = df_data[df_data['CID'].isin(highlight_compounds)]
		plt.scatter(highlight_df[0], highlight_df[1], marker='x',
		            color='black', s=size*0.5, linewidths=0.5)

	plot.set_label("scatter")
	axis = plot
	plot.legend(markerscale=size*0.5, fontsize=size, frameon=False)
	# Remove units from axis
	axis.set(yticks=[])
	axis.set(xticks=[])
	# Add labels
	axis.set_title('')
	axis.set_xlabel('UMAP1',fontsize=size)
	axis.set_ylabel('UMAP2',fontsize=size)
	axis.spines['top'].set_visible(False)
	axis.spines['right'].set_visible(False)
	axis.spines['bottom'].set_visible(False)
	axis.spines['left'].set_visible(False)

	plt.rcParams['pdf.fonttype'] = 42
	plt.rcParams['svg.fonttype'] = 'none'

	if(outfile_prefix):
		plt.savefig(outfile_prefix + '.png', dpi=400)
		plt.savefig(outfile_prefix + '.pdf', dpi=400)

		if highlight_compounds:
			highlight_df.to_csv(outfile_prefix + '_highlight_CIDs.tsv', index=False, sep='\t')


def main(argv=None):
	"""script main.

	"""
	
	plt.ioff()

	if argv is None:
		argv = sys.argv

	parser = argparse.ArgumentParser()

	parser.add_argument('--pubchem-smiles-infile', '-pc',
						dest='pubchem_smiles_inf',
						type=str,
						required=True,
						help='Filename for the pubchem SMILES')

	parser.add_argument('--screen-compound-smiles-infile', '-sc',
						dest='screen_compound_smiles_inf',
						type=str,
						required=True,
						help='Filename for the screen compound SMILES')


	parser.add_argument('--seed',
						dest='seed',
						type=int,
						help='Value for seed used in K-fold stratification')


	parser.add_argument('--output-prefix', '-out',
						dest='output_prefix',
						type=str,
						help='Prefix for the outfiles')

	parser.add_argument('--logfile', '-log',
						dest='logname',
						type=str,
						help=('Filname for log file'))


	parser.add_argument('--n-neighbors', '-nn',
						dest='n_neighbors',
						type=int,
						help=('UMAP Number of neighbours'))


	parser.add_argument('--minimimum-distance', '-md',
						dest='min_dist',
						type=float,
						help=('UMAP Mimumum distance'))

	parser.add_argument('--highlight-compounds',
					dest='highlight_compounds',
					type=float,
					nargs='+',
					default=None,
					help=('CIDs for compounds to highlight'))


	parser.set_defaults(n_neighbors=20,
						min_dist=0.2,
						seed=1985,
						output_prefix='./umap')

	args = parser.parse_args()

	if not args.logname:
                logging.basicConfig(stream=sys.stdout,
                                    format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                                    datefmt='%H:%M:%S')
	else:
		logging.basicConfig(filename=args.logname,
							filemode='w',
							format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
							datefmt='%H:%M:%S',
							level=logging.INFO)

	logging.getLogger('matplotlib.font_manager').disabled = True


	logging.info("Script parameters:")
	for k, v in vars(args).items():
		logging.info("%s: %s" % (k, v))
	logging.info("----")  

	logging.info("Reading in the SMILES")

	pubchem_smiles = pd.read_table(args.pubchem_smiles_inf, header=None)
	pubchem_smiles.columns = ['CID', 'smile']

	# file contains smiles for all compounds
	screen_smiles = pd.read_table(args.screen_compound_smiles_inf, index_col=0)

	screen_smiles['Compound Class'] = screen_smiles['Compound Class'].str.capitalize()

	screen_smiles.dropna(inplace=True)


	all_smiles = list(pubchem_smiles['smile']) + list(screen_smiles["PUBCHEM_isosmiles"])
	all_targets = ['PubChem']*len(pubchem_smiles['smile'])+ list(screen_smiles["Compound Class"])
	
	logging.info("Obtaining the ECFPs")

	cp = Plotter.from_smiles(all_smiles, target=all_targets, target_type="C", sim_type='structural')

	logging.info("UMAP")

	my_umap = get_UMAP_vs_background(cp,
	                                 n_neighbors=args.n_neighbors,
	                                 min_dist=args.min_dist,
	                                 random_state=args.seed)
	
	my_umap = my_umap.reset_index()

	if args.highlight_compounds is not None:
		# Add CID if highlighting of specific CID is required
		retained_smiles = set(cp._Plotter__df_descriptors.index)
		cids = list(pubchem_smiles['CID']) + list(screen_smiles['PubChem CID'])
		cids = [cid for cid, smiles in zip(cids, all_smiles) if smiles in retained_smiles]
		my_umap['CID'] = cids

	logging.info("Plotting")

	plot_UMAP_vs_background(my_umap, outfile_prefix=args.output_prefix,
	                        highlight_compounds=args.highlight_compounds)

	logging.info("Extracting ECFPs")
	all_compounds2smiles = pd.DataFrame(all_targets, all_smiles, columns=['Compound Class'])

	ecfp = cp._Plotter__df_descriptors
	ecfp = ecfp.merge(all_compounds2smiles, left_index=True, right_index=True)

	logging.info("Convert to Jaccard distance")	
	
	screen_rows = ecfp[ecfp['Compound Class']!='PubChem']
	pubchem_rows = ecfp[ecfp['Compound Class']=='PubChem']
	
	jacc_dis = cdist(pubchem_rows.drop('Compound Class', axis=1),
	                 screen_rows.drop('Compound Class', axis=1),
	                 'jaccard')

	logging.info("Get minimum distances between compounds")	
	
	pollutant_classes = ['Pesticide', 'Pesticide metabolite', 'Pesticide-related', 'Industrial chemical', 'Mycotoxin']
	pollutant_indexes = [ix for ix,x in enumerate(screen_rows['Compound Class']) if x  in pollutant_classes]
	
	drug_indexes = [ix for ix,x in enumerate(screen_rows['Compound Class']) if x =='Pharmaceutical drugs']

	min_jacc = np.append(
		np.append(np.min(jacc_dis, 1),
		          np.min(jacc_dis[:,pollutant_indexes], 1)),
		np.min(jacc_dis[:,drug_indexes], 1))
	
	comparison = ((['All',]*len(pubchem_rows.index)) +
	              ['Environmental Pollutants',]*len(pubchem_rows.index) + 
	              ['Pharmaceutical drugs',]*len(pubchem_rows.index))

	min_jacc_df = pd.DataFrame({'min_jacc': min_jacc,
	                           'comparison': comparison})

	logging.info("Saving minimum jaccard distances to file")		
	min_jacc_df.to_csv(args.output_prefix + '_min_jaccard.tsv', index=False, sep='\t')
	
	logging.info("Get minimum distances between pollutants and drugs")	
	
	pollutant_rows = ecfp[ecfp['Compound Class'].isin(pollutant_classes)]
	drug_rows = ecfp[ecfp['Compound Class']=='Pharmaceutical drugs']

	pollutant_drug_jacc_dis = cdist(pollutant_rows.drop('Compound Class', axis=1),
	                                drug_rows.drop('Compound Class', axis=1),
	                                'jaccard')

	min_pollutant_drug_jacc = np.min(pollutant_drug_jacc_dis, 1)
	
	# Dicts required to map from smile back to name and CID
	smile2chem = {x:y for x,y in zip(screen_smiles['PUBCHEM_isosmiles'], screen_smiles['Chemical name'])}
	smile2cid = {x:y for x,y in zip(screen_smiles['PUBCHEM_isosmiles'], screen_smiles['PubChem CID'])}

	min_jac_pol_drug_rows = list()
	
	# Find the drug with the min jacc distance for each pollutant
	for pol_ix, min_jacc in enumerate(min_pollutant_drug_jacc):
	    pol_smile = pollutant_rows.index[pol_ix]
	    
	    # Will break ties since only first match returned
	    min_ix = list(pollutant_drug_jacc_dis[pol_ix]).index(min_jacc)

	    drug_smile = drug_rows.index[min_ix]
	    
	    min_jac_pol_drug_rows.append([smile2cid[pol_smile],
	          smile2chem[pol_smile].upper(),
	          smile2cid[drug_smile],
	          smile2chem[drug_smile].upper(),
	          1-min_jacc])

	min_pollutant_drug_jacc_df = pd.DataFrame.from_records(min_jac_pol_drug_rows)
	min_pollutant_drug_jacc_df.columns = [
		'Pollutant_PubChem_CID', 'Pollutant_name', 'Drug_PubChem_CID', 'Drug_name', 'Tanimoto']
	min_pollutant_drug_jacc_df.sort_values('Tanimoto', ascending=False, inplace=True)
	
	min_pollutant_drug_jacc_df.to_csv(args.output_prefix + '_min_pollutant_drug_jaccard.tsv', sep='\t')

	logging.info("Finished!")

if __name__ == "__main__":
	sys.exit(main())
