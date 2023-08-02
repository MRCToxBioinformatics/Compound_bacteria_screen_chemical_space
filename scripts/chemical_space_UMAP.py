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

import copy

import numpy as np

#from rdkit import DataStructs

from scipy.spatial.distance import pdist, jaccard, squareform
from collections import Counter, defaultdict

# Functions to use UMAP to embed and project and then plot
def get_UMAP_vs_background(plotter_obj,
                           background_name='PubChem',
                           foreground_name=None,
                           n_neighbors = 20,
                           min_dist = 0.25,
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
                            palette=('#EEEEEE',
                                     "#F0E442",
                                     "#56B4E9",
                                     "#D55E00",
                                     "#E69F00",
                                     "#0072B2",
                                     "#000000",
                                     "#CC79A7",
                                     "#009E73"),
                            outfile_prefix = None):
	
	size=20
	hue = 'type'

	palette = palette
	sns.set_style("white")
	sns.set_context("notebook", font_scale=size*0.15)
	fig, ax = plt.subplots(figsize=(size,size))

	plot = sns.scatterplot(x=0, y=1, hue=hue,
                               hue_order=['PubChem',
                                          'Pharmaceutical drugs',
                                          'Pesticide',
                                          'Pesticide metabolite',
                                          'Pesticide-related',
                                          'Industrial chemical',
                                          'Mycotoxin'],
                               palette=palette, data=df_data, s=size*3)

	plot.set_label("scatter")
	axis = plot
	plot.legend(markerscale=size*0.1, fontsize=size*1.2, frameon=False)
	# Remove units from axis
	axis.set(yticks=[])
	axis.set(xticks=[])
	# Add labels
	axis.set_title('',fontsize=size*2)
	axis.set_xlabel('UMAP1',fontsize=size*2)
	axis.set_ylabel('UMAP2',fontsize=size*2)
	axis.spines['top'].set_visible(False)
	axis.spines['right'].set_visible(False)
	axis.spines['bottom'].set_visible(False)
	axis.spines['left'].set_visible(False)

	if(outfile_prefix):
		plt.savefig(outfile_prefix + '.png')
		plt.savefig(outfile_prefix + '.pdf')


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


	parser.set_defaults(n_neighbors=20,
	                    min_dist=0.1,
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

	all_cids = list(pubchem_smiles['CID']) + list(screen_smiles['PubChem CID'])


	all_targets_random = copy.copy(all_targets)

	np.random.seed(args.seed)
	random_selection = np.random.choice(range(0, len(pubchem_smiles['smile'])),
	                                    size=len(screen_smiles["Compound Class"]), replace=False)

	for ix in random_selection:
		all_targets_random[ix] = 'Random'
		
	logging.info("Obtaining the ECFPs")

	cp = Plotter.from_smiles(all_smiles, target=all_targets, target_type="C", sim_type='structural')

	logging.info("UMAP")

	umap = get_UMAP_vs_background(cp, n_neighbors=args.n_neighbors, min_dist=args.min_dist)

	logging.info("Plotting")

	plot_UMAP_vs_background(umap, outfile_prefix=args.output_prefix)

	logging.info("Finished!")

if __name__ == "__main__":
	sys.exit(main())
