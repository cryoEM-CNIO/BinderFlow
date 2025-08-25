#!/usr/bin/env python3 


import pandas as pd 
import argparse
import subprocess

'''
Read the csv with the scorings and the sc files, get all the info, and end the run once enough binders are generated

Input:

--number: Max number of hits that want to be generated

Output:

campaign_done: An empty file that marks that all wanted binders have been generated and stops the run
'''

# Define number of designs desired
parser=argparse.ArgumentParser()
parser.add_argument('--number', '-n' , type=int, default=100, help='Maximum number of successful binder desired')
args=parser.parse_args()

number=args.number



# define Thresholds

pae_interaction_thres=10
plddt_binder_thres=80
# CUTRE_thres=10
# unsat_hbonds_thres=4
# hbond_thres=3
# binder_surface_hyd_thres=0.35
# shape_complementarity_thres=0.55
# dSASA_thres=1000 # is a percentage better ?

# Read the csv file
try:
	df=pd.read_csv('Scoring_Stats.csv')

	# Load the variable 

	hits_number=0

	# Check the conditions

	hits_number = (
	# 	(df['CUTRE'] <= CUTRE_thres) &
	# 	(df['interface_unsat_hbonds'] <= unsat_hbonds_thres) &
	# 	(df['interface_hbonds'] >= hbond_thres) &
	# 	(df['binder_surf_hyd'] <= binder_surface_hyd_thres) &
	# 	(df['Shape_complementarity'] >= shape_complementarity_thres) &
	# 	(df['dSASA'] >= dSASA_thres)
		(df['pae_interaction'] <= pae_interaction_thres) & 
		(df['plddt_binder'] >= plddt_binder_thres)).sum()
	
	command='touch campaign_done'

	if hits_number >= number:
		print('done')
		subprocess.run(command, shell=True)
	else:
		print('continue') 

except FileNotFoundError:
	print('No Run has been completed yet')
