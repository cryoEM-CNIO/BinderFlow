'''
This is where all the functions of the watcher are going to be stored in order to make the file more readable

'''
import pandas as pd
import os
import re
import time 
import numpy as np
import glob 
import argparse
import subprocess
import plotly.express as px
import plotly.graph_objects as go
from Bio import SeqIO

#Count the residues in chain in order to get the length
#Oye, muy bien pensado al que haya escrito esto, chapeau
def count_residues_in_chain(pdb_file, chain_id='A'):
    '''
    Function to compute the lenght of the binders. Since the backbone is composed of Gly, the number of residues is computed as the number of atoms divided by 4 (the number of atoms in a _Gly residue) (chapeau Nayim cause this is faster than the old function and very witty).

    Input:

    pdb_file --> Path to the pdb file whose length is being computed 

    chain_id --> Chain ID whose length we want to compute, always pointing at A which is the binder chain in RFD

    Output:

    residue_count/4 --> The number of residues of the binder  
    '''
    residue_count = 0
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM') and line.split()[4] == chain_id:
                residue_count += 1
    return residue_count/4


def trim_substituted(text):
    '''
    Trim descriptions of silents to get original pdb designs names. This code is used to get the original design label. This is specially useful for the ngl visualization

    Input:

    text --> The description in the merged_df

    Output:

    whole_path --> Path to the original pdb in the output
    '''
    file = text.split('_substituted')[0]+'.pdb'
    whole_path = 'output/'+file.split('_design')[0]+'/'+file
    return whole_path

def merge_csv_files(directory, input_pdb_path):
    '''
    Function to merge all the scoring data (the one from AF2IG and from scoring.py into one df)

    Input:

    directory --> Directory where the data is stored (before the output)
    
    input_pdb_path --> Path to the input PDB is stored (For comparisons in scaffolding and partial diffusion cases, important because )
    
    Output:

    merged_df --> A df that integrates all the metrics from AF2IG and PyRosetta

    '''
    try:
        input_pdb_path=f'{directory}/{input_pdb_path}'
        df_whole=pd.read_csv(f'{directory}/Scoring_Stats.csv')
        df_whole['description']=df_whole['description'].astype(str)  # Ensure 'description' is a string
        df_whole['original_design'] = [trim_substituted(name) for name in df_whole['description'] if type(name) == str]
        input_df=pd.DataFrame(get_input_data(input_pdb_path))
        merged_df=pd.concat([df_whole,input_df], ignore_index=True)
        return merged_df
    except FileNotFoundError:
         return pd.DataFrame(columns=['plddt_binder','pae_interaction', 'CUTRE','dG','dSASA', 'Shape_complementarity', 'Packstat', 'dG_SASA_ratio', 'SAP', 'binder_int_hyd', 'binder_surf_hyd', 'interface_hbonds', 'interface_unsat_hbonds' ])

def check_logs_ended(log_file):
    '''
    Function to check if there is a done in the log file of the job. This is used to check if the job has finished or not
    
    Input:
    file --> The path to the log file of the job
    Output:
    Returns 1 if the job has finished, 0 if it has not finished
    '''
    with open(log_file, 'r') as file:
        lines = file.readlines()
        for line in lines:
            if 'done' in line.lower():
                return 1
            else:
                return 0


# Function to track job status
def track_job_status(directory):
    """
    Tracks job statuses for runs in the given directory.
    Returns a DataFrame with job, GPU, status, and path.
    """
    job_status = []
    states = ['Waiting', 'RFD', 'Filtering', 'pMPNN', 'AF2', 'Finished', 'FAILED']

    gpu_pattern = re.compile(r'.*\d+_gpu(\d+)')
    num_pattern = re.compile(r'\d+')

    for root, dirs, _ in os.walk(directory):
        for dir_name in dirs:
            job_path = os.path.join(root, dir_name)
            slurm_logs_path = os.path.join(job_path, "slurm_logs")
            if not os.path.exists(slurm_logs_path):
                continue

            for subdir in os.listdir(slurm_logs_path):
                subdir_path = os.path.join(slurm_logs_path, subdir)
                if not os.path.isdir(subdir_path):
                    continue

                gpu_match = gpu_pattern.match(subdir)
                if not gpu_match:
                    continue
                gpu_number = gpu_match.group(1)

                # Count completed logs
                counter = sum(check_logs_ended(os.path.join(subdir_path, f))
                              for f in os.listdir(subdir_path) if f.endswith("out"))

                job_record = {
                    "job": dir_name,
                    "gpu": gpu_number,
                    "status": states[counter],
                    "path": subdir_path
                }

                # Check done/sc files
                files_in_job = set(os.listdir(job_path))
                if f"{dir_name}_done" in files_in_job:
                    job_record["status"] = "Finished" if any(f.endswith(".sc") for f in files_in_job) else "FAILED"

                job_status.append(job_record)

    # Sort jobs numerically by extracted number
    status_df = pd.DataFrame(job_status)
    if not status_df.empty:
        status_df = status_df.sort_values(
            by="job",
            key=lambda col: col.map(lambda x: int(num_pattern.search(x).group()))
        )
    # Remove duplicates
    status_df = status_df.drop_duplicates()
    return status_df


#Function to get the data from the input
def get_input_data(input_pdb_path):
    '''
    Reads the last lines of the pdb used as input and stores the information in a df. The input must be the output from the watcher (to have all the metrics)
    or the output of a previous RFD run (only has the "soft filter" metrics)
    This is meant to be used in the case of scaffold or partial diffusion

    Input:

    input_pdb_path --> The path of where the input is stored

    Output:

    binding_analysis_dict --> A dictionary containing all the metrics of the initial input
    '''

    #Get the file working with regular expression using glob
    try:
        input_pdb_path=glob.glob(input_pdb_path)[0]
    except IndexError:
        input_pdb_path=input_pdb_path
    # Load dictionary
    binding_analysis_dict = {
        'pae_interaction':[],
        'pae_binder':[],
        'pae_target':[],
        'plddt_total':[],
        'plddt_binder':[],
        'plddt_target':[],
        'description': [],
        'original_design': [],
        'CUTRE': [],
        'dG': [],
        'dSASA': [],
        'Shape_complementarity': [],
        'Packstat': [],
        'dG_SASA_ratio': [],
        'SAP': [],
        'binder_int_hyd': [],
        'binder_surf_hyd': [],
        'interface_hbonds': [],
        'interface_unsat_hbonds': []
    }
    try:
        with open(f'{input_pdb_path}', 'r') as file:
            lines=file.readlines()
            for key in binding_analysis_dict.keys():
                for line in lines:
                # Iterate over the dictionary keys
                    if line.startswith(key):
                        # Extract the value after the key and any whitespace
                        value = line[len(key):].strip()
                        # Append the value to the corresponding key's list
                        binding_analysis_dict[key].append(float(value))
                        break
        #Assuming the input has the same name pattern as the one we use in microrun:
        binding_analysis_dict['original_design'].append(input_pdb_path)
        binding_analysis_dict['description'].append('Input')
        #check if some of the keys of the dictionary are empty and add a nan
        for key, value in binding_analysis_dict.items():
            if isinstance(value, list) and not value:  # Check if it's an empty list
                value.append(np.nan)
        print("Data loaded successfully.")
    except FileNotFoundError:
        print(f"File {input_pdb_path} not found.")
        for key in binding_analysis_dict.keys():
            binding_analysis_dict[key].append(np.nan)

    return binding_analysis_dict


#Filtering of DataFrame to get only hits according to the metrics
def filtering_df(merged_df,pae_interaction_thres, CUTRE_thres, plddt_binder_thres, dsasa_thres, shape_complementarity_thres, interface_hbond_thres, interface_unsat_hbond_thres, binder_surf_hyd_thres):
    '''
    --> Function to filter the df for the selection between the threshold. Separated to make code more legible

    Input:
    
    merged_df --> The dataframe with all the data
    All the thresholds --> (list of two elements containing the maximum and the minimum thresholds) 

    Output:

    filtered_df --> A df containing only the designs that fulfill all the metrics 
    
    '''
    # Filter the DataFrame based on the thresholds
    filtered_df=merged_df[merged_df['pae_interaction'].between(pae_interaction_thres[0], pae_interaction_thres[1])]
    filtered_df=filtered_df[filtered_df['CUTRE'].between(CUTRE_thres[0], CUTRE_thres[1])]
    filtered_df=filtered_df[filtered_df['plddt_binder'].between(plddt_binder_thres[0], plddt_binder_thres[1])]
    filtered_df=filtered_df[filtered_df['Shape_complementarity'].between(shape_complementarity_thres[0], shape_complementarity_thres[1])]
    filtered_df=filtered_df[filtered_df['dSASA'].between(dsasa_thres[0], dsasa_thres[1])]
    filtered_df=filtered_df[filtered_df['interface_hbonds'].between(interface_hbond_thres[0], interface_hbond_thres[1])]
    filtered_df=filtered_df[filtered_df['interface_unsat_hbonds'].between(interface_unsat_hbond_thres[0], interface_unsat_hbond_thres[1])]
    filtered_df=filtered_df[filtered_df['binder_surf_hyd'].between(binder_surf_hyd_thres[0], binder_surf_hyd_thres[1])]

    # Ensure unique values in the description column
    filtered_df = filtered_df.drop_duplicates(subset='description')
    return filtered_df



def get_working_directories(working_dir):
    """
    Get a list of directories within the parent directory that contain a folder named 'output'.
    
    Input:
    parent_dir (str): Path to the parent directory.
    
    Output:
    list: A list of paths to directories containing a folder named 'output'.
    """
    directories_list = []

    for root, dirs, files in os.walk(working_dir):
        if 'Scoring_Stats.csv' in files:  # Check if 'Scoring_Stats.csv' is one of the directories, since this file is needed for the plotting
            directories_list.append(root)

    return directories_list












