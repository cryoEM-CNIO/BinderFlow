#!/usr/bin/env python3
''''
This code is intended to compute the Rosetta scoring statistics of the designs

Input:

--silent: silent file with the pdb structures to analyze (output of AF2)
--run_number: Run number

Output:

Scoring_Stats.csv: csv with different metrics to evaluate the designs

'''

from pyrosetta import *
import argparse
import pandas as pd
import json
import numpy as np
import glob
import re
import time
import os
from superimpose import superpose_pose_by_chain

def get_sc_scorings(run_number):
    '''Function to get AF2-IG prediction metrics of the designs'''
    working_directory=f'./output/run_{run_number}'
    df_list=[]
    for root, dirs, files in os.walk(working_directory):
        for file in files:
            if file.endswith('.sc'):
                file_path = os.path.join(root, file)
                df = pd.read_table(file_path, sep=r'\s+', encoding='utf-8')
                df_list.append(df)
    concat_df=pd.concat(df_list)          
    return concat_df

def get_close_residues(pose,distance=10 ):
    '''
    Function to get the interacting residfues between binder and target (10A as cutoff distance)
    The notation is a little bit messy: 
    In the first part we are computing the residues close to the binder (close_residues_target because they belong to the target)
    In the second part we are computing residues close to the target (close_residues_binder because they belong to the binder) 
    '''
    NRS_binder=pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector()
    NRS_target=pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector()
    NRS_target.set_distance(distance)
    NRS_binder.set_distance(distance)
    #Get the binder start and ending numeration
    start_residue_binder=pose.conformation().chain_begin(1)
    end_residue_binder=pose.conformation().chain_end(1)

    #Set focus and search neighbors of the target
    NRS_binder.set_focus(f'{start_residue_binder}-{end_residue_binder}')
    neighbors_binder=NRS_binder.apply(pose)

    #Save the neighbours of the target
    close_residues_target=[]
    for i in range(end_residue_binder,len(neighbors_binder)):
        if neighbors_binder[i]:
            close_residues_target.append(i)
    
    
    #Get the target numeration
    start_residue_target=pose.conformation().chain_begin(2)
    end_residue_target=pose.conformation().chain_end(2)

    #Set focus and search
    NRS_target.set_focus(f'{start_residue_target}-{end_residue_target}')
    neighbors_target=NRS_target.apply(pose)

    #Save the neighbors of the binder
    close_residues_binder=[]
    for j in range(start_residue_binder, end_residue_binder):
        if neighbors_target[j]:
            close_residues_binder.append(j)

    length=int(end_residue_binder)

    return close_residues_binder, close_residues_target, length 

def compute_CUTRE(run_number,protein_name, close_residues_binder, close_residues_target):
    '''Function to compute the CUTRE score'''
    #Path to the json file where all the pae info is stored 
    json_file_path=glob.glob(fr'output/run_{run_number}/*{protein_name}.json')[0]

    #Load data from JSON file
    with open(json_file_path, 'r') as json_file:
        data=json.load(json_file)

    pae=data['predicted_aligned_error']
    plddt=data['plddt']
    residues=[]
    CUTRE=[]

    average_number_interface_residues=(len(close_residues_binder)+len(close_residues_target))/2

    if average_number_interface_residues >= 7: 
        for residue_target in close_residues_target: #Please Carlos, I am you, review how the PAE matrix work and put it here so this code it is easier to understand 
            lists=pae[residue_target-1]
            for residue_binder in close_residues_binder:
                CUTRE.append(lists[residue_binder-1]/(plddt[residue_target-1]/100))
        for residue_binder in close_residues_binder:
            lists=pae[residue_binder-1]
            for residue_target in close_residues_target:
                CUTRE.append(lists[residue_target-1]/(plddt[residue_binder-1]/100))

        CUTRE=np.mean(CUTRE)

    else:
        CUTRE=np.nan

    return CUTRE

def compute_rosetta_metrics(pose,interface):
    '''Function to compute binder rosetta metrics'''

    # Rosetta Ops
    interface_analyzer=pyrosetta.rosetta.protocols.analysis.InterfaceAnalyzerMover()
    interface_analyzer.set_interface(interface)
    scorefxn=pyrosetta.get_fa_scorefxn()
    interface_analyzer.set_scorefunction(scorefxn)
    interface_analyzer.set_compute_packstat(True)
    interface_analyzer.set_compute_interface_energy(True)
    interface_analyzer.set_calc_dSASA(True)
    interface_analyzer.set_calc_hbond_sasaE(True)
    interface_analyzer.set_compute_interface_sc(True)
    interface_analyzer.set_pack_separated(True)
    interface_analyzer.apply(pose)
    

    ##Calculate things
    interface_score=interface_analyzer.get_all_data()
    interface_sc=interface_score.sc_value #Shape complementarity
    interface_dsasa=interface_analyzer.get_interface_delta_sasa() #Interface dSASA
    interface_packstat=interface_analyzer.get_interface_packstat()#Interface packstat score (above 0.65 considered good enough)
    interface_hbonds=interface_score.interface_hbonds
    buns_filter=pyrosetta.rosetta.protocols.rosetta_scripts.XmlObjects.static_get_filter('<BuriedUnsatHbonds report_all_heavy_atom_unsats="true" scorefxn="scorefxn" ignore_surface_res="false" use_ddG_style="true" dalphaball_sasa="0" probe_radius="1.1" burial_cutoff_apo="0.2" confidence="0" />')
    interface_delta_unsat_hbonds = buns_filter.report_sm(pose)

    #SAP calculus
    binder_chain_selector=pyrosetta.rosetta.core.select.residue_selector.ChainSelector(1) #Select binder chain

    SAPscm=pyrosetta.rosetta.core.pack.guidance_scoreterms.sap.SapScoreMetric()
    SAPscm.set_sap_calculate_selector(binder_chain_selector)
    SAP_score=SAPscm.calculate(pose)


    return interface_sc, interface_dsasa, interface_packstat, interface_hbonds, interface_delta_unsat_hbonds, SAP_score 

def compute_interface_hydrophobicity(pose,close_residues_binder):
    '''Function to compute the propotion of hydrophobic residues the binder has in the interface'''

    total_number_of_residues=0
    hydrophobic_residues=0
    for residue_index in close_residues_binder:
        res=pose.residue(residue_index)
        if res.is_apolar() == True or res.name() == 'PHE' or res.name() == 'TRP' or res.name() == 'TYR':
            hydrophobic_residues+=1
        total_number_of_residues+=1
    try:
        binder_interface_hydrophobicity= hydrophobic_residues/total_number_of_residues 
    except:
        binder_interface_hydrophobicity=0
    return binder_interface_hydrophobicity

def compute_surface_hydrophobicity(pose):
    '''Function to compute binder surface percentage of hydrophoibic residues'''

    binder_pose = {pose.pdb_info().chain(pose.conformation().chain_begin(i)): p for i, p in zip(range(1, pose.num_chains()+1), pose.split_by_chain())}['A']
    layer_sel=pyrosetta.rosetta.core.select.residue_selector.LayerSelector()
    layer_sel.set_layers(pick_core=False, pick_boundary=False, pick_surface=True)
    surface_res_binder=layer_sel.apply(binder_pose)

    hydrophobic_count = 0
    total_count = 0 
    
    # count apolar and aromatic residues at the surface
    for i in range(1, len(surface_res_binder) + 1):
        if surface_res_binder[i] == True:
            res = binder_pose.residue(i)

            # count apolar and aromatic residues as hydrophobic
            if res.is_apolar() == True or res.name() == 'PHE' or res.name() == 'TRP' or res.name() == 'TYR':
                hydrophobic_count += 1
            total_count += 1

    surface_hydrophobicity = hydrophobic_count/total_count
    return surface_hydrophobicity

def compute_rmsd_agreement(run_number,pose):
    '''Function to compute rmsd agreement between the original design and the AF prediction'''

    #Get the design name 
    protein_name=pose.pdb_info().name()
    pattern=r'.*(run_[0-9]+_gpu_[0-9]+_design_[0-9]+).*'
    cropped_name=re.search(pattern, protein_name).group(1)
    
    #Get design structure
    #This is going to be problematic with seq_diversity since there is no PDB to compare with. Comparing with input in those situations
    original_design=f'output/run_{run_number}/{cropped_name}.pdb'
    
    if os.path.isfile(original_design):
        pose_original=pose_from_pdb(original_design)
    
        #Align the original design and the target structure using EH Baugh code
        superpose_pose_by_chain(pose, pose_original, chain='B') 

        #Compute binder RMSD 
        pose_binder=pose.split_by_chain(1)
        pose_original_binder=pose_original.split_by_chain(1)
        rmsd_binder=pyrosetta.rosetta.core.scoring.CA_rmsd(pose_binder,pose_original_binder) #Compute
        
        return rmsd_binder
    else:
        return 0

def ipsae_function(pae,L):
    if L < 27:
        d0=1
    else:
        d0=1.24*(L-15)**(1/3)-1.8
    return (1/(1+(pae/d0)**2))

def compute_ipsae(protein_name, length):
    '''Function to compute the ipSAE following Dunbrack paper (great paper!)'''
    
    #Get the run_number using re
    pattern=r'.*run_([0-9]+).*'
    run_number=re.search(pattern,protein_name).group(1)
    #Path to the json file where all the pae info is stored 
    json_file_path=glob.glob(fr'output/run_{run_number}/*{protein_name}.json')[0]

    #Load data from JSON file
    with open(json_file_path, 'r') as json_file:
        data=json.load(json_file)

    #Get the pae_binder_list   
    pae_binder=data['predicted_aligned_error'][:length]
    pae_target=data['predicted_aligned_error'][length:]

    pae_cutoff=15 # Cutoff at 15 seems a good tradeoff, can be more strict

    #Compute ipSAE aligning in the binder
    ipsae_ab_list=[]
    ipsae_ba_list=[]
    for alignment in pae_binder:
        masked_pae=[pae for pae in alignment[length:] if pae < pae_cutoff]  #residues that pass the cutoff
        L_pae=len(masked_pae)                                               #Number of residues that pass the cutoff for that alignment
        #To avoid empty lists (The later max avoids this have any effect more than better variable handling)
        if L_pae !=0:   
            ipsae_ij_list=[ipsae_function(pae,L_pae)for pae in masked_pae ]
        else:
            ipsae_ij_list=[0]
        ipsae_ab_list.append(np.mean(ipsae_ij_list))
        ipsae_ab_list.append(0)
    ipsae_ab=max(ipsae_ab_list)

    #Compute ipSAE aligning in the target
    for alignment in pae_target:
        masked_pae=[pae for pae in alignment[:length] if pae < pae_cutoff]
        L_pae=len(masked_pae)
        if L_pae != 0:
            ipsae_ji_list=[ipsae_function(pae,L_pae)for pae in masked_pae]
        else:
            ipsae_ji_list=[0]
        ipsae_ba_list.append(np.mean(ipsae_ji_list))
        ipsae_ba_list.append(0)
    ipsae_ba=max(ipsae_ba_list)
    print(f'The ipsae AB is {ipsae_ab}')
    print(f'The ipsae BA is {ipsae_ba}')
    ipsae=max([ipsae_ab, ipsae_ba,0])

    return ipsae
        

if __name__== '__main__':
    init("-mute all")

    # Protein name
    parser = argparse.ArgumentParser()
    parser.add_argument('--silent', type=str, help='silent file name')
    parser.add_argument('--run_number', type=int, help='run_number, so it can be performed with sequence diversity generation')
    args = parser.parse_args()
    
    #Load files
    file=args.silent
    run_number=args.run_number
    poses=poses_from_silent(args.silent)
    interface='A_B'

    #Load dictionary
    binding_analysis_dict={
        'description':[],
        'close_residues_target':[],
        'close_residues_binder':[],
        'CUTRE':[],
        'dSASA':[],
        'Shape_complementarity':[],
        'Packstat':[],
        'length':[],
        'SAP':[],
        'binder_int_hyd':[],
        'binder_surf_hyd':[],
        'interface_hbonds':[],
        'interface_unsat_hbonds':[],
        'RMSD':[],
        'ipSAE':[]
        }

    for pose in poses:
        t0=time.perf_counter()

        # get protein name
        protein_name=pose.pdb_info().name()

        # Get close residues
        close_residues_binder,close_residues_target,length=get_close_residues(pose)
        
        # Compute CUTRE
        try: CUTRE=compute_CUTRE(run_number, protein_name, close_residues_binder, close_residues_target)
        except:
            CUTRE=0
            print('CUTRE did not work')
        
        # Compute Rosetta operations
        shape_complementarity, dSASA, PackStat,hbonds, unsat_hbonds, SAP_score=compute_rosetta_metrics(pose, interface)

        # Binder Interface Hydrophobicity
        binder_interface_hydrophobicity=compute_interface_hydrophobicity(pose,close_residues_binder)

        # Binder Surface hydrophobicity
        surface_hydrophobicity=compute_surface_hydrophobicity(pose)
        
        #Compute RMSD between binder prediction and design
        RMSD=compute_rmsd_agreement(run_number, pose)
        
        #Compute ipSAE
        try: ipsae=compute_ipsae(protein_name, length)
        except:
            ipsae=0        
            print('ipsae did not work')

        # Info Appending
        close_residues_target_str=' '.join(map(str, list(close_residues_target)))
        close_residues_binder_str=' '.join(map(str, list(close_residues_binder)))
        binding_analysis_dict['close_residues_binder'].append(close_residues_binder_str)
        binding_analysis_dict['close_residues_target'].append(close_residues_target_str)
        binding_analysis_dict['length'].append(length)
        binding_analysis_dict['CUTRE'].append(CUTRE)
        binding_analysis_dict['description'].append(protein_name)
        binding_analysis_dict['dSASA'].append(dSASA)
        binding_analysis_dict['Shape_complementarity'].append(shape_complementarity)
        binding_analysis_dict['Packstat'].append(PackStat)
        binding_analysis_dict['SAP'].append(SAP_score)
        binding_analysis_dict['binder_int_hyd'].append(binder_interface_hydrophobicity)
        binding_analysis_dict['binder_surf_hyd'].append(surface_hydrophobicity)
        binding_analysis_dict['interface_hbonds'].append(hbonds)
        binding_analysis_dict['interface_unsat_hbonds'].append(unsat_hbonds)
        binding_analysis_dict['RMSD'].append(RMSD)
        binding_analysis_dict['ipSAE'].append(ipsae)

    # Save Info
    df=pd.DataFrame(binding_analysis_dict)

    ## Get sc_scorings

    sc_df=get_sc_scorings(run_number).drop('SCORE:', axis=1)
    merged_df=pd.merge(df, sc_df, on='description')
    print(merged_df)
    try:
        existing_data = pd.read_csv('Scoring_Stats.csv')
        updated_data = pd.concat([existing_data,merged_df], ignore_index=True)
        updated_data.to_csv('Scoring_Stats.csv', index=False)
    except FileNotFoundError:
        merged_df.to_csv('Scoring_Stats.csv', index=False)

    t_final=time.perf_counter()
    print(f'{pose} PyRosetta analysis took {t_final-t0} seconds to run')