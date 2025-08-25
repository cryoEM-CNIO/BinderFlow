#!/usr/bin/env python3

'''
Script to align the template structure to the one used for the binder generation and substitute it in order to perform the next steps using the 
template structure (pMPNN, AF2-IG)

Input:
    --input_dir --> directory where the designs are stored
    --chain --> Chain to modify and substitute (Always B for RFD)
    --template --> Path to the template structure
    --reference -> Path to a reference structure

Output:
    It generates one run_X_design_N_substituted.pdb per input pdb

'''
import numpy
from Bio import Align
from Bio import PDB, SeqIO
from Bio.SeqIO import PdbIO
from Bio.PDB import Superimposer
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.StructureAlignment import StructureAlignment
from Bio.Align import Alignment
from Bio.Seq import Seq
from Bio.PDB.DSSP import DSSP
import argparse
import glob
import os
import warnings
from Bio import BiopythonWarning
# from pyrosetta import * 
import re
warnings.simplefilter('ignore', BiopythonWarning)


###############################################
#-------------- PDB ORDER FUNCTION -----------#
###############################################
'''RFD needs the pdb to be ordered to use it. This functions make sure the PDB is correctly ordered'''
def order_pdb(pdb_file):
    res_chains=[]
    chain_dictionary={}
    with open(pdb_file, 'r') as file:
        for line in file.readlines():
            if line.startswith('ATOM'):
                res_chain=line[22]
                if res_chain not in res_chains:
                    res_chains.append(res_chain)
                if res_chain not in chain_dictionary.keys():
                        chain_dictionary[res_chain]=[]
                chain_dictionary[res_chain].append(line)
    res_chains_sorted=res_chains.copy()
    res_chains_sorted.sort()
    if res_chains_sorted != res_chains:
        with open(pdb_file, 'w') as file:
            for chain in res_chains_sorted:
                for line in chain_dictionary[chain]:
                    file.write(line)

###############################################
#-------------- ALIGNMENT FUNCTIONS ----------#
###############################################
'''This set of functions do all the alignments needed to replace your target with a larger target to perform the sequence assignment and scoring'''
#Function to extract the fasta seqs of the hits 
def extract_seq(input_file,chain_selected,structure):
    #Extract sequence 
    with open(input_file, 'r') as pdb_file:
        for record in PdbIO.AtomIterator(pdb_file,structure):
            chain = record.annotations['chain']
            if chain == chain_selected:
                sequence = record.seq
                modified_seq = Seq(str(sequence).replace("X", "-"))
                return modified_seq

#Function to aligne the sequence of the template and target
def seq_alignment(template, moving): 
    aligner=Align.PairwiseAligner()
    alignments =aligner.align(template, moving)
    return alignments

#Function to perform the structural alignment
def structure_alignment(moving, template, chain):

    pdb_parser= PDBParser(PERMISSIVE=1)

    moving_id, template_id = (file.split('/')[-1].split('.')[0] for file in [moving, template])
    moving_model=pdb_parser.get_structure(moving_id, moving)
    moving_model_target=moving_model[0]['B']
    template_model=pdb_parser.get_structure(template_id, template)
    template_model_target=template_model[0]['B']
    
    #Get sequences
    moving_sequence=extract_seq(moving,chain,moving_model)
    template_sequence=extract_seq(template, chain, template_model)

    
    #Align and get the coordinates
    alignments=seq_alignment(template_sequence, moving_sequence)
    alignment_sub=alignments[0]
    coordinates=alignment_sub.coordinates
    
    #Align for the structural alignment
    alignment = Alignment([template_sequence, moving_sequence], coordinates)

    structure_alignment_maps=StructureAlignment(alignment, template_model_target, moving_model_target).get_maps()[1]

    ref_ids=[]
    moving_ids=[]
    ref_atoms=[]
    moving_atoms=[]

    for key,value in structure_alignment_maps.items():
        if value:
            ref_ids.append(value.get_id()[1])
            moving_ids.append(key.get_id()[1])

    for ref_chain in template_model[0]:
        for ref_res in ref_chain:
            if ref_res.get_id()[1] in ref_ids:
                ref_atoms.append(ref_res['CA'])

    for moving_chain in moving_model[0]:
        if moving_chain.get_id() == 'B':
            for moving_res in moving_chain:
                if moving_res.get_id()[1] in moving_ids:
                    moving_atoms.append(moving_res['CA'])

    superimposer=Superimposer()
    superimposer.set_atoms(ref_atoms, moving_atoms)
    superimposer.apply(moving_model.get_atoms())
    
    return moving_model


#################################################
#--------------SUBSTITUTION FUNCTIONS-----------#
#################################################
'''Functions to substitute the target use in RFD with the template for the rest of the process'''

# Function to remove the chain B from the pdb
def remove_chain(structure, chain_id):
    for model in structure:
        for chain in list(model):
            if chain.id == chain_id:
                model.detach_child(chain.id)

#Function to add the template as new chain in the pdb file 
def substitute_chain(pdb_path, chain, template):
    pdb_id = os.path.splitext(os.path.basename(pdb_path))[0]
    
    pdb_parser = PDB.PDBParser(QUIET=True)
    
    # Load the first structure and remove the specified chain
    structure1 = pdb_parser.get_structure(pdb_id, pdb_path)
    remove_chain(structure1, chain)

    # Load the second structure and get the desired chain
    for model in template:
        for template_chain in model:
            if template_chain.id == chain:
                new_chain = template_chain

    # Add the new chain to the first structure
    for model in structure1:
        model.add(new_chain)

    return structure1

#################################################
#--------------- FILTER FUNCTIONS -------------#
#################################################
'''These are a series of functions we use to filter the results from the initial backbone generation and use more efficiently computational resources'''

#This function discards the binder if it makes steric clashes with the template once the chains are substituted
def check_clashes(structure):
    structure_binder=structure[0]['A']
    structure_target=structure[0]['B']
    residues_binder = [r for r in structure_binder.get_residues()] 
    residues_target = [r for r in structure_target.get_residues()] 

    # Calculate the distance between the alpha carbons for a pair of residues
    for residue_i in residues_binder:
        for residue_j in residues_target: 
            residue_i_coord = residue_i["CA"].get_coord()
            residue_j_coord = residue_j["CA"].get_coord()
            distance=numpy.linalg.norm(residue_i_coord-residue_j_coord)
            if distance < 0.5:
                print(f'CLASHES DETECTED, STRUCTURE {structure.get_id()} WONT BE USED')
                return 0
            else:
                continue
    print('NO CLASHES DETECTED')
    return 1

#This function filters the designs based on dssp. If the designs are predicted to have a single helix or a hairpin structure they are discarded
def filter_by_dssp(pdb_file):
    
    #define the pattern
    hairpin=r'^[-ST]*[HIG]+[-ST]+[HIG]+[-ST]*$'
    full_helix=r'^[-ST]*[HIG]+[-ST]*$'

    p=PDBParser()
    structure=p.get_structure('binder', pdb_file)
    model=structure[0]
    dssp=DSSP(model, pdb_file)

    #Get the ss
    ss_string=''
    for key in dssp.keys():
        if key[0] == 'A':
            ss_string+=dssp[key][2]
        else:
            break

    if re.match(hairpin, ss_string):
        print(f'The structure {pdb_file} is a hairpin, not processing it ')
        return 0
    elif re.match(full_helix, ss_string):
        print(f'The structure {pdb_file} is a long helix, not processing it ')
        return 0
    else:
        print(f'The structure {pdb_file} is not a hairpin, processing ... bip,bup,bop (robot sounds) ')
        return 1



#################################################
#----------------- SAVE FUNCTIONS --------------#
#################################################

'''Function to save the new pdb file with the substituted chain. Only those that pass the filtering are saved'''

def save_protein_substituted(structure,output_dir):
    pdb_id=structure.get_id()
    output_path = os.path.join(output_dir, f"{pdb_id}_substituted.pdb")
    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(output_path)
    print("Saved "+output_path)
    return output_path
    
#################################################
#----------------- FIXED FUNCTIONS -------------#
#################################################
def residues_length_added(input_path, template_path, partial_diff='True'):
    '''
    When scaffolding is performed, some residues can be added to the N-ter, changing the numbering of the residues,
    This function should recognize the original sequence (if there is chain A), and count the number of residues added to the N-ter
    If the residues are added to the C-ter, the function should not be used
    '''
    if partial_diff == 'True':
        number_of_residues_added=0
    else:
        #extract structure
        template_structure = PDB.PDBParser(QUIET=True).get_structure('template', template_path)
        input_structure = PDB.PDBParser(QUIET=True).get_structure('input', input_path)
        #Extract sequence 
        with open(input_path, 'r') as pdb_file:
            for record in PdbIO.AtomIterator(pdb_file,input_structure):
                chain = record.annotations['chain']
                if chain == 'A':
                    binder_sequence = record.seq
        with open(template_path, 'r') as template_file:
            for record in PdbIO.AtomIterator(template_file,template_structure):
                chain = record.annotations['chain']
                if chain == 'A':
                    template_sequence = record.seq
        number_of_residues_added=len(str(binder_sequence).split(str(template_sequence))[0])
    return int(number_of_residues_added)


'''Function to add the fixed remark to the pdbs in order to fix some residues during pMPNN'''
def add_fixed_residues(input_path,template_path,output_path, residues,partial_diff):
    residues_to_fix = []
    N_ter_residues=residues_length_added(input_path, template_path,partial_diff)
    if residues is not None and residues != 'None' and residues != '[]':
        residues_list = residues.strip('[').strip(']').strip().split(',')  # No problem if there are no commas in the input
        for resi in residues_list:
            try:
                residues_to_fix.append(str(int(resi)+N_ter_residues)+1)
            except ValueError:
                for resi_range_id in range(int(resi.split('-')[0]), int(resi.split('-')[1]) + 1):
                    residues_to_fix.append(int(resi_range_id)+int(N_ter_residues)+1)
        
        # Read the file contents
        with open(output_path, 'r') as read_file:
            contents = read_file.readlines()
        
        # Find the last occurrence of "END"
        last_end_index = None
        for i in range(len(contents) - 1, -1, -1):
            if contents[i].strip() == "END":
                last_end_index = i+1
                break
        # If "END" is not found, set last_end_index to the last TER
        if last_end_index is None:
            for i in range(len(contents) - 1, -1, -1):
                if contents[i].strip() == "TER":
                    last_end_index = i+1
                    break
        # Insert the fixed remarks before the last "TER"
        if last_end_index is not None:
            for residue_id in residues_to_fix:
                contents.insert(last_end_index, f"REMARK PDBinfo-LABEL:{residue_id: >5} FIXED\n")
        
        # Write the updated contents back to the file
        with open(output_path, 'w') as write_file:
            write_file.writelines(contents)
        
        print(f'Residues {[i for i in residues_to_fix]} fixed in the pdb')

def main():
    parser = argparse.ArgumentParser(description="Substitute a chain in PDB files.")
    parser.add_argument("--chain", required=True, help="Chain ID to substitute from the second structure")
    parser.add_argument("--template", required=True, help="Path to the second PDB structure, this is the template from the microrun")
    parser.add_argument("--run", required=True, help='Run that is running')
    parser.add_argument("--t", required=False, default='' , help='Indicator of the GPU in which it is running')
    parser.add_argument("--residues", required=False, default='None', help='List of residues of the binder to fix (Useful for scaffolding)')
    parser.add_argument("--partial_diff", required=False, default='True', help='It is a partial diffusion run ?')
    args = parser.parse_args()

    '''
    This can be a little bit confusing, but what this script does is: 
        -aligning the template structure to the desings 
        -substitute the target with the aligned template

    so in the template alignign, the template works as moving PDB and the original design as template.
    Probably this nomenclature should be modified to facilitate debugging
    '''
    # Variable definition
    io_path=f"output/run_{args.run}/"
    reference=f"{io_path}/run_{args.run}_gpu_{args.t}_design_0.pdb"
    # This first function checks that the pdb has the order it should have for Biopython to understand it 
    order_pdb(args.template)
    template_aligned=structure_alignment(moving=args.template, template=reference, chain=args.chain)
    for pdb_path in glob.glob(f"{io_path}/run_{args.run}_gpu_{args.t}_design_*.pdb"):
        sub_structure=substitute_chain(pdb_path, args.chain, template_aligned)
        clashes,hairpin=(int(check_clashes(sub_structure)), int(filter_by_dssp(pdb_path)))
        if clashes+hairpin == 2:
            output_path=save_protein_substituted(sub_structure,io_path)
            if args.residues != 'None' and args.residues != '':
                add_fixed_residues(pdb_path,args.template, output_path, args.residues, args.partial_diff)
        else:
            continue




if __name__ == '__main__':
    # init("-mute all")
    main()
