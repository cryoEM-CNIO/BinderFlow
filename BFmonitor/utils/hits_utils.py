import pandas as pd
import os
import re
import warnings
import time 
import numpy as np
import glob 
import argparse
import json
import subprocess
import random
import math
import plotly.express as px
import plotly.graph_objects as go
import warnings
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# SeqIO raises some warnings that not affect the outcome, this silence them
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

# Function to get best hit name
def get_hit_names(filtered_df, xaxis_value):
    '''
    Function to get all the hit names. This is employ to get the hit names for the ngl representation dropdown
    Probably could be fused with the enxt function

    Input: 

    filtered_df --> DF filtered with only those designs that fulfill all the filters

    xaxis_value --> Variable which is employ to order the designs in the dropdown
    
    Output:

    hit_names --> A list of the hit names, ordered in ascending order using the x variable
    '''

    if not filtered_df.empty:
        sorted_df = filtered_df.sort_values(by=xaxis_value, ascending=True)
        hit_names = sorted_df['description'].tolist()
        return hit_names
    else:
        return ["No hits found using current filters"]
    
#Function to get hits file path and name for representation 
def get_design_file_path_and_name(hits_names, directory, input_pdb_path):
    '''
    Gets a design file path and name for its representation in NGL.

    Inputs:
    - hits_names: identifier of the hit name (e.g., run_1_design_2_substituted or run_1_gpu_0_design_5)
    - directory: base directory, e.g., /path/to/output
    - input_pdb_path: path to the input PDB file for viewing the original structure in NGL

    Outputs:
    - data_path: full path to the directory containing the file
    - filename: filename without extension for NGL use
    '''

    description = hits_names

    # Pattern for old format: run_1_design_2_substituted
    pattern_old = r"(run_\d+)_design_(\d+_substituted)"
    # Pattern for new format: run_1_gpu_0_design_5
    pattern_new = r"(run_\d+)_gpu_(\d+)_design_(\d+_substituted)"
    # Input case
    pattern_input = r"Input"

    match_old = re.match(pattern_old, description)
    match_new = re.match(pattern_new, description)
    match_input = re.match(pattern_input, description)

    if match_input:
        input_path = os.path.join(directory, input_pdb_path)
        print(input_path)
        data_path = '/'.join(glob.glob(input_path)[0].split('/')[:-1])
        filename = glob.glob(input_path)[0].split('/')[-1].split('.')[0]
        return data_path, filename

    elif match_new:
        directory = os.path.join(directory, 'output')
        run_part,gpu_part,design_part = match_new.groups()
        data_path = os.path.join(directory, run_part)
        filename = f"{run_part}_gpu_{gpu_part}_design_{design_part}"
        return data_path, filename

    elif match_old:
        directory = os.path.join(directory, 'output')
        run_part, design_part = match_old.groups()
        data_path = os.path.join(directory, run_part)
        filename = f"{run_part}_design_{design_part}"
        return data_path, filename

    else:
        print("Invalid description format")
        return None, None

#Function to add the stats of the hits to the PDB, as they are added in after AF2IG
#Makes PD comparison more comfortable
def add_stats_to_pdb(description,directory):
    '''
    add the stats to the pdb, in the same fashion af2IG does (adding the name and the metric separated by a whitespace)

    Input:

    description --> identifier of the hit to record the metrics inside them

    directory --> directory where we are working and were the Scoring_Stats.csv is

    Output:

    Writes the data at the end of the pdb file    
    '''

    df_rosetta = pd.read_csv(f'{directory}/Scoring_Stats.csv')
    design_metrics = df_rosetta[df_rosetta['description'] == description]
    design_metrics = design_metrics.drop(['close_residues_target', 'close_residues_binder'], axis=1)
    #The fastas has an abbreviated form of the name, so you have to crop it
    description_short=description.split('_af2pred')[0]
    pdb_path = f'{directory}/hits/{description}.pdb'
    fasta_path=f'{directory}/hits/fastas/{description_short}.fasta'
    mw, ip, extinction_coefficient=param_stats(fasta_path)
    with open(pdb_path, 'a') as file:
        for column in design_metrics.columns:
            if column != 'description':  # Skip writing the 'description' column itself
                value = design_metrics[column].iloc[0]  # Extract the first value from the series
                file.write(f'{column} {round(value,4)}\n')
        file.write(f'molecular_weight {round(mw,4)}\n')
        file.write(f'isoelectric_point {round(ip,4)}\n')
        file.write(f'extinction_coefficient {round(extinction_coefficient,4)}\n')

#Function to extract the fasta seqs of the hits 
def extract_fasta_seq(description, fasta_dir):
    '''
    Extract fasta sequences for order

    Input:

    description --> Hit name

    Output:

    {description}.fasta --> Fasta sequence of the hit
    '''

    # Define variables
    input_file=f'{description}.pdb'
    pattern = r'.*(run_\d+.*_dldesign_\d+)'
    directory = os.path.dirname(input_file)

    #Get fasta names 
    try:
        fasta_name = re.search(pattern, input_file).group(1)
    except AttributeError:
        fasta_name = input_file[:-4]
    
    
    #Extract sequence 
    with open(input_file, 'r') as pdb_file:
        for record in SeqIO.parse(pdb_file, 'pdb-atom'):
            chain = record.annotations['chain']
            if chain == 'A':
                sequence = record.seq

    #write things down

    fasta_file = f'{fasta_dir}/{fasta_name}.fasta'
    with open(fasta_file, 'w') as fasta:
        fasta.write(f'>{fasta_name}\n{sequence}\n')

    return fasta_file
# Function to extract the DNA seq with the desired characteristics
def extract_dna_seq_CT(multifasta_file, output_dir, params_dict):
    '''
    Extract optimized dna seqs for order

    Input:

    input --> Fasta file with the aminoacidic sequence

    output --> Folder in which the dna seqs are gonna be stored

    params --> Dict with all the params which are:

        organism --> Organism desired for codon optimization

        met --> Whether to add a methionine at the start of the sequence or not (default: True)

        overhang_5 --> Overhang sequence to add at the 5\' (default=Nothing) 

        overhang_3 --> Overhang sequence to add at the 3\' (default=Nothing)

        length --> Number of random bases that must be added to reach a certain size. Half of the bases are added at each extreme, between the overhangs and the sequence (default=0)        
    
        GC --> GC content of the random sequence, express in % (default=50)

        enzyme --> Restriction enzymes to check in the sequence
    
    Output:

    Nothing, it writes the DNA seq at the folder specified in the output key
    '''

    stdout = ''

    # Setting the input from the command line arguments
    input = multifasta_file
    output = output_dir
    organism = params_dict['organism']
    met = params_dict['add-met']
    overhang_5 = params_dict['five_prime_overhang']
    overhang_3 = params_dict['three_prime_overhang']
    length = params_dict['random_sequence']
    GC = params_dict['GC_content']
    enzyme = params_dict['enzyme']

    stdout += 'Selected parameters:\n'
    stdout += f'Organism: {organism}\n'
    stdout += f'Add Methionine: {met}\n'
    stdout += f'5\' Overhang: {overhang_5}\n'
    stdout += f'3\' Overhang: {overhang_3}\n'
    stdout += f'Length to reach: {length}\n'
    stdout += f'GC content of random sequence: {GC}\n'
    stdout += f'Enzymes to check: {enzyme}\n'


    #Define a dictionary which is later gonna be used to store all the order information
    order_dictionary={
        'design_name':[],
        'organism':[],
        'met_added':[],
        '5\' overhang':[],
        '3\' overhang':[],
        'fasta':[],
        'dna_seq':[],
        'binder_order':[]
    }

    # Collect all the ids and dna sequences
    id_list=[]
    dna_seq_list=[]

    # Extract the sequences from the multifasta
    with open (input, 'r') as fasta_file:
        for line in fasta_file.readlines():
            if line.startswith('>'):
                id = line[1:].strip()
            else:
                sequence=line
        stdout += f'Aminoacid Sequence of {id}:\n{sequence}\n'
        #Add initial met 
        if met == True:
            sequence='M'+sequence
            stdout += 'Methionine added at the start of the sequence\n'
            stdout += f'New sequence: {sequence}\n'

        # Variable preparation
        codon_transformer_output=os.path.join(output, id+'_RevTrans.fasta')

        all_ok=False
        counter=1
        while all_ok == False:
            
            stdout += f'Running CodonTransformer for {id} for {counter} time\n'

            # run CodonTransformer
            command = (
                f"python3 /emdata_fast/cchacon/protein_design/monitoring/utils/dna_extraction/CT/CodonTransformer_seq.py --protein '{sequence}' --organism '{organism}'"
            )

            # Execute the command
            process = subprocess.Popen(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                shell=True,
                text=True
            )
            seq_generated, stderr = process.communicate()

            dna_seq=seq_generated.strip().lower()

            ## Design random sequences for the 3 and 5 ends to reach the desired size
            length_to_add=int(length)-(len(dna_seq)+len(overhang_3)+len(overhang_5))
            if length_to_add > 0:
                stdout += f'Length of the sequence {len(dna_seq)} is smaller than the desired length {length}, ' + \
                          f'A random sequence {length_to_add} residues long will be added'
            random_sequence_5,random_sequence_3, out_info= RandomSequenceGenerator(length_to_add, GC)

            stdout += f'\n{out_info}\n'

            full_seq=f'{overhang_5}' + f'{random_sequence_5}' + f'{dna_seq}'+f'{random_sequence_3}' + f'{overhang_3}'
            full_seq,all_ok=check_enzyme_cut(enzyme,dna_seq)

            if all_ok == False:
                stdout += f'Enzyme restriction site found in {id}, regenerating DNA sequence...\n'
            counter +=1
        stdout += f'CodonTransformer output correctly generated for {id}:\n{full_seq}\n'

    ## Open the RevTrans fasta_file

        with open (codon_transformer_output, 'w') as dna_seq_fasta:
            dna_seq_fasta.write(f'{id}\n{full_seq}\n')


    return stdout


def extract_dna_seq_JB(dnas_dir, multifasta_file, params_dict):
    '''
    Function to optimize functions using JohnBercow method

    Inputs:

    dnas_dir --> Directory where the output files are going to be saved,

    multifasta_file --> Fasta file with the aminoacidic sequences to be optimized,
    
    params_dict --> Dictionary with the parameters for the JohnBercow command, which are:
        
        Mandatory:

        order_name --> Name of the order, used to name the output files

        gg_vector -->  Golden Gate vector employed for cloning from the  

        species --> Organism for which the sequences are optimized,

        design_prefix --> Prefix for the design, used to name the output files

        design_id --> Identifier for the design, used to name the output files

        Optional:

        idt_score --> IDT score for the sequences, default is 7

        starting_kmers_weight --> Weight for the starting kmers, default is 10

        n_domesticator_steps --> Number of domesticator steps, default is 10

        max_attempts --> Maximum number of attempts to find a valid sequence, default is 20

        max_length --> Maximum length of the sequence, default is 1500

        skip_idt_radio --> Whether to skip the IDT score calculation, default is False

        print_heterooligomers --> Whether to print heterooligomers, default is False

        no_layout --> Whether to skip the layout generation, default is False

        no_plasmids --> Whether to skip the plasmid generation, default is False

        verbose --> Whether to print verbose output, default is True

        echo --> To be completed
    
    Output:

    Executes the JohnBercow command, given all the needed files for IDTs DNA sequence ordering
    '''


    # Get the parameters from the params_dict
    order_name = params_dict['order_name']
    design_prefix = params_dict['design_prefix']
    design_id = params_dict['design_id']
    gg_vector = params_dict['golden_gate_vector_JB']
    species = params_dict['species_JB']
    idt_score = params_dict['idt_score']
    starting_kmers_weight = params_dict['starting_kmers_weight']
    n_domesticator_steps = params_dict['n_domesticator_steps']
    max_attemps = params_dict['max_attempts']
    max_length = params_dict['sequence_max_length']
    skip_idt = params_dict['skip_idt_radio']
    print_hto = params_dict['print_heterooligomers']
    no_layout = params_dict['no_layout']
    no_plasmids = params_dict['no_plasmids']
    verbose = params_dict['verbose']
    echo = params_dict['echo']

    #Check the inputs have been provided correctly
    stdout = check_all_conditions(
        order_name, design_prefix, design_id, gg_vector, species
    )
    if stdout:
        return stdout
    
    # Get into a list the flags to expand
    optional_flags = ['--skip_idt', '--print_hto', '--no_layout', '--no_plasmids', '--verbose', '--echo']
    selected_options = [skip_idt, print_hto, no_layout, no_plasmids, verbose, echo]
    expanded_flags = [flag for flag, selected in zip(optional_flags, selected_options) if selected]

    # Define the command to run JohnBercow

    command = (
        f"python3 -u /emdata_fast/cchacon/protein_design/monitoring/utils/dna_extraction/JB/JohnBercow.py --order_fasta '{multifasta_file}' --order_name '{order_name}' "+
        f"--gg_vector '{gg_vector}' --output_folder '{dnas_dir}' --species '{species}' --design_prefix '{design_prefix}' --design_id '{design_id}' --idt_score '{idt_score}' "+
        f"--starting_kmers_weight '{starting_kmers_weight}' --n_domesticator_steps '{n_domesticator_steps}' --max_attempts '{max_attemps}' --max_length '{max_length}' " +
        f"{' '.join(expanded_flags)}" 
    )

    # Execute the command
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        shell=True,
        text=True
    )
    stdout, stderr = process.communicate()

    # Check for errors
    if process.returncode != 0:
        print("Error:", process.stderr)
        raise RuntimeError("JohnBercow process failed")
    return stdout

def get_gg_vectors(path_to_lab_vectors):
    '''
    Returns a list of Golden Gate vectors used for cloning.
    Needed for the dropdown list

    Input:
    path_to_lab_vectors (str): Path to the folder containing the lab vectors.

    Output:
    List of strings representing the Golden Gate vectors.
    '''

    lab_vectors = glob.glob(os.path.join(path_to_lab_vectors, '*.fa'))

    return [os.path.basename(v).split('_')[0] for v in sorted(lab_vectors)]


def RandomSequenceGenerator(length, GC=50):
    '''
    Generate random DNA sequences with the desired GC proportion. The sequence is added symmetrically 
    at the 3' and 5' ends, with half the required length at each terminus.

    Input:
        length (int): Length of the random sequence to add. It is split equally between 5' and 3' ends.
        GC (float): Desired GC content in percentage (default=50).

    Output:
        tuple: (sequence_5, sequence_3) Random sequences for the 5' and 3' ends.
    '''
    if length <= 0:
        return '', '', ''

    # Convert GC content percentage to probabilities
    GC_prob = GC / 100 / 2
    AT_prob = (1 - GC_prob * 2) / 2
    nucleotide_weights = [AT_prob, GC_prob, AT_prob, GC_prob]  # For 'a', 'g', 't', 'c'

    # Calculate the lengths for 5' and 3' ends
    length_5 = (length + 1) // 2  # Round up for the 5' end
    length_3 = length // 2       # Remaining for the 3' end
    
    
    # Initialize sequences
    sequence_5 = []
    sequence_3 = []

    out_info = f'Creating a random sequence {length} residues long'
    
    # Generate sequence for the 3' end
    for _ in range(length_3):
        sequence_3.append(random.choices(['a', 'g', 't', 'c'], weights=nucleotide_weights)[0])

    # Generate sequence for the 5' end
    for _ in range(length_5):
        if len(sequence_5) >= 2 and sequence_5[-2:] == ['a', 't']:  # Avoid ATG motif
            sequence_5.append(random.choices(['a', 'g', 't', 'c'], weights=[AT_prob, 0, AT_prob, GC_prob * 2])[0])
        else:
            sequence_5.append(random.choices(['a', 'g', 't', 'c'], weights=nucleotide_weights)[0])
    
    # Convert list to string and return
    return ''.join(sequence_5), ''.join(sequence_3), out_info




def generate_order_csv(extraction_list, hits_folder):
    '''
    Generates a csv which all the information required for ordering (in our case)

    Input:

    extraction_list --> List of designs whose structure and sequence is going to be extracted

    hits_folder --> Folder where all the hits information is being saved
    Output:

    order.csv -->  A csv file which details the path to the original structure, as well as to the aminoacidic and nucleotidic sequences, 
                   and the design identifier, aminoacidic sequence and nucleotidic sequence
    
    '''
    
        # Load dna_seqs
    data_dictionary={
        'design_path':[],
        'description':[],
        'fasta_path':[],
        'protein_sequence':[],
        'dna_path':[],
        'dna_sequence':[],
        'molecular_weight':[],
        'isoelectric_point':[],
        'molar_extinction_coefficient':[]
    }

    pattern = r'.*(run_\d+_design_\d+.*_dldesign_\d+).*'
    pattern_run=r'.*(run_\d+)_.*'


    #Load protein_seq
    for description in extraction_list:
        protein_path=f'{hits_folder}/fastas/{re.search(pattern,description).group(1)}.fasta' 
        with open(protein_path, 'r') as protein_file:
            for line in protein_file.readlines():
                if not line.startswith('>'):
                    protein_seq=line
        mw,ip,extinction_coefficient=param_stats(protein_path)
        #Load DNA seq
        dna_path=f'{hits_folder}/dna_seqs/{re.search(pattern,description).group(1)}_RevTrans.fasta'
        with open(dna_path, 'r') as dna_file:
            for line in dna_file.readlines():
                if not line.startswith('>'):
                    dna_seq=line

        run_number=re.search(pattern_run, description).group(1)
        path=f'output/{run_number}/{description}'
        data_dictionary['design_path'].append(path)
        data_dictionary['description'].append(re.search(pattern, description).group(1))
        data_dictionary['fasta_path'].append(protein_path)
        data_dictionary['protein_sequence'].append(protein_seq)
        data_dictionary['dna_path'].append(dna_path)
        data_dictionary['dna_sequence'].append(dna_seq)
        data_dictionary['molecular_weight'].append(mw)
        data_dictionary['isoelectric_point'].append(ip)
        data_dictionary['molar_extinction_coefficient'].append(extinction_coefficient)
    
    data_df=pd.DataFrame(data_dictionary)
    data_df.to_csv(f'{hits_folder}/order.csv')


def param_stats(fasta):
    '''
    This function is meant to check if there are aromatic residues in the protein to facilitate protein detection and quantification
    
    Input:

    fasta --> Fasta file with the protein sequence

    Output:

    mw --> molecular weight of the protein
    ip --> isoelectric point of the protein
    extinction_coefficient --> Extinction coefficient of the protein 

    '''

    #Extract the sequence

    with open(fasta, 'r') as file:
        for line in file.readlines():
            if not line.startswith('>'):
                sequence=line

    #open it in biopython suite
    X= ProteinAnalysis(sequence)

    #param calculation

    mw=X.molecular_weight()
    ip=X.isoelectric_point()
    extinction_coefficient=X.molar_extinction_coefficient()[0]

    ## Check the number of W to determine if the extinction coefficient is faithful

    return mw, ip, extinction_coefficient 


def check_enzyme_cut(enzyme, sequence):
    '''
    This function checks possible enzyme restriction cutting sites and replaces problematic regions with equivalent sequences.

    Input:
    enzyme --> List of enzymes to check.
    sequence --> DNA sequence to check.

    Output:
    new_sequence --> Sequence without the restriction sites.
    '''
    amino_acid_codons = {
        'A': ['GCT', 'GCC', 'GCA', 'GCG'],  # Alanine
        'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],  # Arginine
        'N': ['AAT', 'AAC'],  # Asparagine
        'D': ['GAT', 'GAC'],  # Aspartic acid
        'C': ['TGT', 'TGC'],  # Cysteine
        'Q': ['CAA', 'CAG'],  # Glutamine
        'E': ['GAA', 'GAG'],  # Glutamic acid
        'G': ['GGT', 'GGC', 'GGA', 'GGG'],  # Glycine
        'H': ['CAT', 'CAC'],  # Histidine
        'I': ['ATT', 'ATC', 'ATA'],  # Isoleucine
        'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],  # Leucine
        'K': ['AAA', 'AAG'],  # Lysine
        'M': ['ATG'],  # Methionine (Start codon)
        'F': ['TTT', 'TTC'],  # Phenylalanine
        'P': ['CCT', 'CCC', 'CCA', 'CCG'],  # Proline
        'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],  # Serine
        'T': ['ACT', 'ACC', 'ACA', 'ACG'],  # Threonine
        'W': ['TGG'],  # Tryptophan
        'Y': ['TAT', 'TAC'],  # Tyrosine
        'V': ['GTT', 'GTC', 'GTA', 'GTG'],  # Valine
        '*': ['TAA', 'TAG', 'TGA']  # Stop codons
    }
    restriction_enzymes = {
        'EcoRI': 'GAATTC',
        'BamHI': 'GGATCC',
        'HindIII': 'AAGCTT',
        'NotI': 'GCGGCCGC',
        'XhoI': 'CTCGAG',
        'PstI': 'CTGCAG',
        'SacI': 'GAGCTC',
        'KpnI': 'GGTACC',
        'SmaI': 'CCCGGG',
        'XbaI': 'TCTAGA',
        'SpeI': 'ACTAGT',
        'NcoI': 'CCATGG',
        'SalI': 'GTCGAC',
        'ApaI': 'GGGCCC',
        'HaeIII': 'GGCC',
        'AluI': 'AGCT',
        'TaqI': 'TCGA',
        'BglII': 'AGATCT',
        'ClaI': 'ATCGAT',
        'MluI': 'ACGCGT',
        'BsaI': 'GGTCTC'
    }
    # Get restriction sequences for selected enzymes
    selected_re = [restriction_enzymes[y] for y in enzyme]
    # Convert sequence to a mutable list
    list_sequence = list(sequence.upper())
    # Iterate through the sequence
    for i in range(len(sequence) - 5):  # Ensuring enough length for restriction sites
        sliding_window = sequence[i:i+6]  # Most restriction sites are 6 bp long
        found_enzyme = next((enzyme_re for enzyme_re in selected_re if enzyme_re in sliding_window), None)
        if found_enzyme is not None:
            problematic_enzyme=list(restriction_enzymes.keys())[list(restriction_enzymes.values()).index(found_enzyme)]
            print(f'Sequence of cut for {problematic_enzyme} has been found')
            #Find the codon at position i+6
            codon_index=math.floor((i+3)/3) #zero indexed
            codon_identity=sequence[codon_index*3:codon_index*3+3]   
            # Find the amino acid corresponding to the codon
            residue = None
            for aa, codons in amino_acid_codons.items():
                if codon_identity in codons:
                    residue = aa
                    if residue in ['M', 'W']:
                        codon_index-=1
                        codon_identity=sequence[codon_index*3:codon_index*3+3]
                        for aa, codons in amino_acid_codons.items():
                            if codon_identity in codons:
                                residue = aa
                    print(f'Searching for a new codon for {residue}')
                    break
            if residue != '*':
                # Select an alternative codon (excluding the original one)
                alternative_codons = [codon for codon in amino_acid_codons[residue] if codon != codon_identity] 
                if alternative_codons:
                    alternative_codon = random.choice(alternative_codons)  # Pick a random alternative codon
                    list_sequence[codon_index*3:codon_index*3+3] = alternative_codon  # Replace codon
                    print(f'Replacing {"".join(codon_identity)} with {"".join(alternative_codon)} ')
                    sequence=''.join(list_sequence)
                else:
                    print('M and W next to each other, regenerating DNA seq to avoid restriction site')
                    return '',False
            else:
                print('The codon selected is a stop codon, must not be modified, redesigning DNA seq') # Probably I should select other residue
                return '',False
    # Convert list back to string
    new_sequence = ''.join(list_sequence)
    return new_sequence.lower(),True

def extract_pdbs(description, working_dir, pdbs_folder):
    '''
    Function to extract the pdbs.

    Input
    description --> Description of the pdbs to extract
    working_dir --> Directory where the pdbs are located
    pdbs_folder --> Folder containing the pdb files

    Output
    Extract the pdbs into the specified folder.
    '''
    # Check if it is new or old format

    run_number, gpu_number, design_number= None, None, None
    match = re.match(r'run_(\d+)_gpu_(\d+)_design_(\d+).*', description)

    if match:
        run_number, gpu_number, design_number = match.groups()
    
    # If gpu_number does not exist, old format

    if gpu_number is None:
        run_number = re.search(r'run_\d+',description).group()
        design_number = math.floor(int(re.search(r'run_\d+_design_(\d+).*',description).group(1))/10)

        # Check if it is from an even older run and generate the paths for the af2 files

        if design_number != 0:
            af2_file=f'{working_dir}/output/{run_number}/{run_number}_design_{design_number}_input_out_af2.silent'
        
        else:
            if os.path.isfile(f'{working_dir}/output/{run_number}/{run_number}_input_out.silent'):
                af2_file=f'{working_dir}/output/{run_number}/{run_number}_input_out_af2.silent'
            else:
                design_number=int(re.search(r'run_\d+_design_(\d+).*',description).group(1))
                af2_file=f'{working_dir}/output/{run_number}/{run_number}_design_{design_number}_input_out_af2.silent'

    else:
        af2_file=f'{working_dir}/output/run_{run_number}/run_{run_number}_design_{gpu_number}_input_out_af2.silent'


    if os.path.isfile(af2_file):
        command = f'silentextractspecific {af2_file} ' + description + ' > pdbs_extraction.log'
        subprocess.run(command, cwd=pdbs_folder, shell=True)
    else:
        print('The AF2 file required is not being found...')
    
def multifastas(fastas_folder, output_file_name):
    '''
    Function to join all fasta files into a single multifasta file

    Input:

    fastas_folder --> Folder containing the fasta files

    output_file_name --> Name of the output multifasta file

    Output:

    Creates a fasta file with all the fasta sequences
    '''
    sequences={
            'id':[],
            'sequence':[]
            }
    for file in os.listdir(fastas_folder):
       if file.endswith('.fasta'):
           with open(os.path.join(fastas_folder, file), 'r') as infile:
               for line in infile.readlines():
                     if line.startswith('>'):
                          sequences['id'].append(line.strip('>').strip())
                     else:
                          sequences['sequence'].append(line.strip())

    # Write the sequences to a multifasta file
    with open(os.path.join(fastas_folder, output_file_name), 'w') as outfile:
        for seq_id, seq in zip(sequences['id'], sequences['sequence']):
            outfile.write(f'>{seq_id}\n{seq}\n')

def check_all_conditions(order_name, design_prefix, design_id, gg_vector, species):
    '''
    Function to check that all the conditions are met for the JohnBercow command

    Input:

    order_name --> Name of the order, used to name the output files

    design_prefix --> Prefix for the design, used to name the output files

    design_id --> Identifier for the design, used to name the output files

    gg_vector --> Golden Gate vector employed for cloning from the  

    species --> Organism for which the sequences are optimized,

    Output:

    If any of the conditions is not met, it raises an error
    '''
    stdout=''
    if order_name == '':
        stdout += 'An Order name is required\n'
        
    if design_prefix == '':
        stdout += 'A Design prefix is required\n'
    
    if design_id == '':
        stdout += 'A Design ID is required\n'
    else:
        try:
            int(design_id)
        except ValueError:
            stdout += 'Design ID must be an integer\n'
    if gg_vector == '':
        stdout += 'A Golden Gate vector is required\n'
    if species == '':
        stdout += 'An species is required\n'
    
    if stdout != '':
        stdout = 'The following conditions are not met:\n' + stdout + 'JohnBercow DNA extraction not running'
        return stdout
    return stdout
