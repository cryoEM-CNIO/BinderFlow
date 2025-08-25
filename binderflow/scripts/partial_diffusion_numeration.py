'''
Script that takes as input a pdb with the following shape:

Chain A: Binder
Chain B: Target

And change the chain B numeration so it can fits the binder (The target has a continuous numeration, which also follows the numeration of the binder)

'''


import sys
from Bio import PDB

def renumber_chain(pdb_file):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)

    chain_a = structure[0]['A']
    chain_b = structure[0]['B']

    # Get the sequence of chain A
    seq_a = ''.join([residue.get_resname() for residue in chain_a])
    
    # Get the length of chain A
    length_a = len(seq_a)/3 #Dividing by 3 to get the number of residues (each residue is represented by 3 characters)

    # Renumber chain B to match the sequence of chain A
    for i, residue in enumerate(chain_b):
        residue.id = (' ', i + length_a + 1, ' ')

    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save('renumbered_' + pdb_file)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python partial_diffusion_numeration.py <pdb_file>")
        sys.exit(1)

    pdb_file = sys.argv[1]


    renumber_chain(pdb_file)
    print(f"Renumbered chain B to match chain A in {pdb_file}. Output saved as 'renumbered_{pdb_file}'.")