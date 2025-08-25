#!/usr/bin/env python3 

import argparse
from biopython_align import add_fixed_residues

if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument("--fixed", required=False, default=None, help='List of residues of the binder to fix (Useful for scaffolding)')
    parser.add_argument("--pdb_input", required=True, help='Path to the pdb file')
    args=parser.parse_args()

    add_fixed_residues(output_path=args.pdb_input, residues=args.fixed)
