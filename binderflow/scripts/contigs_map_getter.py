#!/usr/bin/env python3

'''
Script to calculate the contigs.contigsmap that RFD needs.

Input:

-i|--input : pdb file path of the input
-pd|--partial_diff :  Is the contigs for partial diffusion or for normal binder generation ? (True or False)

Output:

Prints in the screen the groupped residue number in the contigsmap format

'''
import argparse
import re

def extract_unique_residue_numbers(pdb_file, partial_diff):
    unique_residue_numbers_B = set()
    unique_residue_numbers_A = set()
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                chain=line.split()[4][0]
                if partial_diff==True:
                    if chain=='A':
                        if not bool(re.search(r'\d', line.split()[4])): 
                            residue_number = int(line.split()[5])
                        else:
                            residue_number = int(line.split()[4][1:])
                        unique_residue_numbers_A.add(residue_number)
                if chain == 'B':
                    if not bool(re.search(r'\d', line.split()[4])): 
                        residue_number = int(line.split()[5])
                    else:
                        residue_number = int(line.split()[4][1:])
                    unique_residue_numbers_B.add(residue_number)
    #Disclaimer if the numeration is wrong for partial diffusion
    if partial_diff == True:
        if  sorted(list(unique_residue_numbers_B))[0] != sorted(list(unique_residue_numbers_A))[-1] +1:
            print('This numeration is not valid for partial diffusion. Consider modifying it using Pymol or Biopython') 
    return sorted(list(unique_residue_numbers_A)), sorted(list(unique_residue_numbers_B))

def group_consecutive_numbers(numbers,chain):
    if chain == 'A':
        range_A=[]
        if numbers:
            end=len(numbers)
            range_A.append(f'{end}-{end}/0 ')
        return ''.join(range_A)
    if chain == 'B':
        range_B=[]
        start=end=numbers[0]

        for number in numbers[1:]:
            if number == end+1:
                end=number
            else:
                range_B.append(f'B{start}-{end}/')
                start=end=number
    if start == end:
        range_B.append(str(start))
    else:
        range_B.append(f"B{start}-{end}/")
        return ''.join(range_B)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", "-i", help = "cropped PDB to calculate contigs")
    parser.add_argument("--partial_diff", "-pd", type=bool, default=False, required=False, help='Do you need the output for partial diffusion?' )

    args, unknown = parser.parse_known_args()
    partial_diff=args.partial_diff
    if not args.input:
        print('Input pdb is missing. Run as: python3 contigs_map_getter.py -i <pdbfile>')
    else:
        pdb_file = str(args.input)

    unique_residue_numbers_A, unique_residue_numbers_B = extract_unique_residue_numbers(pdb_file, partial_diff)
    grouped_residue_string_A = group_consecutive_numbers(unique_residue_numbers_A, 'A')
    grouped_residue_string_B = group_consecutive_numbers(unique_residue_numbers_B, 'B')
    print(f'[ {grouped_residue_string_A}{grouped_residue_string_B}0 ]')


    # unique_residue_string = ", ".join(map(str, unique_residue_numbers))
    # print(f"Unique Residue Numbers: {unique_residue_string}")
