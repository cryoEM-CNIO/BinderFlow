#!/usr/bin/env python3 

'''
Script to make a JSON file with all the binder design metadata

Input:

-- All the microrun variables

Output:

input.json: A JSON file with all the metadata
'''

import argparse
import json
import os

def main():
    # Initialize the argument parser
    parser = argparse.ArgumentParser(description="Generate JSON file with specified parameters.")

    # Add arguments for each key in the JSON file
    parser.add_argument("--input", type=str, required=True, help="Path for design output.")
    parser.add_argument("--template", type=str, required=True, help="Name of the binder.")
    parser.add_argument("--max_threads", type=str, required=True, help="Maximum number of threads.")
    parser.add_argument("--rfd_contigs", type=str, required=True, help="Chains involved in the design.")
    parser.add_argument("--rfd_hotspots", type=str, required=True, help="Target hotspot residues.")
    parser.add_argument("--rfd_ndesigns", type=int, required=True, help="Number of designs.")
    parser.add_argument("--pmp_nseqs", type=int, help="Number of final designs.")
    parser.add_argument("--pmp_relax_cycles", type=str, required=True, help="Number of relaxation cycles.")
    parser.add_argument("--partial_diff", type=str, required=True, help="Partial differential setting.")
    parser.add_argument("--noise_steps", type=str, required=True, help="Noise steps.")
    parser.add_argument("--noise_scale", type=str, required=True, help="Noise scale.")
    parser.add_argument("--ckp", type=str, required=True, help="Checkpoint file.")
    parser.add_argument("--node", type=str, required=False, default="", help="Node configuration.")
    parser.add_argument("--residues", type=str, required=True, help='Residues to fix, useful for scaffolding')
    parser.add_argument('--hits_number', type=int, required=True, help='Number of hits to design')
    # Parse the arguments
    args = parser.parse_args()

    # Create the JSON data
    json_data = {
        "input": args.input,
        "template": args.template,
        "max_threads": args.max_threads,
        "rfd_contigs": args.rfd_contigs,
        "rfd_hotspots": args.rfd_hotspots,
        "rfd_ndesigns": args.rfd_ndesigns,
        "pmp_nseqs": args.pmp_nseqs,
        "pmp_relax_cycles": args.pmp_relax_cycles,
        "partial_diff": args.partial_diff,
        "noise_steps": args.noise_steps,
        "noise_scale": args.noise_scale,
        "checkpoint": args.ckp,
        "node": args.node,
        "residues":args.residues,
        "hits_number": args.hits_number
    }

    output_path=os.getcwd()
    # Write the JSON data to a file
    output_file = f"{output_path}/input.json"
    with open(output_file, "w") as f:
        json.dump(json_data, f, indent=4)

    print(f"JSON file generated and saved as {output_file}.")

if __name__ == "__main__":
    main()
