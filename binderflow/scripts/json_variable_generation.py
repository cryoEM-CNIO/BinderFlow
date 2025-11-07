#!/usr/bin/env python3 

'''
Script to make a JSON file with all the binder design metadata

Input:

-- All the binderflow variables

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
    parser.add_argument("--rfd_hotspots", type=str, required=False, help="Target hotspot residues.")
    parser.add_argument("--rfd_ndesigns", type=int, required=False, help="Number of designs.")
    parser.add_argument("--pmp_nseqs", type=int, help="Number of final designs.")
    parser.add_argument("--pmp_relax_cycles", type=str, required=False, help="Number of relaxation cycles.")
    parser.add_argument("--partial_diff", type=str, required=False, help="Partial differential setting.")
    parser.add_argument("--noise_steps", type=str, required=False, help="Noise steps.")
    parser.add_argument("--noise_scale", type=str, required=False, help="Noise scale.")
    parser.add_argument("--ckp", type=str, required=False, help="Checkpoint file.")
    parser.add_argument("--node", type=str, required=False, default="", help="Node configuration.")
    parser.add_argument("--residues", type=str, required=False, help='Residues to fix, useful for scaffolding')
    parser.add_argument('--hits_number', type=int, required=False, help='Number of hits to design')
    parser.add_argument('--sequence_diversity', type=str, required=False, default="False", help='Whether to enforce sequence diversity among designs')
    # Parse the arguments
    args = parser.parse_args()

    # If some arguments are None, set default values
    if args.pmp_nseqs is None:
        args.pmp_nseqs = 1
    if args.node is None:
        args.node = ""
    if args.max_threads is None:
        args.max_threads = "1"
    # Handle rfd_hotspots based on the design strategy
    if args.partial_diff == 'False' and args.sequence_diversity == 'False':
        if args.rfd_hotspots == "" or args.rfd_hotspots is None:
            raise ValueError("During normal diffusion, rfd_hotspots must be specified.")
    elif args.partial_diff == 'True' and args.sequence_diversity == 'True':
        raise ValueError("Partial diffusion and sequence diversity cannot be run at the same time.")
        # For partial diffusion or sequence diversity, ensure rfd_hotspots is a string
    else:
        if args.rfd_hotspots == "" or args.rfd_hotspots is None:
            args.rfd_hotspots = ""
        else:
            raise ValueError("During partial diffusion or sequence diversity, rfd_hotspots must be an empty string.")
    if args.hits_number is None:
        args.hits_number = 999
    if args.residues is None:
        args.residues = "None"
    if args.noise_steps is None:
        args.noise_steps = "20"
    if args.noise_scale is None:
        args.noise_scale = "1"
    if args.partial_diff is None:
        args.partial_diff = "False"
    if args.pmp_relax_cycles is None:
        args.pmp_relax_cycles = "0"
    if args.rfd_ndesigns is None:
        args.rfd_ndesigns = 10
    if args.ckp is None:
        RFD_PATH= os.environ.get("RFD_PATH")
        args.ckp = f"{RFD_PATH}/models/Complex_base_ckpt.pt"
    if args.sequence_diversity == 'True':
        args.rfd_contigs = ''
    else:
        if args.rfd_contigs is None or args.rfd_contigs == "":
            raise ValueError("RFD contigs are required.")
    # And if some of the important arguments are missing, raise an error
    if args.input is None or args.input == "":
        raise ValueError("Input path is required.")
    if args.template is None or args.template == "":
        raise ValueError("Template name is required.")



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
        "hits_number": args.hits_number,
        "sequence_diversity": args.sequence_diversity
    }

    output_path=os.getcwd()
    # Write the JSON data to a file
    output_file = f"{output_path}/input.json"
    with open(output_file, "w") as f:
        json.dump(json_data, f, indent=4)

    print(f"JSON file generated and saved as {output_file}.")

if __name__ == "__main__":
    main()
