#!/bin/bash

while [[ $# -gt 0 ]]; do
    key="$1"

    case $key in
        --directory) directory="$2" ; shift ;; 
        --run) run="$2" ; shift ;;
        --n_seqs) n_seqs="$2" ; shift ;;
        --relax_cycles) relax_cycles="$2" ; shift ;;
        --t) t="$2" ; shift ;;
        *) echo "Unknown option: $1" ; exit 1 ;;
    esac
    shift # Shift past the current argument
done

#Load all variables
source $directory/config.sh
conda activate $PMPNN_ENV

machine=`hostname`
echo "Current machine $machine"
input_silent="output/run_${run}/run_${run}_design_${t}_input.silent"
silent_out=`echo "$input_silent" | sed 's#.silent#_out.silent#'`
silent_point=`echo "$input_silent" | sed 's#.silent#_out.point#'`
fixed_residues_path="output/run_$i/fixed_residues_$t.json"
echo "Running pMPNN on $input_silent"

python3 -u $PMPNN_PATH/mpnn_fr/dl_interface_design.py -silent "$input_silent" -checkpoint_path "$PMPNN_PATH/mpnn_fr/ProteinMPNN/vanilla_model_weights/v_48_030.pt" -outsilent "$silent_out" -relax_cycles "$relax_cycles" -seqs_per_struct "$n_seqs" -checkpoint_name "$silent_point" #requires using our own modified version of mpnn, in which temp.pdb has a different name

echo "done"

