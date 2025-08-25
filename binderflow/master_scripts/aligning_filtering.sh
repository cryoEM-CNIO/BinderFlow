#!/bin/bash
# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    key="$1"

    case $key in
        --template) template="$2" ; shift ;;
        --run) run="$2" ; shift ;;
        --t) t="$2" ; shift ;;
        --residues) residues="$2" ; shift ;;
        --directory) directory="$2" ; shift ;;
        --partial_diff) partial_diff="$2"; shift ;;
        *) echo "Unknown option: $1" ; exit 1 ;;
    esac
    shift 
done

#Load all variables
source $directory/config.sh
conda activate $BINDERFLOW_ENV



# Variables declaring
silent_output="output/run_${run}/run_${run}_design_${t}_input.silent"



# Display the parsed values
machine=`hostname`
echo "Current machine $machine"
# echo "Template aligning and chain substitution using Biopython for folder $input_dir"
# echo "Silent output for next step will be saved as $silent_output at $input_dir"

# Run!
#run align

python3 $BINDERFLOW_PATH/binderflow/scripts/biopython_align.py --template $template --chain "B" --run "$run" --t "$t"  --residues "$residues" --partial_diff "$partial_diff"
#run silent to split in 4 groups each for each gpu
echo "Creating silent files"
/apps/rosetta/dl_binder_design/include/silent_tools/silentfrompdbsparallel "output/run_${run}/run_${run}_gpu_${t}_design_*_substituted.pdb" > "$silent_output"

echo "done"
