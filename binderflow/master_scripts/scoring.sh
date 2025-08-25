#!/bin/bash

while [[ $# -gt 0 ]]; do
    key="$1"

    case $key in
        --t) t="$2" ; shift ;;
        --directory) directory="$2" ; shift ;; 
        --run) run="$2" ; shift ;;
        *) echo "Unknown option: $1" ; exit 1 ;;
    esac
    shift
done

#Load all variables
source $directory/config.sh
conda activate $AF2_ENV

machine=`hostname`
echo "Current machine $machine"

input_silent="output/run_${run}/run_${run}_design_${t}_input.silent"
af2_in=`echo "$input_silent" | sed 's#\.silent#_out.silent#'`
af2_out=`echo "$input_silent" | sed 's#\.silent#_out_af2.silent#'`
af2_score=`echo "$input_silent" | sed 's#\.silent#_out_af2.sc#'`
af2_point=`echo "$input_silent" | sed 's#\.silent#_out_af2.point#'`
af2_json=$(echo "$input_silent" | sed -E 's#run_[0-9]+_design_.*_input.*\.silent#pae#')
input_scoring=`echo "$input_silent" | sed s'#\.silent#_out_af2.silent#'`

echo "Running AF2 on $af2_in"

python3 -u $AF2IG_PATH/af2_initial_guess/predict.py -silent "$af2_in" -outsilent "$af2_out" -scorefilename "$af2_score"  -checkpoint_name "$af2_point" -jsonfilename "$af2_json"

conda deactivate
conda activate $BINDERFLOW_ENV

echo "Running_scoring on $input_scoring"
python3 $BINDERFLOW_PATH/binderflow/scripts/scoring_tools.py --silent "$input_scoring" --run_number "$run"

echo "done"

