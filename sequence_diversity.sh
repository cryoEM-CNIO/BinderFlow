#!/bin/bash

#Get script dir and load all the variables
SCRIPT_DIR="$(dirname "$(realpath "${BASH_SOURCE[0]}")")"
source $SCRIPT_DIR/config.sh 
conda activate $BINDERFLOW_ENV 
# Load defaults 
fr=1
nseqs=1
fixed="None"
node=""

while [[ $# -gt 0 ]]; do
    key="$1"

    case $key in
        --input) input="$2" ; shift ;;
        --max_threads) threads="$2" ; shift ;;
        --total_nseqs) max="$2" ; shift ;;
        --residues) residues="$2" ; shift ;;
        --batch_nseqs) nseqs="$2" ; shift ;;
        --relax_cycles) fr="$2" ; shift ;;
        --node) node="$2" ; shift ;;
        *)
        esac
        shift
done

mkdir -p "output"

last_run_folder=$(ls -d "./output/run_"* 2>/dev/null | sort -V | tail -n 1 | sed 's#./output/run_##')

if [[ -n "$last_run_folder" ]]; then
    i="$last_run_folder"
else
    i=0
fi

## Fix residues

if [ $fixed != "None" ]; then 
    python3 "$BINDERFLOW_PATH/scripts/fixing_residues.py" --residues "$residues" --pdb_input "$input"
fi

## create silent file
"$SILENT_PATH/include/silent_tools/silentfrompdbs"  "$input" > "initial_input.silent"


#while loop to generate diversity
while true;do
    i=$((i+1))
    if [ "$i" -lt "$threads" ]; then
        previous=0
    else
        previous=$((i - $threads))
    fi
    
    # Variable generation and definition

    mkdir -p "output/run_${i}/slurm_logs"

    waitfor="output/run_${previous}/run_${previous}_done"

    if [ "$previous" -ne 0 ]; then 
        
        while [ ! -e "$waitfor" ]; do 
            echo "Waiting for previous job to complete: $waitfor"
            sleep 60
        done
    
    fi

    sbatch -w "$node" --nodes="$NODES" -p "$PARTITION" --open-mode=append --gres="$GRES" --exclusive --cpus-per-gpu="$CPUS_PER_GPU" -o ./output/run_$i/slurm_logs/%j.out -e ./output/run_$i/slurm_logs/%j.err \
            $BINDERFLOW_PATH/binderflow/slurm_submit/submit_sequence_diversity.sh --run "$i" --nseqs "$nseqs" --fr "$fr" --directory "$SCRIPT_DIR"

    total_seqs_generated=$(($i*4*$nseqs)) #This is patatero, we have to change it
    if [ $total_seqs_generated -gt $max ]; then 
        break 
    fi  
done
counter=$((counter+1))

echo " #######################"
echo "- DIVERSITY GENERATION CONCLUDED"
echo "- ${total_seqs_generated} sequences generated in total"
echo " #######################"