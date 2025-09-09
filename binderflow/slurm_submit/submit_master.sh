#!/bin/bash

## Parse command-line arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -i|--input) input="$2" ; shift ;;
        -t|--template) template="$2" ; shift ;;
        -r|--run) run="$2" ; shift ;;
        -c|--rfd_contigs) rfd_contigs="$2" ; shift  ;;    
        -h|--rfd_hotspots) rfd_hotspots="$2" ; shift  ;;    
        -nd|--rfd_ndesigns) rfd_ndesigns="$2" ; shift  ;;    
        -np|--pmp_nseqs) pmp_nseqs="$2" ; shift  ;;    
        -rc|--pmp_relax_cycles) pmp_relax_cycles="$2" ; shift  ;;   
        -pd|--partial_diff) partial_diff="$2" ; shift  ;; 
        -nst|--noise_steps) noise_steps="$2" ; shift  ;;
        -nsc|--noise_scale) noise_scale="$2" ; shift  ;;
        -ck|--ckp) ckp="$2" ; shift  ;; #Add the path to the checkpoint to add weight toward some fold
        -w|--node) node="$2" ; shift  ;; # Provide a specific name of a node to submit to this node with -w. If not provided it will be as usual.
        -hn|--hits_number) hits_number="$2" ; shift ;;
        -re|--residues) residues="$2" ; shift ;; # Residues to fix, useful for scaffolding
        -d|--directory) directory="$2" ; shift ;; # Directory where the script is located
        *) echo "Unknown option: $1" ; exit 1 ;;
    esac
    shift  # Shift past the current argument
done

source $directory/config.sh 


# Get available GPUs from SLURM
GPUS_AVAILABLE=$(nvidia-smi --query-gpu=index --format=csv,noheader | tr '\n' ' ')
echo "GPUs available: $GPUS_AVAILABLE"

t=1


for GPU_ID in $GPUS_AVAILABLE; do
    echo "Using $GPU_ID"
    (
        export CUDA_VISIBLE_DEVICES=$GPU_ID
        LOG_DIR="output/run_${run}/slurm_logs/${SLURM_JOB_ID}_gpu${GPU_ID}"
        mkdir -p "$LOG_DIR"

        # ----------------------------------------
        # 1. RFD
        # ----------------------------------------
        bash $BINDERFLOW_PATH/binderflow/master_scripts/rfd.sh \
            --run "$run" --t "$GPU_ID" \
            --input_pdb "$input" --contigmap_descriptor "$rfd_contigs" \
            --designs_n "$rfd_ndesigns" --noise_steps "$noise_steps" \
            --noise_scale "$noise_scale" --ckp "$ckp" \
            --partial_diff "$partial_diff" --residues "$residues" --hotspots_descriptor "$rfd_hotspots" --directory "$directory"> "$LOG_DIR/rfd.out" 2> "$LOG_DIR/rfd.err"
        wait


        # --------------------------------------------
        # 2. Aligning + filtering
        # --------------------------------------------
        bash $BINDERFLOW_PATH/binderflow/master_scripts/aligning_filtering.sh --template "$template" --run "$run" --t "$GPU_ID" --residues "$residues" --directory "$directory" > "$LOG_DIR/aligning_filtering.out" 2> "$LOG_DIR/aligning_filtering.err"
        wait


        # --------------------------------------------
        # 3. pMPNN
        # --------------------------------------------
        bash $BINDERFLOW_PATH/binderflow/master_scripts/pmpnn.sh --run "$run" --t "$GPU_ID" --n_seqs "$pmp_nseqs" --relax_cycles "$pmp_relax_cycles" --directory "$directory" > "$LOG_DIR/pmpnn.out" 2> "$LOG_DIR/pmpnn.err"
        wait


        # --------------------------------------------
        # 4. Scoring
        # --------------------------------------------
        bash $BINDERFLOW_PATH/binderflow/master_scripts/scoring.sh --run "$run" --t "$GPU_ID" --directory "$directory" > "$LOG_DIR/scoring.out" 2> "$LOG_DIR/scoring.err"
        wait

    ) &
    ((t=t+1))
done
wait
# --------------------------------------------
# 5 Finish binderflow
# --------------------------------------------

bash $BINDERFLOW_PATH/binderflow/master_scripts/ending.sh --number "$hits_number" --run "$run" --directory "$directory"