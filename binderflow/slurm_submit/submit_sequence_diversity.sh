#!/bin/bash

# Display the parsed values
hits_number=10000
## Parse command-line arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -r|--run) run="$2" ; shift ;;
        -ns|--nseqs) pmp_nseqs="$2" ; shift  ;;    
        -fr|--fr) pmp_relax_cycles="$2" ; shift  ;;   
        -hn|--hits_number) hits_number="$2" ; shift ;; 
        -d|--directory) directory="$2" ; shift ;;
        *) echo "Unknown option: $1" ; exit 1 ;;
    esac
    shift  # Shift past the current argument
done

#Load all variables
source $directory/config.sh
conda activate $MICRORUN_ENV 
# Get available GPUs from SLURM
GPUS_AVAILABLE=$(nvidia-smi --query-gpu=index --format=csv,noheader | tr '\n' ' ')
echo "GPUs available: $GPUS_AVAILABLE"

t=1



for GPU_ID in $GPUS_AVAILABLE; do
    echo $GPU_ID
    (
        export CUDA_VISIBLE_DEVICES=$GPU_ID
        LOG_DIR="output/run_$run/slurm_logs/${SLURM_JOB_ID}_gpu${GPU_ID}"
        mkdir -p "$LOG_DIR"

        # --------------------------------------------
        # 1 Generate the silent file
        # --------------------------------------------

        $SILENT_PATH/include/silent_tools/silentrename initial_input.silent "run_${run}_design_${t}" > "output/run_${run}/run_${run}_design_${t}_input.silent" 
        wait
        # --------------------------------------------
        # 2 pMPNN
        # --------------------------------------------

        bash $MICRORUN_PATH/microrun/master_scripts/pmpnn.sh --run "$run" --t "$t" --n_seqs "$pmp_nseqs" --relax_cycles "$pmp_relax_cycles" --directory "$directory" > "$LOG_DIR/pmpnn.out" 2> "$LOG_DIR/pmpnn.err"
        wait
        # --------------------------------------------
        # 3 Scoring(AF2IG + PyRosetta)
        # --------------------------------------------

        bash $MICRORUN_PATH/microrun/master_scripts/scoring.sh --run "$run" --t "$t" --directory "$directory" > "$LOG_DIR/scoring.out" 2> "$LOG_DIR/scoring.err"
        wait
    ) &
    ((t=t+1))
done
wait
# --------------------------------------------
# 4 Finish Microrun
# --------------------------------------------

bash $MICRORUN_PATH/microrun/master_scripts/ending.sh --number "$hits_number" --run "$run"