#! /bin/bash

#Get script dir and load all the variables
SCRIPT_DIR="$(dirname "$(realpath "${BASH_SOURCE[0]}")")"
source $SCRIPT_DIR/config.sh 

#Checking defined enviroments are installed
echo "Checking conda environments..."
installed_envs=$(conda env list)
for env in BINDERFLOW_ENV RFD_ENV PMPENN_ENV AF2_ENV; do
    if [[ $installed_envs != *"${!env}"* ]]; then
        echo "Conda environment ${!env} not found. Please create it before running the script."
        exit 1
    fi
done
echo "All required conda environments are installed."
echo ""
echo "Checking software paths..."
# Check the paths defined indeed exist and are accessible
for var in RFD_PATH PMPNN_PATH AF2IG_PATH SILENT_PATH; do
    dir="${!var}"
    if [ ! -d "$dir" ]; then
        echo "$var directory $dir does not exist. Please check the path in config.sh."
        exit 1
    fi
done
echo "All software paths are valid."
echo ""


# Master script for deep searches in binderflow

# 1: Get all info needed
##Set defaults
partial_diff="False"
noise_steps=20
pmp_nseqs=1
rfd_ndesigns=10
pmp_relax_cycles=1
noise_scale=1
checkpoint="$RFD_PATH/models/Complex_base_ckpt.pt"
node=''
hits_number=999
fixed_residues="None"
json="None"
sequence_diversity="False"

## Parse command-line arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -j|--json) json="$2" ; shift ; break ;; # Path to the JSON file with the variables
        -i|--input) input="$2" ; shift ;;
        -t|--template) template="$2" ; shift ;;
        -m|--max_threads) max_threads="$2" ; shift  ;;    
        -c|--rfd_contigs) rfd_contigs="$2" ; shift  ;;    
        -h|--rfd_hotspots) rfd_hotspots="$2" ; shift  ;;
        -nd|--rfd_ndesigns) rfd_ndesigns="$2" ; shift  ;; #Number of designs for the RFD    
        -np|--pmp_nseqs) pmp_nseqs="$2" ; shift  ;;    
        -rc|--pmp_relax_cycles) pmp_relax_cycles="$2" ; shift  ;;   
        -pd|--partial_diff) partial_diff="$2" ; shift  ;; 
        -nst|--noise_steps) noise_steps="$2" ; shift  ;;
        -nsc|--noise_scale) noise_scale="$2" ; shift  ;;
        -ck|--ckp) checkpoint="$2" ; shift  ;; #Add the path to the checkpoint to add weight toward some fold
        -w|--node) node="$2" ; shift  ;; # Provide a specific name of a node to submit to this node with -w. If not provided it will be as usual.
        -hn|--hits_number) hits_number="$2" ; shift ;;
        -re|--residues) fixed_residues="$2" ; shift ;; #Residues index to fix, useful for scaffolding
        -sd|--sequence_diversity) sequence_diversity="$2" ; shift ;; #Whether to run sequence diversity
        *) echo "Unknown option: $1" ; exit 1 ;;
    esac
    shift  # Shift past the current argument
done
# Making a touch od an empty file to indicate this folder has a design project
touch .binder_design_project

#Getting all info from the JSON file if provided
if [ $json != "None" ]; then 
    eval $(python3 $BINDERFLOW_PATH/binderflow/scripts/input_json_reader.py $json)
    echo "Using JSON file: $json"
    echo "Input: $input"
    echo "Template: $template"
    echo "Max threads: $max_threads"
    echo "RFD Contigs: $rfd_contigs"
    echo "RFD Hotspots: $rfd_hotspots"
    echo "RFD Ndesigns: $rfd_ndesigns"
    echo "pMPNN Nseqs: $pmp_nseqs"
    echo "pMPNN Relax cycles: $pmp_relax_cycles"
    echo "Partial Diffusion: $partial_diff"
    echo "Noise steps: $noise_steps"
    echo "Noise scale: $noise_scale"
    echo "Checkpoint: $checkpoint"
    echo "Fixed residues: $fixed_residues"
    echo "Hits number: $hits_number"
    echo "Sequence diversity: $sequence_diversity"
fi

mkdir -p ./output

last_run_folder=$(ls -d "./output/run_"* 2>/dev/null | sort -V | tail -n 1 | sed 's#./output/run_##')

if [[ -n "$last_run_folder" ]]; then
    i="$last_run_folder"
else
    i=0
fi

## Contigs Getter, automatic for partial diffusion
if [ "$partial_diff" = "True" ]; then 
        rfd_contigs=$(python3 $BINDERFLOW_PATH/binderflow/scripts/contigs_map_getter.py --input "$input" --partial_diff "$partial_diff")
fi
echo "Using contigs map: $rfd_contigs"

## Prepare Folder & Variables
echo "Preparing JSON to save the run variables"

# If the following script raises an error, stop the execution
python3 $BINDERFLOW_PATH/binderflow/scripts/json_variable_generation.py --input "$input" --template "$template" \
                                                                    --max_threads "$max_threads" --rfd_contigs "$rfd_contigs" \
                                                                    --rfd_hotspots "$rfd_hotspots" --rfd_ndesigns "$rfd_ndesigns" \
                                                                    --pmp_nseqs "$pmp_nseqs" --pmp_relax_cycles "$pmp_relax_cycles" \
                                                                    --partial_diff "$partial_diff" --noise_steps "$noise_steps" \
                                                                    --noise_scale "$noise_scale" --ckp "$checkpoint" \
                                                                    --residues "$fixed_residues" --hits_number "$hits_number" \
                                                                    --sequence_diversity "$sequence_diversity" || { echo "Error: json_variable_generation.py failed. Stopping script."; exit 1; }




old_i=1
# RUN
while [ ! -f 'campaign_done' ]; do
    i=$((i+1))


    # Wait to ensure we don't exceed max_threads running jobs
    if [ "$i" -gt $max_threads ]; then
        # Count ALL jobs that are still running (from job 1 to job i-1)
        running_count=$max_threads
        while [ $running_count -ge $max_threads ]; do
            running_count=0
            running_jobs=""
            for j in $(seq 1 $((i - 1))); do
                if [ ! -e "output/run_${j}/run_${j}_done" ]; then
                    running_count=$((running_count + 1))
                    running_jobs="$running_jobs $j"
                fi
            done
            
            if [ $running_count -ge $max_threads ]; then
                echo "Currently $running_count jobs running (runs:$running_jobs). Waiting for one to complete..."
                sleep 60
            fi
        done
        echo "Jobs running: $running_count (below limit of $max_threads). Proceeding with run_${i}"
    fi
    mkdir -p ./output/run_$i/slurm_logs

    sbatch -w "$node" --nodes="$NODES" -p "$PARTITION" --open-mode=append --gres="$GRES" --exclusive --cpus-per-gpu="$CPUS_PER_GPU" -o ./output/run_$i/slurm_logs/%j.out -e ./output/run_$i/slurm_logs/%j.err \
        "$BINDERFLOW_PATH/binderflow/slurm_submit/submit_master.sh" --input "$input" --template "$template" --run "$i" --rfd_contigs "$rfd_contigs" --rfd_ndesigns "$rfd_ndesigns" \
        --pmp_nseqs "$pmp_nseqs" --pmp_relax_cycles "$pmp_relax_cycles" --partial_diff "$partial_diff" --noise_steps "$noise_steps" --noise_scale "$noise_scale" --ckp "$checkpoint" \
        --residues "$fixed_residues" --hits_number "$hits_number" --directory "$SCRIPT_DIR" --sequence_diversity "$sequence_diversity" --rfd_hotspots "$rfd_hotspots"

done

echo "Campaign finshed"
