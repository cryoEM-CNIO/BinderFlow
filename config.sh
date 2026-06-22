#Manually set the paths for every other software you have installed and environment name

export RFD_PATH="/apps/rosetta/RFDifussion"
export PMPNN_PATH="/apps/rosetta/dl_binder_design" #If you are usign the tools from nrbennett repo, it should be the same as SILENT and AF2IG
export AF2IG_PATH="/apps/rosetta/dl_binder_design" #If you are usign the tools from nrbennett repo, it should be the same as PMPNN and SILENT
export BINDERFLOW_PATH="$(dirname "$(realpath "${BASH_SOURCE[0]}")")" #Do not change this
export SILENT_PATH="/apps/rosetta/dl_binder_design" #If you are usign the tools from nrbennett repo, it should be the same as PMPNN and AF2IG

# MANUALLY SET ENVIRONMENTS NAMES (example names set)

export RFD_ENV="SE3nv4090"
export PMPNN_ENV="dl_binder_design"
export AF2_ENV="af2_binder_design"
export BINDERFLOW_ENV="watcher"


#MANUALLY SET SBATCH CONFIGURATIONS (examples in place)

export NODES=1
export PARTITION=RFD
export CPUS_PER_GPU=12
export GRES=gpu:1

source "$(conda info --base)/etc/profile.d/conda.sh" # ACTIVATE CONDA COMMAND, CHANGE FOR YOUR COMMAND IF NEEDED
