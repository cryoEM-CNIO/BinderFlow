#!/bin/bash

#Default option
pkg_manager="conda"


while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -p|--pkg) pkg_manager="$2" ; shift ;;
        *) echo "Unknown option: $1" ; exit 1 ;;
    esac
done

#Output selected package manager

echo -e "Package manager selected: $pkg_manager"

# set paths needed for installation and check for conda installation
install_dir=$(pwd)
CONDA_BASE=$(conda info --base 2>/dev/null) || { echo -e "Error: conda is not installed or cannot be initialised."; exit 1; }
echo -e "Conda is installed at: $CONDA_BASE"

### Start binderflow install
echo -e "Installing BinderFlow environment\n"
$pkg_manager create --name binderflow python=3.10 -y || { echo -e "Error: Failed to create binderflow conda environment"; exit 1; }
conda env list | grep -w 'binderflow' >/dev/null 2>&1 || { echo -e "Error: Conda environment 'binderflow' does not exist after creation."; exit 1; }

# Load newly created binderflow environment
echo -e "Loading binderflow environment\n"
source ${CONDA_BASE}/bin/activate ${CONDA_BASE}/envs/binderflow || { echo -e "Error: Failed to activate the binderflow environment."; exit 1; }
[ "$CONDA_DEFAULT_ENV" = "binderflow" ] || { echo -e "Error: The binderflow environment is not active."; exit 1; }
echo -e "binderflow environment activated at ${CONDA_BASE}/envs/binderflow"

# install required conda packages
echo -e "Instaling conda requirements\n"
$pkg_manager install pyrosetta --channel https://conda.graylab.jhu.edu -y  || { echo -e "Error: Failed to install conda packages."; exit 1; }
pip install CodonTransformer
pip install dash
pip install dash-bio
pip install biopython
pip install dnachisel
pip install openpyxl

# make sure all required packages were installed
required_packages=(pyrosetta codontransformer dash dash-bio biopython dnachisel openpyxl)
missing_packages=()

# Check each package
for pkg in "${required_packages[@]}"; do
    conda list "$pkg" | grep -w "$pkg" >/dev/null 2>&1 || missing_packages+=("$pkg")
done

# If any packages are missing, output error and exit
if [ ${#missing_packages[@]} -ne 0 ]; then
    echo -e "Error: The following packages are missing from the environment:"
    for pkg in "${missing_packages[@]}"; do
        echo -e " - $pkg"
    done
    exit 1
fi

echo -e "Successfully finished binderflow enviroment installation!\n"
echo -e "Activate environment using command: \"$pkg_manager activate binderflow\""
