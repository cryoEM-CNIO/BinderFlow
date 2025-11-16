# BinderFlow

## BinderFlow in detail

In BinderFlow, we have integrated backbone design, filtering of suboptimal backbones, sequence inference and score calculations in a continuous workflow, automating the input and output handling between each step. As a result, binder design can be executed end-to-end in batches of a small number of designs per job. 

A design campaign is split into multiple, parallel BinderFlow instances, each addressing a single batch of designs and independently executed on available GPUs. Each BinderFlow instance executes the following scripts sequentially per GPU: 

- **RFD.sh**: runs RFD to produce n binder backbones. 

- **align_filtering.sh**: replaces the cropped target chain with a complete version, providing more context for the sequence assignment and scoring algorithms. To avoid subsequent waste of computational resources, it also filters out designs with steric clashes and backbones composed of a single, long helix or a single hairpin, which are difficult to produce experimentally. 

- **pMPNN.sh**: runs pMPNN to assign sequences to the binder backbone in the context of the target.  

- **scoring.sh**: predicts the binder-target complex structure using a modified version of AF2IG to obtain the per-residue predicted Local Distance Difference Test (pLDDT) and Predicted Aligned Error (PAE) scores. Then, it calculates physics-based metrics using PyRosetta to further filter and characterise the designs. 

Before finishing, every BinderFlow instance appends the scores associated with each design to a .csv file for live monitoring of the campaign and evaluates the number of designs that meet the *in silico* conditions for a binder to be considered a hit (e.g., PAE_interaction < 10, pLDDT_binder > 80). This workflow is then repeated until the desired number of hits is obtained. This architecture is fully modular, so each step can be adapted to a different software that performs similar functions, and further steps can be added or substituted as new tools become available.  


## Setting BinderFlow
<details>
    <summary> PDB preparation, Input contigs and Folder organization </summary>

Once you have [installed](../README.md), you have to modify **config.sh**. This file contains all the paths and environment names **binderflow.sh** looks at, as well as the SLURM settings. Change it accordingly to your situation, some example variables have been setted up

### Preparing input

Before running BinderFlow, is important to select the protein target for the binder generation. Since RFD backbone generation scales with the number of residues, you can use a cropped version of your target as input, being careful hydrophobic patches are not exposed in the vicinity of the targeted surface. If you think some of the context you remove can be important for subsequent sequence assignment and design scoring, you can kept a template structure inside your working folder. BinderFlow will replace the target for the template structure after RFD. 

We recommend setting the input structure as a single chain structure, as chain 'B' and with numeration starting in 1000 or above. We have notice that if binder and target residues share numeration, AF2IG scoring is worse. 

For hotspot selection, the pipeline yields better results if the hotspot region is hydrophobic and concave (but not too narrow). We recommend selecting a maximum of three closeby residues (they can be hydrophobic or not). If you select residue far away the model will tend to generate elongated structures which are more difficult to produce.

### Getting input contigs

RFD needs a special description of the target to work properly. Moreover, during this target description we have to provide the lenght of our binder which can be added as a range. We can get the contigs of the target using **contigs_map_getter.py**, which is stored in the scripts dir. The command to run will be:

```bash
conda activate watcher
python3 /path/to/contigs_map_getter.py --input /path/to/input.pdb
```
This command will print on screen the contigs of the target, e.g for `PDL1_trimmed.pdb` will be `[ B1018-1132/0 ]`. You can add the desired length as `[65-155/0 B1018-1132]`

If you want to use some initial protein structure as scaffold (named with the chain A, with the target being the chain B), you can run the same command with the scaffold-target pdb file. It will only return the contigs for the chain B. Then you can detail the chain A residues the same way, and specified whether to grow the structure from the N-ter (left), C-ter (right) or both. e.g A valid contig for `scaffold.pdb` will be `[ 15-20/AWWWW-YYYY/15-20/0 BXXXX-ZZZZ/0]`

If you want the contigs for a Partial Diffusion campaign, you must prepare the input for the run (two chains, binder A and target B, with continuous numeration) and run again the script, but this time with an extra flag

```bash
conda activate watcher
python3 /path/to/contigs_map_getter.py --input /path/to/partial_diff_input.pdb --partial_diff True
```

This command will print on screen the contigs of the target, e.g for `partial_diff_input.pdb` will be `[ VV-VV BXX-ZZZ/0 ]`


### Preparing folder organization

BinderFlow will create the following directory structure once running:

```
project_dir/
├── outputs/                                                                     # Logs generated by RFD
├── output/                                                                      # Directory with the designs
    ├── run_1/
        ├── run_1_gpu_0_design_0.pdb                                             # Initial structure generation (Gly backbone)
        ...                                         
        ├── run_1_gpu_0_design_N.pdb                                             # Initial structure generation (Gly backbone)
        ├── run_1_gpu_1_design_0.pdb                                             # Initial structure generation (Gly backbone)
        ...
        ├── run_1_gpu_1_design_N.pdb                                             # Initial structure generation (Gly backbone)
        ├── run_1_gpu_K_design_N.pdb
        ├── run_1_gpu_0_design_0_substituted.pdb                                 # Initial structure with the target substituted with the template
        ...
        ├── run_1_gpu_0_design_N_substituted.pdb                                 # Initial structure with the target substituted with the template
        ├── run_1_gpu_K_design_N_substituted.pdb                                 # Initial structure with the target substituted with the template
        ...
        ├── run_1_input_design_1_input.silent                                    # Silent input of pMPNN
        ├── run_1__design_1_input_out.silent                                     # Silent output of pMPNN and input of AF2-IG
        ├── run_1_design_1_input_out_af2.silent                                  # Silent output of AF2-IG and input of scoring.py
        ├── run_1_design_1_input_out_af2.sc                                      # Scoring of AF2-IG
        ├── pae_run_1_gpu_0_design_N_substituted_dldesign_0_cycle1_af2pred.json  #JSON with the PAE and pLDDT info of the AF2 prediction
        ...
        ├──traj/                                                                 #Directory where the trajectories are stored (for cool movies)
        └── slurm_logs/                                                          # logs of each of the jobs send to the cluster
            ├── [job_number].out                                                 # outfile of the job 
            ├── [job_number].err                                                 # errorfile of the job    
            └── [job_number]_gpuK/
                ├── rfd.out
                ├── rfd.err
                ├── aligning_filtering.out
                ├── aligning_filtering.err
                ├── pmpnn.out
                ├── pmpnn.err
                ├── scoring.out
                └── scoring.err
    ├── run_2/
        ...
    ...
    └── run_I/
        ...
├── Scoring_Stats.csv                                                            # CSV file with the scoring of all designs
├── [project_name].log                                                           # Log file of the project with all the slurm job information
└── input.json                                                                   # JSON with metadata of the run to ensure its reproducibility
```

To keep everything nice and tidy we recommend running BinderFlow from a dedicated project directory, where we store the input and template structures in a folder called input
</details>

## BinderFlow input

BinderFlow can take inputs in two different ways. The easiest way, and the one we recommend is through a JSON file that contains all the information about the run. This JSON file has the same structure as the one BinderFlow creates when it starts a run, facilitating reproducibility between campaings.

The JSON file structure is:

```json
{
    "input": "input/PDL1_trimmed.pdb", 
    "template": "input/PDL1_trimmed.pdb",
    "max_threads": "4",
    "rfd_contigs": "[ 65-155/0 B1018-1132/0 ]",
    "rfd_hotspots": "[B1054,B1068]",
    "rfd_ndesigns": 10,
    "pmp_nseqs": 1,
    "pmp_relax_cycles": 0,
    "partial_diff": "False",
    "noise_steps": 20,
    "noise_scale": 1,
    "checkpoint": "${RFD_PATH}/models/Complex_base_ckpt.pt",
    "node": "",
    "residues": "None"
}
```

The other option is through flags. There are few which are mandatory and have to be provided, while others have default values.

**The mandatory flags are:**
- `--input`: Path to the target structure
- `--template`: Path to the template structure (uncropped or less cropped version of the input)
- `--rfd_contigs`: Contigs map description of the input structure and binder length. This is optional for Partial Diffusion processes
- `--rfd_hotspots`: Hotspot residues for the binder design (The binder structure is going to be generated close to these residues)
- `--max_threads`: Number of *instances* to be executed in parallel (each *instance* is one node occupied)

**The optional flags are:**
- `--pmp_nseqs`: Number of sequences to be generated by pMPNN. **Default=1**
- `--pmp_relax_cycles`: Number of Fast Relax cycles to be performed by pMPNN. **Default=1**

`Note: Due to pMPNN incompatibilities, only one sequence can be generated if FR is used `

- `--partial_diff`: Boolean flag to determine if you want to perform normal binder generation or partial diffusion. **Default=False**
- `--noise_steps`: Number of steps to perform in the partial diffusion protocol. **Default=20**. This flag is ignored outside of partial diffusion.
- `--noise_scale`: Amount of noise to add at each step during partial diffusion. **Default=1**. This flag is ignored outside of partial diffusion.
- `--ckp`: Checkpoint path to bias the structure generation. **Default="**Complex_base_ckpt.pt"
- `--node`: Node in SLURM to which the job wants to be sent. **Default=''**
- `--hits_number`: Number of designs which pass the filtering metrics, after which the process stops. **Default=100**
- `--residues`: List of residues from the design you want to fix. It should be provided between brackets, separating residues with commas and defining 
ranges with -. If you are doing scaffolding, the template structure must have the chain A you are using as scaffold (with the same length you are using for scaffolding).For example: "[1,3,10-20]" **Default="None"**
- `--json`: Json file with all the info above. **Default="None"**

## Run BinderFlow

### Initial Generation

As mentioned in the initial [readme](../README.md), the recommended way to run BinderFlow is using the JSON file, which contains all the needed information.  The command is:

```bash
nohup /path/to/binderflow.sh --json /path/to/input.json > project_name.log 2>&1 &
```

If you prefer to use the flags, you have to run:

```bash
nohup /path/to/binderflow.sh --input /path/to/input.pdb --template /path/to/template.pdb --max_threads 2 --rfd_contigs "[65-155/0 BXXXX-YYYY]" --rfd_hotspots "[BXXXX, BXXXX] > project_name.log 2>&1 &"
```

You can add some of the optional flags from the input section if you want to modify the conditions.

### Partial Diffusion

We have already implemented in BinderFlow some known *in-silico* refinement strategies.**Partial Diffusion** is a method in which you use the same tools as in the initial generation to gently explore the nearby structural space of a previous binder. To run it, you just have to set `partial_diff` to True, whether in the JSON or in the flags.

```json
{
    "input": "input/partial_diff_input.pdb",
    "template": "input/partial_diff_input.pdb",
    "max_threads": "2",
    "rfd_contigs": "[ 117-117/0 B118-232/0 ]",
    "rfd_hotspots": " ",
    "rfd_ndesigns": 10,
    "pmp_nseqs": 2,
    "pmp_relax_cycles": "0",
    "partial_diff": "True",
    "noise_steps": "20",
    "noise_scale": "1",
    "checkpoint": "${RFD_PATH}/models/Complex_base_ckpt.pt",
    "node": "",
    "fixed_residues": "",
    "hits_number": 48
}

```
```bash
nohup /path/to/binderflow.sh --json /path/to/input.json > project_name.log 2>&1 &
```

or using the flags:

```bash
nohup /path/to/binderflow.sh --input /path/to/input.pdb --template /path/to/template.pdb --max_threads 2 --rfd_contigs "[ZZZ-ZZZ/0 BXXXX-YYYY]" --partial_diff "True" > project_name.log 2>&1 &
```

You can modify the amount of noise chainging the flags `noise_steps` and `noise_scale`

### Sequence Diversity

The other refinement strategy we have implemented is **Sequence Diversity**, in which you can explore the sequence space around an initial binder backbone structure. Sequence Diversity is implemented as an option inside *binderflow.sh*. To run it with flags, you have to run:

```bash
nohup /path/to/binderflow.sh --input /path/to/input.pdb --template /path/to/template.pdb --max_threads 2 --pmp_nseqs 100 --pmp_relax_cycles 0 --sequence_diversity "True" > project_name.log 2>&1 &
```

And to run it with a JSON file, set `sequence_diversity` to True in the file
```json
{
    "input": "input/sequence_diversity_input.pdb",
    "template": "input/sequence_diversity_input.pdb",
    "max_threads": "2",
    "rfd_contigs": " ",
    "rfd_hotspots": " ",
    "rfd_ndesigns": 10,
    "pmp_nseqs": 100,
    "pmp_relax_cycles": "0",
    "partial_diff": "False",
    "noise_steps": "20",
    "noise_scale": "1",
    "checkpoint": "${RFD_PATH}/models/Complex_base_ckpt.pt",
    "node": "",
    "fixed_residues": "",
    "hits_number": 48,
    "sequence_diversity": "True"
}

```

```bash
nohup /path/to/binderflow.sh --json /path/to/input.json > project_name.log 2>&1 &
```


## Examples

### Initial binder generation

We provide an example JSON file for initial binder generation against PDL1 in `BinderFlow/Examples/Initial_binder_generation`.

To run it, copy the folder into your working directory and from the Initial_binder_generation folder run:

```bash
nohup /path/to/binderflow.sh --json /path/to/Initial_binder_generation/input_initial_generation.json > PDL1_initial_binder.log 2>&1 &
```
In a few minutes results could be observed live running from Initial_binder_generation folder:

```bash
python3 /path/to/BFmonitor/BFmonitor.py --port 8050
```

### Partial Diffusion

We provide an example JSON file for Partial Diffusion refinement of a previous binder in `BinderFlow/Examples/Partial_diffusion`.
To run it, copy the folder into your working directory and from the Partial_diffusion folder run:

```bash
nohup /path/to/binderflow.sh --json /path/to/Partial_diffusion/input_partial_diff.json > PDL1_partial_diffusion.log 2>&1 &
``` 
In a few minutes results could be observed live running from Partial_diffusion folder:

```bash
python3 /path/to/BFmonitor/BFmonitor.py --port 8050
``` 

### Sequence Diversity

We provide an example for Sequence Diversity refinement of a previous binder in `BinderFlow/Examples/Sequence_diversity`.
To run it, copy the folder into your working directory and from the Sequence_diversity folder run:
```bash
nohup /path/to/binderflow.sh --json /path/to/Sequence_diversity/input_sequence_diversity.json > PDL1_sequence_diversity.log 2>&1 &
``` 
In a few minutes results could be observed live running from Sequence_diversity folder: 

```bash
python3 /path/to/BFmonitor/BFmonitor.py --port 8050
``` 