#!/bin/bash
export HYDRA_FULL_ERROR=1


# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --run) run="$2" ; shift ;;
        --t) t="$2" ; shift ;; 
        -ip|--input_pdb) input_pdb="$2" ; shift ;;
        -cd|--contigmap_descriptor) contigmap_descriptor="$2" ; shift ;;
        -hd|--hotspots_descriptor) hotspots_descriptor="$2" ; shift ;;
        -dn|--designs_n) designs_n="$2" ; shift ;;
        -c|--ckp) ckp="$2" ; shift ;;
        -pd|--partial_diff) partial_diff="$2" ; shift ;;
        -nst|--noise_steps) noise_steps="$2" ; shift ;;
        -nsc|--noise_scale) noise_scale="$2" ; shift ;;
        --r|--residues) residues="$2"; shift ;; 
        --directory) directory="$2" ; shift ;; 

        *) echo "Unknown option: $1" ; exit 1 ;;
    esac
    shift 
done

#Load all variables
source $directory/config.sh
conda activate $RFD_ENV

design_startnum=0

#Define variables
output_prefix="output/run_${run}/run_${run}_gpu_${t}_design"
# Display the parsed values
machine=`hostname`
echo "Current machine $machine"
echo "You chose the following input values:"
echo "  - output prefix: $output_prefix"
echo "  - input PDB: $input_pdb"
echo "  - contigmap descriptor: $contigmap_descriptor"
echo "  - hotspots descriptor: $hotspots_descriptor"
echo "  - number of designs: $designs_n"
echo "  - checkpoint used: $ckp"
echo "  - Partial Diffusion: $partial_diff"
echo "  - noise scale (only if PD): $noise_scale"
echo "  - noise_steps (only if PD): $noise_steps"
echo "  - residues fixed: $residues"

# Run!
if [ $partial_diff = "True" ]; then
    if [ "$residues" = "None" ]; then 
        $RFD_PATH/scripts/run_inference.py inference.design_startnum="$design_startnum" inference.output_prefix="$output_prefix" inference.input_pdb="$input_pdb" contigmap.contigs="$contigmap_descriptor"  inference.num_designs="$designs_n" diffuser.partial_T="$noise_steps" denoiser.noise_scale_ca="$noise_scale" denoiser.noise_scale_frame="$noise_scale"
    else
        $RFD_PATH/scripts/run_inference.py inference.design_startnum="$design_startnum" inference.output_prefix="$output_prefix" inference.input_pdb="$input_pdb" contigmap.contigs="$contigmap_descriptor"  inference.num_designs="$designs_n" diffuser.partial_T="$noise_steps" denoiser.noise_scale_ca="$noise_scale" denoiser.noise_scale_frame="$noise_scale" contigmap.provide_seq="$residues"
    fi
else
    if [ -z "${hotspots_descriptor+x}" ]; then
        $RFD_PATH/scripts/run_inference.py inference.design_startnum="$design_startnum" inference.output_prefix="$output_prefix" inference.input_pdb="$input_pdb" contigmap.contigs="$contigmap_descriptor" inference.num_designs=$designs_n denoiser.noise_scale_ca=0 denoiser.noise_scale_frame=0 inference.ckpt_override_path="$ckp"
    else
        $RFD_PATH/scripts/run_inference.py inference.design_startnum="$design_startnum" inference.output_prefix="$output_prefix" inference.input_pdb="$input_pdb" contigmap.contigs="$contigmap_descriptor" ppi.hotspot_res="$hotspots_descriptor" inference.num_designs=$designs_n denoiser.noise_scale_ca=0 denoiser.noise_scale_frame=0 inference.ckpt_override_path="$ckp"
    fi
fi
echo "done"