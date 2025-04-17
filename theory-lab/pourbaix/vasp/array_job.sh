#!/bin/bash
#SBATCH --job-name=vasp
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --threads-per-core=1
#SBATCH --time=1-00:00:00

echo "$(date "+%Y-%m-%d %H:%M:%S") Job started"

CONFIG_ID=$(printf "%05d" "$SLURM_ARRAY_TASK_ID")
CONFIG_DIR="$PWD/configs/$CONFIG_ID"

python generate_vasp_inputs.py --index=$SLURM_ARRAY_TASK_ID --output_dir=$CONFIG_DIR

SCRATCH_BASE="/scratch/$USER"
SCRATCH_JOBDIR="$SCRATCH_BASE/vasp_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p $SCRATCH_JOBDIR

# Copy input files to scratch
rsync -a "$CONFIG_DIR/" "$SCRATCH_JOBDIR/"
pushd $SCRATCH_JOBDIR
srun /home/lucas/VASP.5.4.4/vasp.5.4.4/bin/vasp_std
popd 

# Copy results back
rsync -a "$SCRATCH_JOBDIR/" "$CONFIG_DIR/"
rm -rf "$SCRATCH_JOBDIR"

echo "$(date "+%Y-%m-%d %H:%M:%S") Job finished"
