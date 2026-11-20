# Source in every job: paths and modules for UMich Great Lakes. Edit ACCOUNT once.
export ACCOUNT=${ACCOUNT:-bioinf590f26_class}      # EDIT: your SLURM account (my.arc.umich.edu)
export SCRATCH=${SCRATCH:-/scratch/${ACCOUNT}_root/${ACCOUNT}/$USER}
export REPO=${REPO:-$HOME/BIOINF-590-Final-Team-Project}
mkdir -p "$SCRATCH"/{raw,harmonised,models,results,logs}
module load python3.11-anaconda/2024.02 2>/dev/null || module load python 2>/dev/null
source activate bioinf590 2>/dev/null || conda activate bioinf590
export PYTHONUNBUFFERED=1 HF_HOME=$SCRATCH/models/hf
