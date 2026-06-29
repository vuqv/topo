#!/bin/bash
# ===========================================================================
# O'Brien Continuous Synthesis Protocol (v6) -- full-length 4c5c reproduction.
#   Submit:   sbatch submit_csp.slurm
#   Resume:   just sbatch again -- the script auto-detects prior progress and
#             flips the control file to restart=1 (v6 resumes from the last
#             completed nascent length).
#
# >>> EDIT the two lines marked TODO for your allocation (GPU partition/account).
# ===========================================================================
#SBATCH --job-name=csp_4c5c
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1                 # one GPU (v6 has use_gpu = 1)
#SBATCH --mem=16G
#SBATCH --time=96:00:00              # full length is long; resubmit to resume
#SBATCH --partition=standard         # TODO: your GPU partition (e.g. sla-prio / gpu)
#SBATCH --account=epo2_cr_default         # TODO: your GPU allocation/account

set -euo pipefail

# --- run from the directory that holds this script / the inputs --------------
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"

CTRL=cont_synth_full.cntrl

# --- environment: conda env with openmm + parmed ----------------------------
source /storage/group/epo2/default/qzv5006/miniconda3/etc/profile.d/conda.sh
conda activate bioenv

# v6 shells out to parse_cg_prm.py and ribosome_traffic by bare name, so the
# tutorial folder (which holds both) must be on PATH.
export PATH="$PWD:$PATH"

echo "=== job $SLURM_JOB_ID on $(hostname) @ $(date) ==="
nvidia-smi -L || true
python -c "import openmm; print('OpenMM', openmm.__version__)"

# --- auto-resume: if a prior run left output/, switch the control to restart=1
# (v6 scans output/*.out for the last finished length and continues from there).
if ls output/*.out >/dev/null 2>&1; then
    echo "Found prior output/ -> resuming (restart=1)."
    sed -i 's/^restart = .*/restart = 1/' "$CTRL"
else
    echo "Fresh run (restart=0)."
    sed -i 's/^restart = .*/restart = 0/' "$CTRL"
fi
grep -E '^scale_factor|^total_nascent_chain_length|^ribosome_traffic|^restart|^use_gpu' "$CTRL"

# --- run --------------------------------------------------------------------
echo "=== starting v6 @ $(date) ==="
python continuous_synthesis_v6.py -f "$CTRL"
echo "=== finished @ $(date) (exit $?) ==="
