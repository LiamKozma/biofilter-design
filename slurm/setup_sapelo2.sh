#!/bin/bash
# One-time environment setup on UGA Sapelo2.
# Run on a login node from the project root after transferring the repo.
#
#   bash slurm/setup_sapelo2.sh
#
# Builds a venv with mpi4py compiled against the cluster's MPI so the MPI stages
# (calibration, Sobol) run across nodes. Adjust the module versions to what
# `module spider Python` and `module spider OpenMPI` report on your system.
set -euo pipefail

module load Python/3.11
module load OpenMPI

python -m venv .venv
source .venv/bin/activate
pip install --upgrade pip

# Install everything except mpi4py first...
pip install numpy scipy pandas openpyxl pyyaml emcee arviz SALib corner matplotlib schwimmbad

# ...then build mpi4py against the loaded MPI (mpicc must be on PATH).
MPICC="$(which mpicc)" pip install --no-binary=mpi4py mpi4py

echo
echo "Verifying:"
python -c "import mpi4py; print('mpi4py', mpi4py.__version__)"
python -c "import emcee, arviz, SALib, scipy; print('science stack OK')"
echo
echo "Next: regenerate tidy data, then submit jobs:"
echo "  python scripts/00_extract_data.py"
echo "  sbatch slurm/calibrate.sub"
echo "  sbatch slurm/sobol.sub"
echo "  sbatch slurm/sbc_array.sub   # then scripts/06_sbc_reduce.py"
