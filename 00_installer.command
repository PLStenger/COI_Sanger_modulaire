#!/bin/bash
set -euo pipefail
command -v conda >/dev/null || { echo 'Installer Miniforge/Conda puis relancer.'; exit 1; }
if conda env list | awk '{print $1}' | grep -qx coi_sanger; then
 conda install -y -n coi_sanger -c conda-forge -c bioconda --strict-channel-priority 'python=3.11' 'biopython=1.85' 'matplotlib>=3.8' numpy mafft iqtree
else
 conda create -y -n coi_sanger -c conda-forge -c bioconda --strict-channel-priority 'python=3.11' 'biopython=1.85' 'matplotlib>=3.8' numpy mafft iqtree
fi
conda run -n coi_sanger python -c 'import Bio, matplotlib; print("Biopython",Bio.__version__,"Matplotlib",matplotlib.__version__)'
echo 'Environnement prêt. Utiliser conda run --no-capture-output -n coi_sanger ...'
