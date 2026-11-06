#!/bin/bash
#SBATCH --job-name=stage-card-data
#SBATCH --partition=day
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=/nfs/roberts/project/cpsc4710/cpsc4710_ch2499/genome-robustness/logs/stage-card-%j.out
#SBATCH --error=/nfs/roberts/project/cpsc4710/cpsc4710_ch2499/genome-robustness/logs/stage-card-%j.err

set -euo pipefail
ROOT=/nfs/roberts/project/cpsc4710/cpsc4710_ch2499/genome-robustness
cd "$ROOT"
source .venv/bin/activate
export HF_HOME="$ROOT/.cache/huggingface"
mkdir -p "$HF_HOME" data/card_amr_binary logs
python -c 'import torch' 2>/dev/null || python -m pip install --no-cache-dir torch==2.5.1 --index-url https://download.pytorch.org/whl/cu124
python -c 'import transformers, peft, Bio, sklearn, accelerate' 2>/dev/null || python -m pip install --no-cache-dir numpy==1.24.4 einops==0.6.1 transformers==4.29.2 peft==0.3.0 omegaconf==2.3.0 evaluate==0.4.0 accelerate==0.20.3 scikit-learn==1.2.2 biopython
python scripts/prepare_card_dataset.py --raw-dir data/card_raw --output-dir data/card_amr_binary
python - <<'PY'
import torch
assert torch.cuda.is_available() is False
for split in ("train", "dev", "test"):
    with open(f"data/card_amr_binary/{split}.csv") as handle:
        assert sum(1 for _ in handle) > 1, split
print("Staging complete")
PY
