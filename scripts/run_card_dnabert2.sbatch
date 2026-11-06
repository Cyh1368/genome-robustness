#!/bin/bash
#SBATCH --job-name=card-dnabert2
#SBATCH --partition=gpu_h100
#SBATCH --gres=gpu:h100:1
#SBATCH --cpus-per-task=12
#SBATCH --mem=96G
#SBATCH --time=08:00:00
#SBATCH --output=/nfs/roberts/project/cpsc4710/cpsc4710_ch2499/genome-robustness/logs/card-dnabert2-%j.out
#SBATCH --error=/nfs/roberts/project/cpsc4710/cpsc4710_ch2499/genome-robustness/logs/card-dnabert2-%j.err

set -euo pipefail
ROOT=/nfs/roberts/project/cpsc4710/cpsc4710_ch2499/genome-robustness
cd "$ROOT"
source .venv/bin/activate
export HF_HOME="$ROOT/.cache/huggingface"
export TRANSFORMERS_CACHE="$HF_HOME"
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
mkdir -p "$HF_HOME" data/card_amr_binary outputs logs

python -c 'import torch, transformers, peft, Bio, sklearn, accelerate'
test -s data/card_amr_binary/train.csv

# The cluster Python lacks Python.h, so Triton's optional flash-attention JIT
# cannot build. DNABERT-2 has a supported PyTorch attention fallback.
for layer_file in "$HF_HOME"/modules/transformers_modules/zhihan1996/DNABERT-2-117M/*/bert_layers.py; do
  if [ -f "$layer_file" ]; then
    sed -i '/^logger = logging.getLogger/a flash_attn_qkvpacked_func = None' "$layer_file"
  fi
done

cd DNABERT_2/finetune
python train.py \
  --model_name_or_path zhihan1996/DNABERT-2-117M \
  --data_path "$ROOT/data/card_amr_binary" \
  --kmer -1 \
  --run_name card_amr_binary_dnabert2 \
  --model_max_length 512 \
  --per_device_train_batch_size 8 \
  --per_device_eval_batch_size 16 \
  --gradient_accumulation_steps 4 \
  --learning_rate 3e-5 \
  --num_train_epochs 5 \
  --fp16 \
  --save_steps 200 \
  --eval_steps 200 \
  --evaluation_strategy steps \
  --warmup_steps 50 \
  --logging_steps 25 \
  --save_total_limit 2 \
  --save_model True \
  --output_dir "$ROOT/outputs/card_amr_binary_dnabert2" \
  --log_level info \
  --find_unused_parameters False
