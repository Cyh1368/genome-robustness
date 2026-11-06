# CARD DNABERT-2 Fine-Tuning Report

## Summary

DNABERT-2 was fine-tuned to classify DNA coding sequences as CARD antibiotic-resistance sequences or non-resistance controls. The complete run finished successfully on the Yale Bouchet cluster in Slurm job `32221728`.

The final test results were:

| Metric | Value |
| --- | ---: |
| Accuracy | 0.9566 |
| Macro F1 | 0.9565 |
| Precision | 0.9580 |
| Recall | 0.9566 |
| Matthews correlation | 0.9145 |
| Test loss | 0.1520 |

The model selected by validation loss is saved in `outputs/card_amr_binary_dnabert2/`, and the selected best checkpoint is `checkpoint-600`.

## Data sources and task definition

### Positive examples

Positive examples came from the current CARD automated data archive:

`https://card.mcmaster.ca/latest/data`

The pipeline reads all CARD nucleotide FASTA files whose names begin with `nucleotide_fasta_`. It extracts the stable `ARO:<number>` identifier from each FASTA header, uppercases the sequence, removes non-ACGT characters, and retains sequences from 90 through 2,500 nucleotides. Duplicate `(ARO, sequence)` pairs are removed, and no more than 20 sequences are sampled for any one ARO term using random seed 42.

### Control examples

Controls are coding sequences extracted from four NCBI RefSeq reference genomes:

- *Escherichia coli* K-12
- *Salmonella enterica* LT2
- *Pseudomonas aeruginosa* PAO1
- *Bacillus subtilis* 168

The pipeline downloads each genome FASTA and GFF3 file, extracts CDS features, reverse-complements minus-strand features, keeps sequences from 90 through 2,500 nucleotides, removes duplicates, and samples controls to match the number of positives in each split. The controls are used as the negative class. They are coding-sequence controls, not a claim that every control gene has no biological relationship to resistance.

### Splitting and counts

Splits are made by ARO term, not by individual sequence. This prevents alleles from the same CARD family from appearing in multiple splits. A seeded shuffle assigns approximately 70% of ARO terms to training, 15% to development, and 15% to testing.

| Split | CARD positives | Controls | Total rows including header |
| --- | ---: | ---: | ---: |
| Train | 4,401 | 4,401 | 8,803 |
| Dev | 943 | 943 | 1,887 |
| Test | 944 | 944 | 1,889 |

The dataset builder is [scripts/prepare_card_dataset.py](scripts/prepare_card_dataset.py).

## Environment and cluster execution

The cluster project directory is:

`/nfs/roberts/project/cpsc4710/cpsc4710_ch2499/genome-robustness`

A Python virtual environment is stored at `.venv`. The environment uses PyTorch 2.5.1 with CUDA 12.4 support and the DNABERT-2 dependencies, including Transformers 4.29.2, PEFT 0.3.0, Accelerate 0.20.3, scikit-learn 1.2.2, Biopython 1.85, and NumPy 1.24.4.

The CPU staging job installed packages and built the dataset before GPU training. The successful staging job was `32221391`. The successful GPU job was `32221728`, using one NVIDIA A5000 GPU, 12 CPUs, 96 GB requested memory, and an 8-hour time limit. It completed in 10 minutes 22 seconds. No job from this run remains active.

The batch scripts are:

- [scripts/stage_card_dnabert2.sbatch](scripts/stage_card_dnabert2.sbatch)
- [scripts/run_card_dnabert2.sbatch](scripts/run_card_dnabert2.sbatch)

## Model and training configuration

- Base model: `zhihan1996/DNABERT-2-117M`
- Task: binary sequence classification
- Tokenizer: DNABERT-2 BPE tokenizer
- Maximum tokenized length: 512
- Per-device training batch size: 8
- Per-device evaluation batch size: 16
- Gradient accumulation: 4 steps
- Effective training batch size: 32
- Learning rate: `3e-5`
- Epochs: 5
- Warmup steps: 50
- Optimizer: AdamW through the Transformers Trainer
- Mixed precision: FP16
- Evaluation and checkpointing: every 200 steps
- Maximum retained checkpoints: 2
- Random seed: 42, from the training script defaults

The run performed 1,375 optimization steps over 8,802 training examples. The lowest development loss was 0.1875168 at epoch 2.18, corresponding to `checkpoint-600`. The Trainer restored that checkpoint before the final test evaluation and final model save.

## Troubleshooting and fixes

1. The first staging attempt passed archive members to Biopython in binary mode. The dataset script now wraps extracted FASTA members in `io.TextIOWrapper`.
2. The first parser correction used an over-escaped ARO regular expression. It now uses `r"ARO:\d+"` in the Python source, which matches CARD headers such as `ARO:3002999`.
3. The initial package installation combined scikit-learn 1.2.2 with NumPy 2.0 and failed with a binary incompatibility. The staging script pins NumPy to 1.24.4.
4. DNABERT-2's optional Triton flash-attention path attempted to compile a CUDA helper and failed because the cluster Python installation does not provide `Python.h`. The GPU script forces DNABERT-2's supported PyTorch attention fallback. This changes the attention implementation for compatibility, not the model weights or task.
5. The first GPU request targeted an unavailable L40S node. The pending request was changed in place to an A5000 GPU to preserve queue age and obtain an available resource.

## Results and artifacts

The final evaluation file is:

`outputs/card_amr_binary_dnabert2/results/card_amr_binary_dnabert2/eval_results.json`

The final model directory contains `config.json`, `pytorch_model.bin`, tokenizer files, training arguments, and Trainer state. Intermediate checkpoints include `checkpoint-600` and `checkpoint-1200`; older checkpoints were removed by `save_total_limit=2`.

The training log is `logs/card-dnabert2-32221728.out`. It records the epoch-by-epoch development metrics and the final test evaluation. The corresponding error log contains no traceback or fatal error after the attention fallback was enabled.

## Reproduction

From the cluster project directory:

```bash
source .venv/bin/activate
sbatch scripts/stage_card_dnabert2.sbatch
# Submit the GPU script after staging completes successfully.
sbatch --partition=gpu --gres=gpu:a5000:1 scripts/run_card_dnabert2.sbatch
```

For a dependent submission, use `--dependency=afterok:<staging_job_id>`. Before submitting a new GPU run, compare current Slurm start estimates across compatible GPU partitions as the estimates can change.
