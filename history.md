# History.md

=====

[2026-11-06 16:59:33]
Reported the completed CARD DNABERT-2 fine-tuning run and verified its saved model and test metrics.

- Fine-tuning job 32221728 completed successfully with exit code 0 in 10 minutes 22 seconds on one A5000 GPU.
- The run trained for five epochs over 8,802 examples and completed 1,375 optimizer steps.
- The selected best checkpoint was checkpoint-600, based on the lowest validation loss of 0.1875.
- Final test results were accuracy 0.9566, macro F1 0.9565, precision 0.9580, recall 0.9566, and Matthews correlation 0.9145.
- The final model and evaluation JSON are saved under outputs/card_amr_binary_dnabert2; no active job from this run remains.
=====

[2026-11-06 17:02:23]
Saved the detailed CARD DNABERT-2 methodology and results locally and committed the report artifacts.

- Added CARD_DNABERT2_REPORT.md with data provenance, split rules, training settings, troubleshooting, and reproducibility instructions.
- Added a local JSON copy of the final test evaluation metrics.
- Added the dataset preparation and Slurm batch scripts to the repository.
- Added Python cache exclusions and verified JSON parsing and Git whitespace checks.
- Committed all changes as 70ffbe3 with message "document CARD DNABERT2 fine-tuning".
=====

[2026-11-05 16:58:12]
Continued the cluster run, repaired two environment and model-runtime failures, and resubmitted fine-tuning.

- Diagnosed the first staging failure as binary FASTA input passed to Biopython and fixed it with a text wrapper.
- Diagnosed the second staging failure as a regular-expression escape error and corrected CARD ARO identifier parsing.
- Confirmed the prepared dataset contains 4,401/943/944 CARD positives with equally sized control sets for train/dev/test.
- Fixed NumPy and scikit-learn incompatibility by pinning NumPy 1.24.4 in the cluster staging script.
- Disabled only DNABERT-2's optional Triton attention path because the cluster Python lacks Python.h; job 32221728 is pending on an A5000 resource.
=====
[2026-11-05 16:43:11]
Set up a CARD-based DNABERT-2 fine-tuning workflow and submitted it to the Bouchet cluster.

- Created local and cluster virtual environments, with package installation moved to an allocated CPU staging job to stay within login-node limits.
- Downloaded the current CARD data archive and verified it includes the nucleotide reference FASTA files required for training.
- Built a reproducible binary dataset pipeline using CARD resistance genes and coding-sequence controls from four bacterial reference genomes.
- Split positives by stable CARD ARO identifier so alleles from the same term cannot appear in more than one data split.
- Compared Slurm start estimates for H100, H200, L40S, and RTX 6000 GPUs, then submitted staging job 32220540 and dependent L40S fine-tuning job 32220541.
=====


