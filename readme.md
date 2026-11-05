# Contrastive Robustness for Genome Language Models

DNA is a string of nucleotides A, C, T and G. Nucleotide triplets are called \textit{codons}, most of which correspond to an amino acid, the fundamental element of proteins. There are 64 different codons; however, they only make 20 distinct types of amino acids. Because some codons make the same amino acid, it may be harmless to modify a nucleotide in a DNA sequence if the final protein turns out to be identical.

Genome Language Models (gLMs) such as DNABERT-2 are transformers pretrained on DNA sequences. They can be fine-tuned for downstream tasks, such as classifying what DNA sequences lead to antibiotic-resistant proteins, and which do not.

Our goal for this project is to identify whether existing gLMs are robust to input nucleotide modifications in the task of classifying for antibacterial resistance. DNABERT-2 and the Nucleotide Transformer will be our baselines. Because we know which codons have synonyms and which do not, we could cheaply create a variety of input sequences.