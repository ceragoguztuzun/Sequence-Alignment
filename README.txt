Use "make" command for compilation, followed by the command line examples below. 

Command line examples: Each line represents different runs for different modes:
• ./allalign --mode global --input sequences.fasta --gapopen -5
• ./allalign --mode aglobal --input sequences.fasta --gapopen -5 --gapext -2
• ./allalign --mode local --input sequences.fasta --gapopen -5
• ./allalign --mode alocal --input sequences.fasta --gapopen -5 --gapext -2

• Same input file ”sequences.fasta”, which contains two DNA sequences.
Parameters:
• −−mode: It will be selected from one of the followings:
• global: Needleman-Wunsch with naı̈ve gap scoring
• local: Smith-Waterman with naı̈ve gap scoring
• aglobal: Needleman-Wunsch with affine gap scoring
• alocal: Smith-Waterman with affine gap scoring
• −−input: Input FASTA file for sequences
• −−gapopen: Gap opening penalty for affine gap model, or unit gap cost for naı̈ve
model
• −−gapext: Gap extension penalty for affine gap model

---------------
Çerağ Oğuztüzün
21704147
20.11.2020
---------------