# Simulence
A Python2 script for generating simulated RNA-seq data with realistic read coverage profiles

## Usage
```
python Simulence.py FASTQ input.fa 100 10 output.fq
```
Sys arguments: 
* FASTQ       Output format. Options: FASTA/FASTQ 
* input.fa    Input file. File must be in FASTA format
* 100         Read length. No size limitations, but read length must be shorter than sequences provided.
* 10          Fold coverage
* output.fq   Output file name

## Info
* Simulated reads are mutated (substitutions/INDELs) to replicate Illumina sequencing.
* Substitution rate per nucleotide: 0.0035
* INDEL rate per nucleotide: 3.95e-6
* Substitution and INDEL rates taken from: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4787001/

## Contributors

* Paul Donovan, PhD (email: pauldonovandonegal@gmail.com)

## License

This project is licensed under the MIT License.
