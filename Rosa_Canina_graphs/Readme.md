# Graphs for Rosa Canina

## CENH3 ChIP-seq of Rosa canina
Sequences downloaded from here: https://www.ncbi.nlm.nih.gov/sra/SRX27753549%5Baccn%5D

```bash

sudo apt install sra-toolkit
prefetch SRX27753549
fastq-dump --threads 8 --split-files SRR32424402.sra 
```
## Genome assembly

Genome was downloaded from here: https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_051942975.1/
