# the assembly

1. the reads were cleaned by using Trim Galore(version 0.4.4) removing 
Nextera adaptors, clipping 15 bp in 5’-end and 1 bp in 3’-end and 
trimming low-quality ends (phred score < 30)

2. the assembly was carried out by using SPAdes (version 3.9.1; 
options: careful mode, automatic k-mers)

3. the genomic contigs smaller than 1 kb were not considered

