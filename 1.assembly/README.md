# genome assembly


## the data

* The sequence data are available from the BioProject PRJNA734792 (NCBI, SRA)

## the initial assembly

- the reads were cleaned by using Trim Galore(version 0.4.4) removing 
Nextera adaptors, clipping 15 bp in 5’-end and 1 bp in 3’-end and 
trimming low-quality ends (phred score < 30)

- the assembly was carried out by using SPAdes (version 3.9.1; 
options: careful mode, automatic k-mers)

- the genomic contigs smaller than 1 kb were not considered

## identification of “apicomplexa” vs. contamination contigs

The contigs were analyzed by using a principal component analysis (PCA) 
based on their 5-mer composition, which allowed classifying them into 6 
groups by using a hierarchical clustering method (HCA) based on the 
Ward criterion.

- k-mers were counted with the Python script: 
`get_kmers_occurences_memory_optimization_freq.py`  (-k 5)

- k-mers including non-ATGC characters were removed with the following 
Shell script:

```
awk ' NR == 1 {for (ii=2;ii<=NF;ii++) 
                  if ($ii !~ "N" ) 
                      {xx++; zz[xx]=ii} 
              } 
              {printf "%s", $1; 
               for (kk=1;kk<=xx; kk++) 
                   {printf " %s", $zz[kk]}; 
                print ""
              }' genomes.sup1kb.kmers5 > genomes.sup1kb.kmers5_noN

```

- PCA and HCA analysis were conducted with the following R code:


```
data=read.table("genomes.sup1kb.kmers5_noN", head=T, row=1)
pca=prcomp(data)

png("pca.png", width = 800, height=1200)
par(mfrow=c(3,1))
plot(pca$x[,1:2], cex=0.5, pch=16)
plot(pca$x[,c(1,3)], cex=0.5, pch=16)
plot(pca$x[,c(2,3)], cex=0.5, pch=16)
dev.off()

png("pca_axes.png")
plot(pca)
dev.off()

d=dist(pca$x)
h=hclust(d, met="ward.D2")

png("h.png")
plot(h, labels=FALSE)
dev.off()

hc=cutree(h, 6)

png("k6.png", width = 800, height = 1200)
par(mfrow=c(3,1))
plot(pca$x[,1:2], cex=0.5, pch=16, col=hc)
plot(pca$x[,c(1,3)], cex=0.5, pch=16, col=hc)
plot(pca$x[,c(2,3)], cex=0.5, pch=16, col=hc)
dev.off()

write.table(hc, file="k6.txt", quote=FALSE)
```

- the putative protein coding genes were predicted by using Augustus 
(version 3.3; gene model: T. gondii)

- the predicted proteins were compared with the NCBI non-redondant 
protein database by using BLAST 


## identification of genomes A and B among “apicomplexa” contigs 

The contigs of the “apicomplexa” cluster were splitted into genomes A 
and B by using the difference of coverage observed for each of the 4 
gDNA libraries.

- the reads of each individual libraries were mapped against the “
apicomplexa” contigs  by using `Bowtie2`.

- the bam files were sorted with `Samtools`

- the bedgraph files were generated with the `Bedtools` (tool: `genomeCoverageBed`)

- the median coverages were computed with the following R script 
(have to be applied for each individual libraries):

```
data = read.table("library_1.bedgraph")
scaffolds = unique(data[,1])
medians = c()
for (i in 1:length(scaffolds))
	{
    s=as.character(scaffolds[i])
    print(paste(s, "-", i, sep=" " ))
    d = data[data[,1] == s,]
    cov = d[,3]
    med = c(s, median(cov))
    medians = rbind(medians, med)
	}
colnames(medians) = c("scaffolds", "median_coverage")
write.table(medians,"library_1.cov", row.names=F, sep="\t", quote=F)
```


