# identification of “apicomplexa” vs. contamination contigs

The contigs were analyzed by using a principal component analysis (PCA) 
based on their 5-mer composition, which allowed classifying them into 6 
groups by using a hierarchical clustering method (HCA) based on the 
Ward criterion.


1. k-mers were counted with the Python script: 
`get_kmers_occurences_memory_optimization_freq.py`  (-k 5)

2. k-mers including non-ATGC characters were removed with the following 
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

3. PCA and HCA analysis were conducted with the following R code:


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

4. the putative protein coding genes were predicted by using Augustus 
(version 3.3; gene model: T. gondii)

5. the predicted proteins were compared with the NCBI non-redondant 
protein database by using BLAST 



