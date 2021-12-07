# identification of genomes A and B among “apicomplexa” contigs 

The contigs of the “apicomplexa” cluster were splitted into two sets, genomes A 
and B, by using the difference of coverage observed for each of the 4 
genomic DNA (gDNA) libraries (JS-470, JS-482, JS-488 and JS-489).


## median coverage for each individual libraries

1. the reads of each individual libraries were mapped against the “
apicomplexa” contigs  by using `bowtie2` 
(this step is repeated for the 4 gDNA libraries).

2. the bam files were sorted with `samtools` (tool: `sort`,this step is 
repeated for the 4 gDNA libraries)

3. the bedgraph files were generated with the `bedtools` 
(tool: `genomeCoverageBed`, this step is repeated for the 4 gDNA libraries)

4. the median coverages were calculated for each contig by using the 
following R script (this step is repeated for all libraries):

```
data=read.table("library_1.bedgraph")
lib1_median=as.data.frame(tapply(data$V3, data$V1,median) )
colnames(lib1_median)="coverage_lib1"
write.table(lib1_median, "lib1_median.csv", quote=F)
```

5. the median coverages obtained for each gDNA libraries were merged into a tabulated file 


## identification of contigs from genomes A and B

1. A PCA was computed from the 4 median coverage values observed for each contigs (in R):

```
library(MASS)

data=read.table("median_lib1_lib2_lib3_lib4.csv",row=1)
#  principal component analysis
pca=prcomp(data,scale=T)
```

2. An initial binary classification was based on k-means algorithm:

```
#  initial classification
mykmeans=kmeans(pca$x, centers=2)
```

3.  A linear discriminant method (training and classification) was 
iteratively repeated 3 times until convergence (manually checked):

```
#  first model
myldak=lda(pca$x, mykmeans$cluster)
mypredk=predict(myldak, pca$x)
table(mypredk$class, mykmeans$cluster)

# second model
myldak2=lda(pca$x, mypredk$class)
mypredk2=predict(myldak2, pca$x)
table(mypredk$class, mypredk2$class)

# third model
myldak3=lda(pca$x, mypredk2$class)
mypredk3=predict(myldak3, pca$x)
table(mypredk2$class, mypredk3$class)

# fourth model (to check convergence)
myldak4=lda(pca$x, mypredk3$class)
mypredk4=predict(myldak4, pca$x)
table(mypredk3$class, mypredk4$class)

```

## identification of contigs corresponding to genome A/B hybrids

1. The median coverage have been computed for all 1kb non-overlapping 
windows along the contigs

2. The calculated values for the windows were analyzed using the same 
protocol as described above for the contigs
     1. A PCA computed from the 4 median coverage values observed for 
     each window
     2. An initial binary classification based on k-means algorithm
     3. A linear discriminant method (training and classification) 
iteratively repeated until convergence (manually checked)

3. The contigs with regions identified as A and B were splitted




