# identification of genomes A and B among “apicomplexa” contigs 

The contigs of the “apicomplexa” cluster were splitted into genomes A 
and B by using the difference of coverage observed for each of the 4 
gDNA libraries.


## median coverage for each individual libraries

1. the reads of each individual libraries were mapped against the “
apicomplexa” contigs  by using `Bowtie2` 
(have to be repeated for all libraries).

2. the bam files were sorted with `Samtools` (tool: `sort`,have to be 
repeated for all libraries))

3. the bedgraph files were generated with the `Bedtools` 
(tool: `genomeCoverageBed`, have to be repeated for all libraries)

4. the median coverages for each contigs were calculated with the 
following R script (have to be repeated for all libraries):

```
data=read.table("library_1.bedgraph")
lib1_median=as.data.frame(tapply(data$V3, data$V1,median) )
colnames(lib1_median)="coverage_lib1"
write.table(lib1_median, "lib1_median.csv", quote=F)
```

5. the median coverages were merged into a tabulated file 


## identification of genomes A and B

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

# fourth model
myldak4=lda(pca$x, mypredk3$class)
mypredk4=predict(myldak4, pca$x)
table(mypredk3$class, mypredk4$class)

```







