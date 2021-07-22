# identification of genomes A and B among “apicomplexa” contigs 

The contigs of the “apicomplexa” cluster were splitted into genomes A 
and B by using the difference of coverage observed for each of the 4 
gDNA libraries.




## median coverage for each individual libraries

1. the reads of each individual libraries were mapped against the “
apicomplexa” contigs  by using `Bowtie2`.

2. the bam files were sorted with `Samtools`

3. the bedgraph files were generated with the `Bedtools` (tool: `genomeCoverageBed`)

4. the median coverages were computed with the following R script 
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

5. the median coverages were merged in a tabulated file 

