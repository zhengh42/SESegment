---
title: "Segmenting CpGs based on beta values using Spearman correlation and euclidean distance"
output: html_document
---

### Load scripts and example dataset

The example data set contains 188 samples (in row) and 1000 CpG sites (in column) from one chromosome. This script is run per-chromosome.

```{r setup}
source("SESegment.R")
load("datasets/coordinate.Robj") # coordinate for single chromosome
load("datasets/mat.Robj") # beta matrix for single chromosome
load("datasets/matall.Robj") # coordinate and  beta matrix for all chromosome
```

### Set thresholds

Spearman correlation and euclidean distance between two CpG sites/segments are calculated across samples, in this example data set, 188 samples.

If two neighbouring CpG sites/segments have higher correlation than the correlation threshold and lower Euclidean distance than the distance threshold, the two CpG sites/blocks are clustered together in one segment. The methylation level (beta values) of a CpG segment is averaged across the CpG sites within a segment.

The function `PairwiseDistanceDistribution` and `PairwiseCorrelationDistribution` can be used to calculate pair-wise distance and correlation for neighbouring CpG sites/blocks.

```{r set_threshold}
corr.spearman.thres=0.3
dist.euclid.thres=2
summary(PairwiseDistanceDistribution(mat))
summary(PairwiseCorrelationDistribution(mat))
```

### Run clustering (per chromosome)

- mat: the input beta value matrix, samples in row and CpGs in column.
- coordinates, the coordinates on the chromosome for each CpG site/segment. The start and end postion are both specified.
- type: how the distance/correlaton, if a CpG block with multiple CpG sites are involved. Possible value ["single","maximum","average"]
- maxblocksize, the maxumun size of a CpG blocks, in terms of corordiates on the genome.

```{r run_cluster}
out<-list()
out[[1]]<-list("matrix"=mat,"coordinates"=coordinate)

for(i in 1:10){
  out[[i+1]]<-SESegment(mat=out[[i]]$matrix,coordinates=out[[i]]$coordinates,corr.spearman.thres=corr.spearman.thres,dist.euclid.thres=dist.euclid.thres,type="average",maxblocksize=10000)
}

sapply(1:11,function(i){mean(PairwiseDistanceDistribution(out[[i]]$matrix)<=dist.euclid.thres)})
sapply(1:11,function(i){mean(PairwiseCorrelationDistribution(out[[i]]$matrix)>=corr.spearman.thres)})
sapply(1:11,function(i){dim(out[[i]]$matrix)})
```

The algorithm is run iteratively 10 times. The output of the previous run is supplied to the next run. The final clustering result is chosen if no further CpG sites/blocks can be clustered, in this case, the 6th iteration. Thus, the final result:

The matrix: `out[[6]]$matrix`  
The coordinates: `out[[6]]$coordinates`

### Run clustering (all chromosomes)

```{r run_cluster_all,echo=T,eval=F}
### Initial objects for output
matall.out<-NULL
coordinate.out<-NULL

### Run SESegment for each chromosome sequentially, and aggregate the result
for (i in paste0("chr",c(1:22,"X","Y"))){
  mat<-matall %>% filter(CpG_chrm==i) %>% droplevels()
  coordinate<-mat[,c("CpG_beg","CpG_end")]
  colnames(coordinate)<-c("start","end")
  mat[,1:4]<-NULL
  mat<-t(mat)

  colsExclude<-which(apply(mat,2,function(x){sd(x,na.rm = T)})==0)
  if(length(colsExclude)>0){
    mat=mat[,-colsExclude]
    coordinate=coordinate[-colsExclude,]
  }
  
  out<-list()
  out[[1]]<-list("matrix"=mat,"coordinates"=coordinate)
  for(j in 1:8){
    out[[j+1]]<-SESegment(mat=out[[j]]$matrix,coordinates=out[[j]]$coordinates,corr.spearman.thres=corr.spearman.thres,dist.euclid.thres=dist.euclid.thres,type="average",maxblocksize=10000)
  }
  
  coordinate.chr<-data.frame(chr=i, start=out[[9]]$coordinates[,1],end=out[[9]]$coordinates[,2])
  matall.out<-cbind(matall.out,out[[9]]$matrix)
  coordinate.out<-rbind(coordinate.out,out[[9]]$coordinate.chr)
}

```

The final beta value matrix for CpG sites/segments and their coordinates can be found in matall.out and coordinate.out.
