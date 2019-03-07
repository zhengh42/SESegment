

SESegment<-function(mat,coordinates,corr.spearman.thres,dist.euclid.thres,type="average",maxblocksize=10000){
  
  cluster<-1:ncol(mat)
  clusterdist<-rep(0,ncol(mat))
  clusterB.start=1
  clusterB.end=1
  
  ### clustering
  for(i in 2:ncol(mat)){
    if(clusterdist[i]==Inf) {next;}
    
    corr.spearman = CalcDistClusters(mat[,clusterB.start:clusterB.end],mat[,i],type = type, dist.type = "spearman")
    dist.euclid = CalcDistClusters(mat[,clusterB.start:clusterB.end],mat[,i],type = type, dist.type = "euclid")
    
    if(corr.spearman>=corr.spearman.thres && dist.euclid<=dist.euclid.thres && (coordinates[i,2] - coordinates[clusterB.start,1]) <= maxblocksize){
      cluster[i]=cluster[clusterB.start]
      clusterB.end = i
      clusterdist[(clusterB.start+1):clusterB.end] = Inf
    } else{
      clusterB.start=i
      clusterB.end=i
    }
  }
  
  ### matrix after clustering
  mataftercluster<-sapply(unique(cluster),function(x){rowMeans(as.data.frame(mat[,which(cluster==x)]))})
  
  ### coordinate after clustering
  coordinates$cluster=cluster
  coordinatesaftercluster<-data.frame(
    start=sapply(unique(coordinates$cluster),function(x){ min(coordinates[coordinates$cluster==x,][,"start"] ) }),
    end=sapply(unique(coordinates$cluster),function(x){ max(coordinates[coordinates$cluster==x,][,"end"] ) })
  )
  
  return(list("cluster"=cluster,"matrix"=mataftercluster,"coordinates"=coordinatesaftercluster))
}

PairwiseDistanceDistribution<-function(mat=mat){
  return(sapply(1:(ncol(mat)-1), function(x){sqrt(sum((mat[,x]-mat[,x+1])^2))}))
}

PairwiseCorrelationDistribution<-function(mat=mat){
  return(unlist(sapply(1:(ncol(mat)-1), function(x){cor(mat[,x],mat[,x+1],method = "spearman")})))
}

CalcDistClustPoint <- function(clust, pnt, type = "average", dist.type = "spearman"){
  stopifnot(is.element(type, c("single", "complete", "average")), is.element(tolower(dist.type), c("spearman", "euclid")), dim(clust)[1] == length(pnt))
  clust <- as.matrix(clust)
  size <- ncol(clust)
  if (size == 1){
    if (is.element(dist.type , c("spearman"))) return( 1 - abs(cor(clust, pnt, use = "complete.obs", method = dist.type))) else
      if (dist.type == "euclid") {
        inds.rm <- union(which(is.na(clust)), which(is.na(pnt)))
        if (length(inds.rm) > 0) return(sqrt(sum((clust[-inds.rm] - pnt[-inds.rm])^2))) 
        else return(sqrt(sum((clust - pnt)^2)))
      }
  }else{
    distances <- apply(clust, 2, function(x){
      if (is.element(dist.type , c("spearman"))) return( 1 - abs(cor(x, pnt, use = "complete.obs", method = dist.type))) else
        if (dist.type == "euclid") {
          inds.rm <- union(which(is.na(x)), which(is.na(pnt)))
          if (length(inds.rm) > 0) return(sqrt(sum((x[-inds.rm] - pnt[-inds.rm])^2))) 
          else return(sqrt(sum((x - pnt)^2)))
        }
    })
    
    if (type == "single") return(min(distances)) else
      if (type == "complete") return(max(distances)) else
        return(mean(distances)) 
  }
}

CalcDistClusters <- function(clust.1, clust.2, type = "single", dist.type = "spearman"){
  stopifnot(is.element(type, c("single", "complete", "average")), is.element(tolower(dist.type), c("spearman","euclid")))
  clust.1 <- as.matrix(clust.1)
  size.1 <- ncol(clust.1)
  if (size.1 == 1) return(CalcDistClustPoint(clust.2, clust.1, type = type, dist.type = dist.type))
  distances <- apply(clust.1, 2, function(x){
    CalcDistClustPoint(clust.2, x, type = type, dist = dist.type)
  })
  if (type == "single") return(min(distances)) else
    if (type == "complete") return(max(distances)) else
      return(mean(distances))
}

