
# SET partitioning framework proposed by:
# Schmera, D., Podani, J., & Legendre, P. (2020). 
# What do beta diversity components reveal from presence-absence community data? 
# Let us connect every indicator to an indicandum! 
# Ecological Indicators, 117(May), 106540. https://doi.org/10.1016/j.ecolind.2020.106540

# Script provided in Supplementary

# Partitioning of Weiher-Boylen beta diversity using SET partitioning
# Input:
# presence-absence matrix, where rows are sites, columns are species
# Output:
# first distance matrix: Weiher-Boylen beta diversity
# second distance matrix: Intersection of nestedness and beta diversity
# third distance matrix: Relative complement of nestedness in beta diversity
setpart.wb<-function(mat)
{
  mat <- as.matrix(mat)
  n <- nrow(mat)
  mat.b <- ifelse(mat>0, 1, 0)
  a <- mat.b %*% t(mat.b)
  b <- mat.b %*% (1 - t(mat.b))
  c <- (1 - mat.b) %*% t(mat.b)
  min.bc <- pmin(b,c)
  BD.wb<-(b+c)			
  I.wb<-matrix(0,n,n)
  RC.wb<-matrix(0,n,n)
  for (i in 2:n){
    for (j in 1:(i-1)){
      aa=a[i,j];bb=b[i,j];cc=c[i,j]
      if (aa==0) I.wb[i,j]<-0 else I.wb[i,j]<-abs(bb-cc)
      if (aa==0) RC.wb[i,j]<-(bb+cc) else RC.wb[i,j]<-2*min(bb,cc)
    }}
  # adding a few lines, as row and columnnames are getting lost in loop
  row.names(I.wb) <- row.names(BD.wb); row.names(RC.wb) <- row.names(BD.wb)
  colnames(I.wb) <- colnames(BD.wb); colnames(RC.wb) <- colnames(BD.wb)
  
  res<-list(as.dist(BD.wb),as.dist(I.wb),as.dist(RC.wb))
  res
}


# Partitioning of Jaccard dissimilarity using SET partitioning
# Input:
# presence-absence matrix, where rows are sites, columns are species
# Output:
# first distance matrix: Jaccard dissimilarity (beta diversity)
# second distance matrix: Intersection of nestedness and beta diversity
# third distance matrix: Relative complement of nestedness in beta diversity
setpart.j<-function(mat)
{
  mat <- as.matrix(mat)
  n <- nrow(mat)
  mat.b <- ifelse(mat>0, 1, 0)
  a <- mat.b %*% t(mat.b)
  b <- mat.b %*% (1 - t(mat.b))
  c <- (1 - mat.b) %*% t(mat.b)
  min.bc <- pmin(b,c)
  BD.j<-(b+c)/(a+b+c)			
  I.j<-matrix(0,n,n)
  RC.j<-matrix(0,n,n)
  for (i in 2:n){
    for (j in 1:(i-1)){
      aa=a[i,j];bb=b[i,j];cc=c[i,j]
      if (aa==0) I.j[i,j]<-0 else I.j[i,j]<-abs(bb-cc)/(aa+bb+cc)
      if (aa==0) RC.j[i,j]<-(bb+cc)/(aa+bb+cc) else RC.j[i,j]<-2*min(bb,cc)/(aa+bb+cc)
    }}
  
  # adding a few lines, as row and columnnames are getting lost in loop
  row.names(I.j) <- row.names(BD.j); row.names(RC.j) <- row.names(BD.j)
  colnames(I.j) <- colnames(BD.j); colnames(RC.j) <- colnames(BD.j)
  
  res<-list(as.dist(BD.j),as.dist(I.j),as.dist(RC.j))
  res
}


# Partitioning of Sørensen dissimilarity using SET partitioning
# Input:
# presence-absence matrix, where rows are sites, columns are species
# Output:
# first distance matrix: Sørensen dissimilarity (beta diversity)
# second distance matrix: Intersection of nestedness and beta diversity
# third distance matrix: Relative complement of nestedness in beta diversity
setpart.s<-function(mat)
{
  mat <- as.matrix(mat)
  n <- nrow(mat)
  mat.b <- ifelse(mat>0, 1, 0)
  a <- mat.b %*% t(mat.b)
  b <- mat.b %*% (1 - t(mat.b))
  c <- (1 - mat.b) %*% t(mat.b)
  min.bc <- pmin(b,c)
  BD.s<-(b+c)/(2*a+b+c)			
  I.s<-matrix(0,n,n)
  RC.s<-matrix(0,n,n)
  for (i in 2:n){
    for (j in 1:(i-1)){
      aa=a[i,j];bb=b[i,j];cc=c[i,j]
      if (aa==0) I.s[i,j]<-0 else I.s[i,j]<-abs(bb-cc)/(2*aa+bb+cc)
      if (aa==0) RC.s[i,j]<-(bb+cc)/(2*aa+bb+cc) else RC.s[i,j]<-2*min(bb,cc)/(2*aa+bb+cc)
    }}
  
  # adding a few lines, as row and columnnames are getting lost in loop
  row.names(I.s) <- row.names(BD.s); row.names(RC.s) <- row.names(BD.s)
  colnames(I.s) <- colnames(BD.s); colnames(RC.s) <- colnames(BD.s)
  
  res<-list(as.dist(BD.s),as.dist(I.s),as.dist(RC.s))
  res
}
