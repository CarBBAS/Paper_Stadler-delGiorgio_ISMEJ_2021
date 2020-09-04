rm(asv.tab, met.df, tax.tab)

comm <- otu_mat(pb)
#comm <- comm[,1:1000]

group <- sample_df(pb)$DnaType
permutations = 100
parallel = 2

perm<- permutations
N <- nobs

rm(pb)

getPermuteMatrix <- function(perm, N,  strata = NULL) {
    ## 'perm' is either a single number, a how() structure or a
    ## permutation matrix
    if (length(perm) == 1) {
      perm <- how(nperm = perm)
    }
    ## apply 'strata', but only if possible: ignore silently other cases
    if (!missing(strata) && !is.null(strata)) {
      if (inherits(perm, "how") && is.null(getBlocks(perm)))
        setBlocks(perm) <- strata
    }
    ## now 'perm' is either a how() or a matrix
    if (inherits(perm, "how"))
      perm <- shuffleSet(N, control = perm)
    else { # matrix: check that it *strictly* integer
      if(!is.integer(perm) && !all(perm == round(perm)))
        stop("permutation matrix must be strictly integers: use round()")
    }
    ## now 'perm' is a matrix (or always was). If it is a plain
    ## matrix, set minimal attributes for printing. This is a dirty
    ## kluge: should be handled more cleanly.
    if (is.null(attr(perm, "control")))
      attr(perm, "control") <-
        structure(list(within=list(type="supplied matrix"),
                       nperm = nrow(perm)), class = "how")
    perm
}


big.simper <- function (comm, group, permutations = 0, trace = FALSE, parallel = getOption("mc.cores"), 
          ...) 
{
  EPS <- sqrt(.Machine$double.eps)
  if (any(rowSums(comm, na.rm = TRUE) == 0)) 
    warning("you have empty rows: results may be meaningless")
  
  pfun <- function(x, comm, comp, i) {
    groupp <- group[perm[x, ]]
    ga <- comm[groupp == comp[i, 1], , drop = FALSE]
    gb <- comm[groupp == comp[i, 2], , drop = FALSE]
    ls.a <- lapply(seq_len(nrow(ga)), function(x) ga[x, , drop = F])
    ls.b <- lapply(seq_len(nrow(gb)), function(x) gb[x, , drop = F])
    n.a <- nrow(ga)
    n.b <- nrow(gb)
    for(j in seq_len(n.b)){
      for(k in seq_len(n.a)){
        md <- abs(ls.a[[k]][, , drop = F]) - ls.b[[j]][, , drop = F]
        me <- ls.a[[k]][, , drop = F] + ls.b[[j]][, , drop = F]
        
        contr.ls[[(j - 1) * n.a + k]] <- md/sum(me)
      }
    }
    contr <- do.call(rbind, contr.ls)
    colMeans(contr)
  }
  
  comm <- as.matrix(comm)
  comp <- t(combn(unique(as.character(group)), 2))
  outlist <- NULL
  P <- ncol(comm)
  nobs <- nrow(comm)
  perm <- getPermuteMatrix(permutations, nobs)
  if (ncol(perm) != nobs) 
    stop(gettextf("'permutations' have %d columns, but data have %d rows", 
                  ncol(perm), nobs))
  nperm <- nrow(perm)
  if (nperm > 0) 
    perm.contr <- matrix(nrow = P, ncol = nperm)
  if (is.null(parallel)) 
    parallel <- 1
  hasClus <- inherits(parallel, "cluster")
  isParal <- hasClus || parallel > 1
  isMulticore <- .Platform$OS.type == "unix" && !hasClus
  if (isParal && !isMulticore && !hasClus) {
    parallel <- makeCluster(parallel)
  }
  for (i in seq_len(nrow(comp))) {
    group.a <- comm[group == comp[i, 1], , drop = FALSE]
    group.b <- comm[group == comp[i, 2], , drop = FALSE]
    ls.a <- lapply(seq_len(nrow(group.a)), function(x) group.a[x, , drop = F])
    ls.b <- lapply(seq_len(nrow(group.b)), function(x) group.b[x, , drop = F])
    n.a <- nrow(group.a)
    n.b <- nrow(group.b)
    
    contr.ls <- list()
    
    for(j in seq_len(n.b)){
      for(k in seq_len(n.a)){
          md <- abs(ls.a[[k]][, , drop = F]) - ls.b[[j]][, , drop = F]
          me <- ls.a[[k]][, , drop = F] + ls.b[[j]][, , drop = F]
          
          contr.ls[[(j - 1) * n.a + k]] <- md/sum(me)
        }
    }
    
    contr <- do.call(rbind, contr.ls)
    average <- colMeans(contr)
    
    if (nperm > 0) {
      if (trace) 
        cat("Permuting", paste(comp[i, 1], comp[i, 2], 
                               sep = "_"), "\n")
      contrp.ls <- list()
      #contrp <- matrix(ncol = P, nrow = n.a * n.b)
      if (isParal) {
        if (isMulticore) {
          perm.contr <- mclapply(as.list(seq_len(nperm)), function(d) pfun(d, 
                                                                  comm, comp, i), mc.cores = parallel)
          perm.contr <- do.call(cbind, perm.contr)
        }
        else {
          perm.contr <- parSapply(parallel, seq_len(nperm), 
                                  function(d) pfun(d, comm, comp, i, contrp))
        }
      }
      else {
        perm.contr <- sapply(1:nperm, function(d) pfun(d, 
                                                       comm, comp, i, contrp))
      }
      p <- (rowSums(apply(perm.contr, 2, function(x) x >= 
                            average - EPS)) + 1)/(nperm + 1)
    }
    else {
      p <- NULL
    }
    overall <- sum(average)
    sdi <- apply(contr, 2, sd)
    ratio <- average/sdi
    ava <- colMeans(group.a)
    avb <- colMeans(group.b)
    ord <- order(average, decreasing = TRUE)
    cusum <- cumsum(average[ord]/overall)
    out <- list(species = colnames(comm), average = average, 
                overall = overall, sd = sdi, ratio = ratio, ava = ava, 
                avb = avb, ord = ord, cusum = cusum, p = p)
    outlist[[paste(comp[i, 1], "_", comp[i, 2], sep = "")]] <- out
  }
  if (isParal && !isMulticore && !hasClus) 
    stopCluster(parallel)
  attr(outlist, "permutations") <- nperm
  attr(outlist, "control") <- attr(perm, "control")
  class(outlist) <- "simper"
  outlist
}
#<bytecode: 0x55dde3a55320>
#  <environment: namespace:vegan>