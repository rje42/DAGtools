##' Get topological order
##'
##' @param graph adjacency matrix for directed graph
##'
topOrder <- function(graph) {

  p <- nrow(graph)
  left <- rep(TRUE, p)

  out <- c()
  cs <- colSums(graph)

  while (any(left)) {
    wh <- which(left & cs == 0)
    if (length(wh) == 0) stop("Cyclic graph")
    out <- c(out, wh)
    left[wh] = FALSE
    cs <- colSums(graph*left)
  }

  out
}

##' Fix Gaussian distribution by dividing by p(a | b)
##'
##' @param Sigma covariance matrix
##' @param a,b sets
##' @param cov_only logical (currently \code{TRUE} is only acceptable value)
##'
##' Returns covariance matrix of same dimensions as Sigma, with
##' dependence of variables in \code{a} on variables in \code{b}
##' removed.
##'
##' @export
intervene <- function(Sigma, a, b, cov_only = TRUE) {
  if (!cov_only) stop("Conditional means not implemented yet")

  p <- nrow(Sigma)
  C <- seq_len(p)[-c(a,b)]  # remaining variables

  out_inv <- solve.default(Sigma)
  ## set relevant correlations to zero
  if (length(b) > 0) {
    Sigma[a,a] = schur(Sigma, a, z=b)
    Sigma[a,b] = Sigma[b,a] = 0
  }
  if (length(C) == 0) return(Sigma)

  ## construct matrix with correct marginal and conditional parts
  out_inv[c(a,b),c(a,b)] <- solve.default(Sigma[c(a,b),c(a,b)]) +
    out_inv[c(a,b), C, drop=FALSE] %*% solve.default(out_inv[C, C]) %*% out_inv[C, c(a,b), drop=FALSE]

  solve.default(out_inv)
}

##' Generate sample from posterior distribution of (Gaussian) BN
##'
##' @param A adjacency matrix for graph
##' @param xbar sample mean of data
##' @param S sample variance of data (not normalized by \code{n})
##' @param n number of samples
##' @param alpha prior parameter, default to dimension of data
##' @param T0 prior precision, defaults to identity matrix
##'
##'
genSamp <- function(A, xbar, S, n, alpha, T0) {
  # require(igraph)

  p <- length(xbar)
  if (missing(alpha)) alpha = p
  if (missing(T0)) T0 = diag(p)

  nu <- 1   ## tuning parameter (same as 'am' in BiDAG package)

  Tn <- T0 + S + n*nu/(n + nu)*(xbar %*% t(xbar))

  ## organise entries topologically
  # topOrd <- topological.sort(graph_from_adjacency_matrix(A, "directed"))
  topOrd <- topOrder(A)
  invOrd <- order(topOrd)
  A <- A[topOrd,topOrd]
  if (any(A[lower.tri(A)] == 1)) stop("Topological sort failed")
  Tn <- Tn[topOrd,topOrd]

  ## generate distribution from saturate model
  W <- rWishart(1, alpha+n, solve.default(Tn))[,,1]   # inverse covariance
  Sig <- solve.default(W)   # covariance
  Out <- matrix(0, p, p)

  for (i in seq_len(p)) {
    pa_i <- which(A[,i] > 0)
    all_i <- c(pa_i, i)
    to_i <- seq_len(i)

    # sub <- solve.default(Sig[to_i,to_i])
    to_add <- solve.default(Sig[all_i,all_i])
    to_add <- to_add[nrow(to_add),]

    Out[i,i] <- 1
    curr <- solve.default(Out[to_i,to_i,drop=FALSE])
    curr[i,all_i] <- curr[all_i,i] <- to_add
    curr[-i,-i] <- curr[-i,-i,drop=FALSE] + curr[-i,i,drop=FALSE] %*% curr[i,-i,drop=FALSE] / curr[i,i]

    Out[to_i, to_i] <- solve.default(curr)
    #
    # if (length(pa_i) > 0) {
    #   sub[-1,-1] <- solve.default(curr) + sub[-1,1,drop=FALSE] %*% sub[1,-1,drop=FALSE]/sub[1,1]
    # }
    # Out[all_i,all_i] <- solve(sub)
  }

  ## return to original vertex ordering
  Out <- Out[invOrd,invOrd]

  return(Out)
}



##' Get matrix of marginal causal effects for every pair of variables
##'
##' @param A adjacency matrix of DAG
##' @param Sigma covariance matrix (assumed Markov with respect to \code{A})
##'
##' @details Gaussian data only.
##'
##' @return gives matrix whose i,j th entry is total
##' effect of i on j.
##'
##' @export
causalEffects <- function(A, Sigma) {

  p <- nrow(A)
  out <- out2 <- diag(p)
  if (p == 0L) return(list(total=out, direct=out))

  # topOrd <- topological.sort(graph_from_adjacency_matrix(A, "directed"))
  ## see if graph is non-topological; if so, then re-sort vertices
  if (any(A[lower.tri(A)] > 0)) {
    topOrd <- topOrder(A)
    invOrd <- order(topOrd)

    A <- A[topOrd,topOrd]
    if (any(A[lower.tri(A)] != 0)) {
      stop("topological error")
    }

    Sigma <- Sigma[topOrd,topOrd]
  }
  else {
    topOrd <- invOrd <- seq_len(p)
  }

  ## get SEM parameters from Sigma
  tmp <- params(Sigma)

  ## go through each variable, intervene given its parents,
  ## and record total effects on other variables.
  for (i in seq_len(p)) {
    # pa_i <- which(A[,i] > 0)
    ch_i <- which(A[i,] > 0)
    # newSig <- intervene(Sigma, i, pa_i)

    tmpB <- tmp$B
    tmpB[,i] <- 0  # set incoming edges to zero
    IBi <- solve.default(diag(p) - tmpB)
    newSig2 <- t(IBi) %*% tmp$Om %*% IBi  # reconstruct new Sigma

    # print(i)
    # print(round(newSig-newSig2,6))

    out[i,-i] = newSig2[i,-i]/newSig2[i,i]

    for (j in ch_i) {
      pa_j <- which(A[,j] > 0)
      tmp2 <- solve(newSig2[c(pa_j,j),c(pa_j,j)])
      tmp2 <- tmp2/sqrt(outer(diag(tmp2), diag(tmp2)))
      idx <- which(pa_j == i)

      out2[i,j] <- -tmp2[idx,length(pa_j)+1]
    }
  }

  out <- out[invOrd,invOrd]
  out2 <- out2[invOrd,invOrd]

  return(list(total=out, direct_normal=out2))
}

# plotCausalEffects <- function(x, ...) {
#
#   tol <- 1e-8
#   pz <- sum(abs(x) < tol)/length(x)
#   x <- x[abs(x) >= tol]
#
#   out <- hist(x, breaks=seq(from=-1,to=1,by=0.02), plot = FALSE)
#   out$density <- out$density*(1-pz)/sum(out$density)
#
#   ylim <- c(0,0.05) # range(out$density)
#   # ylim[2] <- max(ylim[2], pz)
#   xlim <- range(c(-0.6,0.8,out$breaks))
#
#   plot(out, col=2, freq=FALSE, xlim=xlim, ylim=ylim, main="", axes=FALSE, ...)
#   axis(side=1, at=round(seq(-0.6,0.8,by=0.2),1))
#   axis(side=2)
#   # lines(c(0,0), c(0,pz), lwd=2)
#   # points(0,pz,pch=20)
#   # text(-0.1, pz - ylim[2]/10, labels=as.character(round(pz,2)))
# }

##' Augment MCMC
##'
##' Adds samples and causal effects to partition MCMC output
##' for Gaussian data.
##'
##' @param object (Gaussian) output of BiDAG MCMC
##' @param dat data frame used in sampling
##' @param ... other arguments
##' @param verbose logical indicating whether to show progress
##' @param force logical - can this be used on previously augmented chains?
##'
##' @export
augment_mcmc <- function(object, dat, ..., verbose=TRUE, force=FALSE) {
  if (class(object) == "aug_chain") {
    if (force) {
      object <- object[-c("Sigma", "CauEff", "betas", "Om")]
    }
    else {
      message("Chain already augmented, exiting")
      return(object)
    }
  }
  B <- object$B
  n <- nrow(dat)
  p <- object$p

  if (is.null(object$varnames)) object$varnames <- names(dat)

  ## use raw moments since mean already centred
  samp_means <- colMeans(dat, na.rm = TRUE)
  S <- (n-1)*cov(dat, use = "pairwise") + n*outer(samp_means, samp_means)
  Sigma <- CauEff <- CauEff2 <- betas <- Om <- list()

  ## go through sampled graphs, sampling from posterior
  ## distribution of each model.  Record causal effects
  for (i in seq_len(B)) {
    gr <- as.matrix(object$incidence[[i]])
    Sigma[[i]] <- genSamp(gr, rep(0,p), S, n)
    tmp <- params(Sigma[[i]], gr)
    betas[[i]] <- tmp$B
    Om[[i]] <- tmp$Om
    tmp2 <- causalEffects(gr, Sigma[[i]])
    CauEff[[i]] <- tmp2$total
    CauEff2[[i]] <- tmp2$direct_normal

    if (verbose) rje::printPercentage(i, B)
  }

  betas <- array(unlist(betas), c(p,p,B))
  Om <- array(unlist(Om), c(p,p,B))
  CauEff <- array(unlist(CauEff), c(p,p,B))
  CauEff2 <- array(unlist(CauEff2), c(p,p,B))
  dimnames(CauEff) <- dimnames(CauEff2) <- dimnames(betas) <- dimnames(Om) <- list(object$varnames, object$varnames, NULL)

  object <- c(object, list(Sigma=Sigma, CauEff=CauEff, CauEff2=CauEff2, betas=betas, Om=Om))
  class(object) <- c("aug_chain")
  object
}


##' Merge two MCMC chains into one
##'
##'
##' @param object chain or list of chains
##' @param ... other chains
##'
##' Combine separate chains into a single sample.  Should
##' work for \code{MCMCchain}, \code{aug_chain} or \code{summary_MCMCchain}.
##'
##' @export
merge_chain <- function(object, ...) {

  ## for a single chain, put in a list
  if (class(object) %in% c("aug_chain", "summary_MCMCchain", "MCMCchain")) {
    object = c(list(object), list(...))
  }
  else object = c(object, list(...))

  out <- purrr::transpose(object)

  ## sort out info
  out$info <- purrr::transpose(out$info)
  out$info$iterations <- sum(unlist(out$info$iterations))
  out$info$samplesteps <- sum(unlist(out$info$samplesteps))
  out$info$DBN <- out$info$DBN[[1]]
  out$info$algo <- out$info$algo[[1]]
  out$info$spacealgo <- out$info$spacealgo[[1]]
  out$info$sampletype <- out$info$sampletype[[1]]

  ## flatten lists of doubles
  wh <- !(names(out) %in% c("adj", "info"))
  out[wh] <- Map(purrr::flatten, out[wh])
  out$partitionscores <- unlist(out$partitionscores)

  ## unique ones should just be given once
  uq <- names(out) %in% c("varnames", "varnames_short", "p", "n")
  out[uq] <- object[[1]][uq]

  ## deal with adjacency and sparsity
  if (!is.null(out$adj)) {
    w_adj <- mapply(`*`, out$adj, out$B, SIMPLIFY = FALSE)
    out$adj <- reduce(w_adj, `+`)/sum(unlist(out$B))
  }
  if (!is.null(out$sparsity)) out$sparsity <- sum(mapply(`*`, out$sparsity, out$B))/sum(unlist(out$B))
  if (!is.null(out$B)) out$B <- sum(unlist(out$B))

  ## set dimensions correctly for arrays
  wh_arr <- names(out) %in% c("Sigma", "CauEff", "CauEff2", "betas", "Om")
  if (any(wh_arr)) {
    out[wh_arr] <- Map(unlist, out[wh_arr])
    out[wh_arr] <- Map(function(x) `dim<-`(x, value=c(out$p,out$p,out$B)), out[wh_arr])
  }

  # array(unlist(adj), c(rep(p,2)))

  ## flatten things into single lists
  # wh <- purrr::flatten_chr(Map(class, object[[1]])) != purrr::flatten_chr(Map(class, out))
  # out[wh] <- Map(flatten, out[wh])

  # out$n <- out$n[1]
  # out$p <- out$p[1]

  class(out) = class(object[[1]])

  out
}
