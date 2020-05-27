##' Identify SEM parameters for a DAG
##'
##' @param Sigma covariance matrix (either in topological order or with \code{graph} specified)
##' @param graph optionally, adjacency matrix object for underlying graph
##'
##' @details The order of variables is assumed toplogical if the graph
##' is not specified.
##'
##' @return list of \code{B} (upper triangular) and \code{Om} (diagonal)
##' such that if IB = (I-B)^{-1} then \code{t(IB) \%*\% Om \%*\% IB} gives \code{Sigma}.
##'
##' @export
params <- function(Sigma, graph) {
  p <- nrow(Sigma)

  ## if no graph specifed, assume order is already topological
  if (!missing(graph) && any(graph[lower.tri(graph)] == 1)) {
    ## otherwise, re-arrange and start again
    topOrd <- topOrder(graph)
    invOrd <- order(topOrd)

    out <- Recall(Sigma[topOrd,topOrd])
    out$B <- out$B[invOrd, invOrd]
    out$Om <- out$Om[invOrd, invOrd]
    return(out)
  }

  IB <- solve.default(chol.default(Sigma))
  Om_inv <- diag(IB)^2
  IB <- IB/rep(diag(IB), each=p)
  # IB[abs(IB) < 1e-8] = 0

  # print(solve(IB %*% diag(Om_inv) %*% t(IB)) - Sigma)

  # Li <- solve(L)
  ## print(t(Li) %*% Li - Sigma) # zero

  Om <- 1/Om_inv
  B <- diag(p) - IB

  # Om <- diag(U)^2
  # B <- diag(p) - U/diag(U)
  return(list(B=B, Om=diag(Om, nrow=length(Om))))
}

##' Gaussianize Data
##'
##' @param data data set (variables as columns)
##' @param jitter logical: if \code{TRUE}, ties are broken by jittering
##'
##' @export
gaussianize <- function(data, jitter=FALSE) {

  p <- ncol(data)
  n <- nrow(data)

  for (i in seq_len(ncol(data))) {
    if (jitter) {
      diffs <- min(diff(sort.default(unique.default(data[,i]))))
      data[,i] <- data[,i] + rnorm(n, sd=diffs/5)
    }
    data[,i] <- qnorm(rank(data[,i])/(n + 1))
  }

  data
}

# ##' Get edges from index of entry in adjacency matrix
# ##'
# ##' @param idx entry in array
# ##' @param nms character vector of variable names
# ##'
# ##' @examples
# ##' edge_names(4, LETTERS[1:3])
# ##' @export
# edge_names <- function(idx, nms) {
#   p <- length(nms)
#   if (any(idx > p^2)) stop("Indices too large")
#   ## assumes we only look at scales
#   edg <- cbind(row(diag(p))[idx], col(diag(p))[idx])
#   cbind(nms[edg[,1]], nms[edg[,2]])
# }

# ##' highlight differences between graphs
# ##'
# ##' @param g1,g2 adjacency matrices
# ##' @export diff_graph
# diff_graph <- function(g1, g2, show_common=TRUE) {
#   comm <- g1 > 0 & g2 > 0
#   cat(sum(comm)/2, "edges in common")
#   if (show_common) {
#     cat(":\n")
#     common_edges <- edge_names(which(comm & upper.tri(comm)))
#     cat(apply(common_edges, 1, function(x) paste(x, collapse = " -- ")), sep="\n")
#   }
#   else cat("\n")
#   cat("\n")
#
#   g1n2 <- g1 > 0 & g2 == 0
#   g1n2_edges <- edge_names(which(g1n2 & upper.tri(comm)))
#   cat(nrow(g1n2_edges), "edges in g1, but not g2:\n")
#   cat(apply(g1n2_edges, 1, function(x) paste(x, collapse = " -- ")), sep="\n")
#   cat("\n")
#
#   g2n1 <- g1 == 0 & g2 > 0
#   g2n1_edges <- edge_names(which(g2n1 & upper.tri(comm)))
#   cat(nrow(g2n1_edges), "edges in g2, but not g1:\n")
#   cat(apply(g2n1_edges, 1, function(x) paste(x, collapse = " -- ")), sep="\n")
#
#   invisible(NULL)
# }
#
# ####
#
# ## old code
# plot_strength <- function(x, ...) {
#
#   A <- x$adj
#
#   if (!is.null(x$varnames)) rownames(A) <- colnames(A) <- x$varnames
#   else rownames(A) = colnames(A) = as.character(1:x$p)
#
#   plot.new()
#   nAttrs <- list(label=rownames(A))
#   names(nAttrs$label) = rownames(A)
#   attrs <- list(node=list(shape="ellipse", fixedsize=FALSE))
#
#   g2 <- graphAM(A, c("directed"), values = list(weight=A))
#   eAttrs <- list()
#   eAttrs$color <- lapply(edgeWeights(g2), function(x) rgb(0,0,1,x))
#
#   Rgraphviz::plot(g2, nodeAttrs=nAttrs, edgeAttrs=eAttrs, attrs=attrs)
# }
#
