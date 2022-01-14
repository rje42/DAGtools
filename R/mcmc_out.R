##' Set up MCMC Chain
##'
##' @param data_sets list of datasets
##' @param scoretype type of score: "bge" or "bde"
##' @param iterations number of iterations
##' @param ... other arguments to pass to \code{partitionMCMC}
##'
##' @details This will run \code{fit_mcmc} on multiple chains,
##' augment the output for each, and then join the outputs
##' together.  Intended for cases where we have multiple
##' imputations of the dataset.
##'
##'
##' @export
fit_multiple <- function(data_sets, scoretype="bge", iterations,  ...)  {

  n_data <- length(data_sets)
  partition_out <- list()

  ##
  for (i in seq_len(n_data)) {
    message(paste("Running MCMC on data set ", i, sep=""))
    # if (is.null(colnames(data_sets[[i]]))) colnames(data_sets[[i]]) <- paste0("X",seq_len(ncol(data_sets[[i]])))
    # myScore <- scoreparameters(scoretype=scoretype, data=data_sets[[i]])
    partition_out[[i]] <- fit_mcmc(data_sets[[i]], scoretype=scoretype, iterations=iterations, ...)
    # partition_out[[i]]$varnames <- names(data_sets[[i]])
  }

  for (i in seq_len(n_data)) {
    message(paste("Augmenting MCMC for chain ", i, sep=""))
    partition_out[[i]] <- augment_mcmc(partition_out[[i]], data_sets[[i]])
  }

  message("Merging chains")
  out <- merge_chain(partition_out)

  out
}

##' Set up MCMC Chain
##'
##' @param data dataset to fit
##' @param scoretype type of score: "bge" or "bde"
##' @param iterations number of iterations
##' @param ... other arguments to pass to \code{partitionMCMC}
##'
##' @details This is basically a wrapper for
##' BiDAG's \code{scoreparameters()} and \code{partitionMCMC()},
##' but which augments the chain with an array of sampled adjacencies, the
##' number of samples and variables, and the average adjacency.  Variable
##' names are also extracted from the original data set.
##'
##' @return An object of class
##' \code{summary_MCMCchain}.
##'
##'
##' @export
fit_mcmc <- function(data, scoretype="bge", iterations,  ...)  {

  myScore <- scoreparameters(scoretype=scoretype, data=data)
  out <- partitionMCMC(myScore, iterations=iterations, verbose=FALSE, ...)
  out$varnames <- names(data)

  if (!is.null(out)) {
    B <- out$info$samplesteps
    p <- nrow(out$DAG)

    ## get adjacencies
    adj <- unlist(out$traceadd$incidence)
    dim(adj) <- c(p,p,B)

    meanAdj <- .rowMeans(adj, p^2, B)
    dim(meanAdj) <- c(p,p)
    rownames(meanAdj) = colnames(meanAdj) = out$varnames
    sparsity <- sum(meanAdj)/choose(p,2)

    out <- c(out, list(adj=meanAdj, sparsity=sparsity, B=B, p=p))
    out$DAGscores <- unlist(out$trace)
    out$partitionscores <- unlist(out$traceadd$partitionscores)
    out <- c(out, out$traceadd)
    out <- out[names(out) != "traceadd"]

  }

  class(out) <- "summary_MCMCchain"

  out
}

# ##' Summarise output from MCMC
# ##'
# ##' @param object output from \code{BiDAG} package
# ##' @param ... other arguments
# ##'
# ##' @details Constructs the average adjacency matrix
# ##' and mean sparsity level.  Replaces list of DAG and partition
# ##' scores with vectors.  Also has a print and plot
# ##' method.
# ##'
# ##' @export
# summary.MCMCtrace <- function(object, ...) {
#   B <- length(object$incidence)
#   p <- nrow(object$incidence[[1]])
#
#   adj <- unlist(object$incidence)
#   dim(adj) <- c(p,p,B)
#
#   meanAdj <- .rowMeans(adj, p^2, B)
#   dim(meanAdj) <- c(p,p)
#   rownames(meanAdj) = colnames(meanAdj) = object$varnames
#   sparsity <- sum(meanAdj)/choose(p,2)
#
#   out <- c(object, list(adj=meanAdj, sparsity=sparsity, B=B, p=p))
#   out$DAGscores <- unlist(out$DAGscores)
#   out$partitionscores <- unlist(out$partitionscores)
#   class(out) <- "summary_MCMCchain"
#
#   out
# }

##' @export
print.summary_MCMCchain <- function(x, digits = max(2L, getOption("digits") - 4L), ...) {
  cat("Partition MCMC with ",x$B," iterations over ",x$p," variables\n", sep="")
  cat("Mean edge density: ", signif(x$sparsity, digits), "\n", sep="")
  print.default(round(x$adj, digits), ...)

  invisible(x)
}

##' @export
print.aug_chain <- function(x, ...) {
  print.summary_MCMCchain(x, ...)
}

##' Plot average adjacency from MCMC output
##'
##' @param x object of class \code{summary.MCMCchain}
##' @param p_thresh threshold for plotting an edge (default is 0.5)
##' @param ... other arguments
##'
##' @details This plots an undirected and directed summary
##' of the output.  The method will only plot an edge if it appears
##' in at least a proportion \code{p_thresh} of the samples.
##' Consequently, an undirected edge may appear but no
##' corresponding directed edge.
##'
##'
##' @export
plot.summary_MCMCchain <- function(x, p_thresh=0.5, ...) {

  undAdj <- 1*(x$adj + t(x$adj) > p_thresh)
  if (!is.null(x$varnames)) rownames(undAdj) = colnames(undAdj) = x$varnames
  else rownames(undAdj) = colnames(undAdj) = as.character(1:x$p)

  # plot.new()
  par(mfrow=c(1,2))
  nAttrs <- list(label=rownames(undAdj))
  names(nAttrs$label) = rownames(undAdj)
  attrs <- list(node=list(shape="ellipse", fixedsize=FALSE))
  plot1 <- Rgraphviz::plot(graphAM(undAdj), nodeAttrs=nAttrs, attrs=attrs)

  node_pos <- transpose(getNodeXY(plot1))
  names(node_pos) <- rownames(undAdj)

  w_adj <- x$adj
  w_adj[w_adj < p_thresh] <- 0
  # rownames(w_adj) = colnames(w_adj) = rownames(undAdj)
  g2 <- graphAM(w_adj, c("directed"), values = list(weight=w_adj))

  # stop()
  # print(edges(g2))
  # eAttrs <- list(label=rownames(undAdj))
  # nAttrs <- list(label=rownames(undAdj), pos=getNodeXY(plot1))  ## giving position doesn't seem to work
  graph::nodes(g2) <- rownames(undAdj)
  Rgraphviz::plot(g2, nodeAttrs=nAttrs, attrs=attrs)

  invisible(x)
}

##' @export
plot.aug_chain <- function(x, ...) {
  plot.summary_MCMCchain(x, ...)
}
