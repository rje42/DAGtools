##' Generate Tikz Code for Plot
##'
##' @param x adjacency matrix from MCMC summary
##' @param file file to output code (defaults to stout)
##' @param z optionally, a matrix of z-scores
##' @param n_pos optional matrix of node positions (requires Rgraphviz if not specified)
##' @param ud_thresh threshold for drawing edges as undirected
##' @param cutoff threshold for including edges
##' @param bend edge angles
##' @param func vectorized function to transform probabilities
##' @param bend_double bend automatically if arrows in two different directions?
##' @param rev logical - should thickness and density plotting be reversed?
##' @param col colours for undirected, directed and positive and negative edges
##' @param labels character vector of labels for nodes
##' @param nodesOnTop logical: should nodes be on top of everything else
##' @param edgeScale scales edges to a desired width
##' @param ... other arguments
##'
##' @details If more than \code{ud_thresh} of the edges present
##' are in a particular direction, then the two directions are drawn
##' separately, otherwise an undirected edge is drawn.  Edges present
##' in less than \code{cutoff} sampled graphs are ignored entirely. \code{func}
##' is a transformation that can be applied to the probability (e.g.
##' squaring it to reduce the visual impact of weak edges).  This must
##' take a matrix and return a matrix of the transformed values.
##'
##' If \code{z} is specified, then the function will plot the thickness
##' of the line to be related to this quantity.
##' If \code{rev = TRUE} then the thickness of the line is related to the
##' size of the effect, and the density of colour to the proportion of
##' samples containing that effect.
##'
##' \code{col} is given as a character vector of up to four colours.
##' The order is undirected positive and negative, followed by
##' directed positive and negative.  Colours are just repeated if
##' necessary.
##'
##' If \code{labels} is not provided these are taken from the
##' rownames of \code{x}.
##'
##' @export
TikzCode <- function(x, file="", z, n_pos, bend,
                     ud_thresh=0.5, cutoff=0.1, func,
                     bend_double=TRUE, rev=FALSE,
                     col=c("black","black","blue","red"), labels=NULL,
                     nodesOnTop=TRUE, edgeScale = 1, ...) {
  if (missing(func)) func <- function(x) x
  p <- nrow(x)
  if (p == 0) invisible(NULL)

  if (length(col < 4)) col <- rep_len(col, 4)

  ## determine whether edge should be undirected
  # ud_edge <- (x > ud_thresh) & (t(x) > ud_thresh)
  ratio <- x/(x+t(x))
  ratio <- pmin(ratio, 1-ratio)
  ratio[is.na(ratio)] = FALSE
  ud_edge <- (ratio > 1-ud_thresh)

  ## get values for color intensity and width
  if (!rev) {
    A <- func(x)
    A[ud_edge] <- A[ud_edge]*(abs(A)[ud_edge] + t(abs(A))[ud_edge] > cutoff)
    A[!ud_edge] <- A[!ud_edge]*(abs(A)[!ud_edge] > cutoff)

    if (!missing(z)) {
      W <- abs(z)
    }

  }
  else {
    W <- func(x)
    W[ud_edge] <- W[ud_edge]*(abs(W)[ud_edge] + t(abs(W))[ud_edge] > cutoff)
    W[!ud_edge] <- W[!ud_edge]*(abs(W)[!ud_edge] > cutoff)

    if (!missing(z)) {
      A <- abs(z)
    }
  }

  ## if missing n_pos, then use Rgraphviz to retrieve them
  if (missing(n_pos)) {
    tmp <- (x > 0.5) %>% graphAM(edgemode = "directed") %>% agopen(name = "blah") %>% graphLayout
    n_pos <- tmp %>% getNodeXY %>% transpose %>% unlist %>% matrix(nrow=2) %>% t.default / 100
    n_pos[,1] <- n_pos[,1] - mean(n_pos[,1])
    n_pos[,2] <- n_pos[,2] - mean(n_pos[,2])
  }

  ## matrix of bend commands
  if (missing(bend)) {
    bend <- matrix(0, p, p)
    if (bend_double) bend[!ud_edge & x > cutoff & t(x) > cutoff] = 10
  }
  else if (bend_double) {
    bend[bend == 0 & !ud_edge & x > cutoff & t(x) > cutoff] = 10
  }
  # bend_text <- matrix("", p, p)
  bend_text = paste("[bend left=", bend, "]", sep="")
  dim(bend_text) <- dim(bend)
  # bend_text[bend < 0] = paste("[bend right=", bend[bend < 0], "]", sep="")

  ## variable names
  if (missing(labels)) nms <- rownames(x)
  else nms <- labels

  ## start a new file
  cat("", file=file)

  ## node commands
  if (!nodesOnTop) {
    for (i in 1:p) {
      cat("\\node[rv] (", i, ") at (", n_pos[i,1], ",", n_pos[i,2], ") {", nms[i], "};\n", sep="", file=file, append = TRUE)
    }
  }
  else {
    for (i in 1:p) {
      cat("\\node[rv0] (", i, ") at (", n_pos[i,1], ",", n_pos[i,2], ") {};\n", sep="", file=file, append = TRUE)
    }
  }
  cat("\n", file=file, append = TRUE)

  ## edge commands
  for (i in seq_len(p)) for (j in seq_len(p)[-i]) {
    if (ud_edge[i,j] && A[i,j]+A[j,i] > cutoff) {
      ## draw undirected edge (once only)
      if (i < j && missing(z)) cat("\\draw[ua, ", col[1],", opacity=", (A[i,j]+A[j,i]), "] (", i, ") to", bend_text[i,j], " (", j, ");\n", sep="", file=file, append = TRUE)
      else if (i < j && W[i,j]+W[j,i] > 0) cat("\\draw[ua, ", col[1],", line width = ", edgeScale*(W[i,j]+W[j,i]), "mm, opacity=", (A[i,j]+A[j,i]), "] (", i, ") to", bend_text[i,j], " (", j, ");\n", sep="", file=file, append = TRUE)
    }
    else if (ud_edge[i,j] && A[i,j]+A[j,i] < -cutoff) {
      ## draw undirected edge (once only)
      if (i < j && missing(z)) cat("\\draw[ua, ", col[2],", opacity=", abs(A[i,j]+A[j,i]), "] (", i, ") to", bend_text[i,j], " (", j, ");\n", sep="", file=file, append = TRUE)
      else if (i < j && W[i,j]+W[j,i] > 0) cat("\\draw[ua, ", col[2],", line width = ", edgeScale*(abs(W[i,j]+W[j,i])), "mm, opacity=", abs(A[i,j]+A[j,i]), "] (", i, ") to", bend_text[i,j], " (", j, ");\n", sep="", file=file, append = TRUE)
    }
    else if (A[i,j] > cutoff) {
      ## draw directed edge
      if (missing(z)) cat("\\draw[da, ", col[3],", opacity=", A[i,j], "] (", i, ") to", bend_text[i,j], " (", j, ");\n", sep="", file=file, append = TRUE)
      else if (W[i,j] > 0) cat("\\draw[da, ", col[3],", line width = ", edgeScale*W[i,j], "mm, opacity=", A[i,j], "] (", i, ") to", bend_text[i,j], " (", j, ");\n", sep="", file=file, append = TRUE)
    }
    else if (A[i,j] < -cutoff) {
      ## draw directed edge
      if (missing(z)) cat("\\draw[da, ", col[4],", opacity=", abs(A[i,j]), "] (", i, ") to", bend_text[i,j], " (", j, ");\n", sep="", file=file, append = TRUE)
      else if (W[i,j] > 0) cat("\\draw[da, ", col[4],", line width = ", edgeScale*W[i,j], "mm, opacity=", abs(A[i,j]), "] (", i, ") to", bend_text[i,j], " (", j, ");\n", sep="", file=file, append = TRUE)
    }
  }

  if (nodesOnTop) {
    for (i in 1:p) {
      cat("\\node[rv] at (", i, ") {", nms[i], "};\n", sep="", file=file, append = TRUE)
    }
    cat("\n", file=file, append = TRUE)
  }

  invisible(NULL)
}


