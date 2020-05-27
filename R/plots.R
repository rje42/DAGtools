##' Plot graph from adjacency matrix
##'
##' @param amat adjacency matrix
##' @param p_thresh minimum value to plot
##' @param edgemode either \code{"directed"} or \code{"undirected"}
##'
##' @export
plot_graph <- function(amat, p_thresh=0.5, edgemode="directed") {
  # require(Rgraphviz)
  # require(graph)

  if (edgemode == "directed") A <- 1*(amat > p_thresh)
  else {
    A <- 1*(amat + t(amat) > p_thresh)
  }

  if (is.null(dimnames(A))) {
    varnames <- paste("x", seq_len(nrow(amat)), sep="")
    dimnames(A) <- list(varnames, varnames)
  }

  Rgraphviz::plot(graph::graphAM(adjMat=A, edgemode=edgemode))
}


##' Plot grid summarizing causal effects
##'
##' @param x object of class \code{aug_chain}
##' @param subset subset of variables to use
##' @param alpha significance level
##'
##' @details Given the output of \code{augment_chain},
##' this plots pairwise causal effect histograms in
##' a grid.
##'
##' @export
plot_causal_effects <- function(x, subset, alpha=0.05) {

  # require(ggplot2)
  # require(dplyr)
  # require(tidyr)

  CauEff <- x$CauEff
  dimnames(CauEff)[[1]] <- x$varnames_short
  dimnames(CauEff)[[2]] <- x$varnames_short
  # dimnames(CauEff) <- list(x$varnames,x$varnames,1:B)
  if (!missing(subset)) CauEff <- CauEff[subset,subset,]
  tmp <- as.data.frame.table(CauEff)

  p <- nrow(CauEff)

  # tmp <- tmp %>% filter(Var1 != Var2)

  ## determine which effects are significant at level alpha
  signif <- apply(CauEff, 1:2, function(x) mean(x < -1e-8) > 1-alpha || mean(x > 1e-8) > 1-alpha)
  dimnames(signif) = dimnames(CauEff)[1:2]
  signif <- as.data.frame.table(signif)
  names(signif)[3] <- "significant"

  # ACE <- apply(CauEff, 1:2, function(x) mean(x[abs(x) > 1e-8]))
  ACE <- apply(CauEff, 1:2, function(x) c(mean(x[abs(x) > 1e-8]), quantile(x[abs(x) > 1e-8], c(alpha/2,1-alpha/2))))
  ACE[is.na(ACE)] = 0
  names(dimnames(ACE)) <- c("stat", "Var1", "Var2")
  dimnames(ACE)[[1]][1:3] <- c("mn", "lq", "uq")

  ## matrix for annotations
  mns <- matrix(as.character(round(ACE[1,,], 2)), p, p)
  diag(mns) = ""
  # mns[upper.tri(mns) | lower.tri(mns)] <- ACE$caption[3*seq_len(p) - 2] # [p^2-p]
  # mns <- t(mns)
#  ann <- annotate("text", x=-.35, y=20*dim(CauEff)[3], label=c(t(mns)), cex=3)

  ACE <- as.data.frame.table(ACE)
  ACE <- ACE[ACE$Var1 != ACE$Var2,]
  names(ACE)[4] = "value"
  ## add text for caption
  ACE <- spread(ACE, "stat", "value") %>%
    mutate(caption=paste(round(mn,2))) %>% #, " (", round(lq,2), ", ", round(uq,2), ")", sep="")) %>%
    gather(stat, "value", c("mn", "lq", "uq"))

  ## put into useful format for ggplot
  tmp2 <- inner_join(inner_join(tmp, signif), ACE) %>% spread("stat", "value")
  tmp2 <- tmp2[abs(tmp2$Freq) > 1e-8, ]

  print(names(tmp2))

  ## create plot
#  pdf("pairwise.pdf", width=12, height=12)
  # filter(tmp2, Var1 %in% x$varnames[1:6], Var2 %in% x$varnames[1:6]) %>%
  ggplot(tmp2, aes(tmp2$Freq, ..count.., colour=factor(significant))) +
    geom_density() + # coord_cartesian(ylim=c(0, 2e5)) +
    geom_vline(xintercept = 0, color="red", lty=2) + facet_grid(Var1 ~ Var2) +
    scale_color_manual(values=c("FALSE"="black", "TRUE"="red")) +
    # annotate("text", x=-.35, y=20*dim(CauEff)[3], label=c(t(mns)), cex=3) +
    theme(legend.position="none")
 # dev.off()

}
