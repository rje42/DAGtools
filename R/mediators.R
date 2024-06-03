##' Standard summary to use
##'
##' @param x numerical vector
##' @param alpha one minus the coverage
##'
my_sum <- function(x, alpha=0.1) {
  c(mean(x, na.rm = TRUE), quantile(x, c(alpha/2, 1-alpha/2), na.rm = TRUE), mean(abs(x) > 1e-8, na.rm = TRUE))
}


##' Summarize causal relationship between two variables
##'
##' @param object augmented output of MCMC algorithm
##' @param i,j variables whose relationship is to be analysed
##' @param digits how to round numbers in printed output (see details)
##' @param alpha mass to remove from confidence interval (so defaults to 90\% CI)
##' @param show logical: should ACE be printed?
##' @param indirect logical: if true includes indirect effects
##'
##' @details Get summary of causal connections and mediations.
##' Prints summary and returns list of \code{total_eff},
##' \code{direct_eff} and \code{indirect_eff}, each being a list
##' with entries \code{to} and \code{from}.
##'
##' The output is printed as a mixture of percentages and
##' causal effects, each of which is rounded to 2 significant
##' figures by default.  This value can be modified by
##' using the \code{digits} argument.
##'
##' @export
sum_causal_paths <- function(object, i, j, digits=2, alpha=0.1, show=TRUE, indirect=TRUE) {

  betas <- object$betas
  CauEff <- object$CauEff
  p <- nrow(betas)
  B <- dim(betas)[3]

  nms <- object$varnames
  if (is.null(nms)) nms <- paste("x", 1:p, sep="")

  wh_to <- (abs(CauEff[i,j,]) >= 1e-8)
  wh_from <- (abs(CauEff[j,i,]) >= 1e-8)

  total_eff <- my_sum(CauEff[i,j,wh_to], alpha=alpha)
  total_eff_r <- my_sum(CauEff[j,i,wh_from], alpha=alpha)
  total_eff[4] <- mean(wh_to)
  total_eff_r[4] <- mean(wh_from)

  if (show) {
    cat(signif(100*mean(wh_to), digits), "% instances with ", nms[i], " causing ", nms[j], "\n", sep="")
    cat("Average Causal Effect ", round(total_eff[1],digits),
        " (", round(total_eff[2],digits), ", ", round(total_eff[3],digits), ")", "\n\n", sep="")
  }

  dir_eff <- my_sum(betas[i,j,wh_to], alpha=alpha)
  dir_eff_r <- my_sum(betas[j,i,wh_from], alpha=alpha)

  if (show) {
    cat("Direct effect present ", signif(100*dir_eff[4],digits), "% of those cases\n", sep="")
    cat("Average Direct Effect ", round(dir_eff[1],digits),
        " (", round(dir_eff[2],digits), ", ", round(dir_eff[3],digits), ")", "\n\n", sep="")
  }

  if (indirect) {
    ## indirect effects
    indir_eff <- indir_eff_r <- matrix(0, p, 4)
    rownames(indir_eff) <- rownames(indir_eff_r) <- nms
    colnames(indir_eff) <- colnames(indir_eff_r) <- c("mean", "lq", "uq", "prop")

    for (v in seq_len(p)) {
      if (v==i || v==j) next
      tmp <- CauEff[i,v,wh_to]*CauEff[v,j,wh_to]
      tmp_r <- CauEff[j,v,wh_from]*CauEff[v,i,wh_from]
      indir_eff[v,] <- my_sum(tmp, alpha=alpha)
      indir_eff_r[v,] <- my_sum(tmp_r, alpha=alpha)
    }

    indir_eff <- indir_eff[-c(i,j),,drop=FALSE]
    indir_eff <- indir_eff[order(-abs(indir_eff[,1])),,drop=FALSE]
    indir_eff_r <- indir_eff_r[-c(i,j),,drop=FALSE]
    indir_eff_r <- indir_eff_r[order(-abs(indir_eff_r[,1])),,drop=FALSE]

    if (show) {
      for (v in seq_len(nrow(indir_eff))) {
        if (!is.na(indir_eff[v,4]) && indir_eff[v,4] > 0) {
          cat(nms[i], " --> ", rownames(indir_eff)[v], " --> ", nms[j], "  :  (")
          cat(signif(100*indir_eff[v,4], digits), "%)  ", round(indir_eff[v,1], digits), " (", round(indir_eff[v,2],digits), ", ", round(indir_eff[v,3],digits), ")\n", sep="")
        }
      }
      cat("\n")
    }
  }

  if (show) {
    cat(signif(100*mean(wh_from),digits), "% instances with ", nms[j], " causing ", nms[i], "\n", sep="")
    cat("Average Causal Effect ", round(total_eff_r[1],digits),
        " (", round(total_eff_r[2],digits), ", ", round(total_eff_r[3],digits), ")", "\n\n", sep="")

    cat("Direct effect present ", signif(100*dir_eff_r[4],digits), "% of those cases\n", sep="")
    cat("Average Direct Effect ", round(dir_eff_r[1],digits),
        " (", round(dir_eff_r[2],digits), ", ", round(dir_eff_r[3],digits), ")", "\n\n", sep="")

    if (indirect) {
      for (v in seq_len(nrow(indir_eff_r))) {
        if (!is.na(indir_eff_r[v,4]) && indir_eff_r[v,4] > 0) {
          cat(nms[j], " --> ", rownames(indir_eff_r)[v], " --> ", nms[i], "  :  (")
          cat(signif(100*indir_eff_r[v,4],digits), "%)  ", round(indir_eff_r[v,1],digits), " (", round(indir_eff_r[v,2],digits), ", ", round(indir_eff_r[v,3],digits), ")\n", sep="")
        }
      }
    }
    cat("\n")
  }

  if (indirect) invisible(list(total_eff=list(to=total_eff, from=total_eff_r), dir_eff=list(to=dir_eff, from=dir_eff_r), indir_eff=list(to=indir_eff, from=indir_eff_r)))
  else invisible(list(total_eff=list(to=total_eff, from=total_eff_r), dir_eff=list(to=dir_eff, from=dir_eff_r), indir_eff=NULL))
}

