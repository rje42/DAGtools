##' Implement Spearmen's rank for ordinal data
##'
##' @param x a numeric or ordinal vector
##' @param ... further numeric or ordinal vectors
##' @param use data to use
##'
##' @export
sprmn <- function (x, ..., use="pairwise") {

  if (is.data.frame(x)) {
    vars <- x
    nms <- names(x)
  }
  else {
    vars <- list(x, ...)
  }
  ln <- length(vars)

  if (ln < 2) return(NULL)

  ## check
  for (i in seq_along(vars)) {
    if (is.ordered(vars[[i]])) {
      vars[[i]] <- as.numeric(vars[[i]])
    }
    else if (!is.numeric(vars[[i]])) stop("All variables should be numeric or ordinal")
  }

  ## now compute ranks
  out <- diag(ln)
  if (exists("nms")) dimnames(out) <- list(nms, nms)

  for (i in seq_len(ln)[-1]) for (j in seq_len(i-1)) {
    out[i,j] <- out[j,i] <- cor(x=vars[[i]], y=vars[[j]], use=use, method="spearman")
  }

  return(out)
}
