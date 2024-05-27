# Written originally by Devin Johnson

# Function to generate a prediction matrix for spline terms
# and re-parameterize for mixed model fitting
s2rPred <- function(sm, re, data) {
  X <- PredictMat(sm, data)   # Get prediction matrix for new data
  # Transform to random effect parameterization
  if (!is.null(re$trans.U)) X <- X %*% re$trans.U
  X <- t(t(X) * re$trans.D)
  # Re-order columns according to random effect re-ordering
  X[, re$rind] <- X[, re$pen.ind != 0] 
  # Re-order penalization index in the same way  
  pen.ind <- re$pen.ind; pen.ind[re$rind] <- pen.ind[pen.ind > 0]
  # Start return object
  r <- list(rand = list(), Xf = X[, which(re$pen.ind == 0), drop = FALSE])
  for (i in 1:length(re$rand)) { # Loop over random effect matrices
    r$rand[[i]] <- X[, which(pen.ind == i), drop = FALSE]
    attr(r$rand[[i]], "s.label") <- attr(re$rand[[i]], "s.label")
  }
  names(r$rand) <- names(re$rand)
  r
}
