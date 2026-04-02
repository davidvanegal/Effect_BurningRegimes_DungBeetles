# ==============================================================================
# DetAbu.R
# Detected species relative abundance estimation
# Part of JADE (Joint species-rank Abundance Distribution/Estimation)
# Source: Chao et al. (2015) — iNEXT framework
# ==============================================================================
# Description:
#   Estimates the relative abundance of detected species based on
#   individual-based (abundance) data.
#
# Usage:
#   DetAbu(x, zero = FALSE)
#
# Arguments:
#   x     - A numeric vector of species abundance frequencies
#   zero  - Logical. Whether to retain zero-frequency species. Default: FALSE
#
# Returns:
#   A numeric vector of estimated relative abundances
# ==============================================================================

DetAbu <- function(x, zero = FALSE) {
  x  <- unlist(x)
  n  <- sum(x)
  f1 <- sum(x == 1)
  f2 <- sum(x == 2)
  f3 <- sum(x == 3)

  # Handle edge case: no doubletons
  if (f2 == 0) {
    f1 <- max(f1 - 1, 0)
    f2 <- 1
  }

  # Probability adjustment coefficients
  A1 <- f1 / n * ((n - 1) * f1 / ((n - 1) * f1 + 2 * max(f2, 1)))
  A2 <- f2 / choose(n, 2) * ((n - 2) * f2 / ((n - 2) * f2 + 3 * max(f3, 1)))^2

  # Remove zero-abundance species unless requested
  if (zero == FALSE) x <- x[x > 0]

  # Solve for detection probability adjustment factor q
  q.solve <- function(q) {
    e   <- A1 / sum(x / n * exp(-q * x))
    out <- sum((x / n * (1 - e * exp(-q * x)))^2) -
           sum(choose(x, 2) / choose(n, 2)) + A2
    abs(out)
  }

  q <- tryCatch(
    optimize(q.solve, c(0, 1))$min,
    error = function(e) { 1 }
  )

  # Compute adjusted relative abundances
  e <- A1 / sum(x / n * exp(-q * x))
  o <- x / n * (1 - e * exp(-q * x))
  o
}
