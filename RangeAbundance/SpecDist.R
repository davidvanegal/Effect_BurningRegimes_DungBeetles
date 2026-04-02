# ==============================================================================
# SpecDist.R
# Species rank-abundance distribution estimation
# Part of JADE (Joint species-rank Abundance Distribution/Estimation)
# Source: Chao et al. (2015) — iNEXT framework
# ==============================================================================
# Description:
#   Estimates species rank-abundance (RAD) or rank-incidence (RID) distributions
#   for both detected and undetected species, based on abundance or incidence data.
#
# Dependencies:
#   DetAbu.R  — detected species relative abundance (abundance data)
#   UndAbu.R  — undetected species relative abundance (abundance data)
#
# Usage:
#   SpecDist(x, datatype = "abundance")
#
# Arguments:
#   x        - Numeric vector of species abundances or incidence frequencies.
#              For incidence data, the first entry must be the total number of
#              sampling units, followed by per-species incidence frequencies.
#   datatype - Character. Either "abundance" or "incidence". Default: "abundance"
#
# Returns:
#   A data.frame with columns:
#     probability  - Estimated relative abundance or incidence probability
#     method       - "detected" or "undetected"
#   Rows are sorted in descending order of probability.
# ==============================================================================

SpecDist <- function(x, datatype = "abundance") {

  # Validate datatype argument
  TYPE <- c("abundance", "incidence")
  if (is.na(pmatch(datatype, TYPE)))    stop("invalid datatype")
  if (pmatch(datatype, TYPE) == -1)     stop("ambiguous datatype")
  datatype <- match.arg(datatype, TYPE)

  if (datatype == "abundance") {
    out <- rbind(
      data.frame("probability" = DetAbu(x, zero = TRUE), "method" = "detected"),
      data.frame("probability" = UndAbu(x),               "method" = "undetected")
    )
  } else if (datatype == "incidence") {
    out <- rbind(
      data.frame("probability" = DetInc(x, zero = TRUE), "method" = "detected"),
      data.frame("probability" = UndInc(x),               "method" = "undetected")
    )
  }

  # Return sorted by descending probability
  out[order(-out[, 1]), ]
}
