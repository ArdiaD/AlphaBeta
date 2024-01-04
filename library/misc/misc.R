
f_idx_strat <- function(strat, strat.i) {
      if (strat.i == "ALL") {
            idx <- rep(TRUE, length(strat))
      } else {
            idx <- strat == strat.i
      }
      idx
}
