
#' Non-parametric Gaussian kernel density estimator
f_np_pdf <- function(data, h = NULL, n_mesh, lower, upper, type = 1) {
  
  if (is.null(h)) {
    #h  <- stats::bw.nrd(data) # silverman's
    h <- 1.06 * (length(data)^(-1/5)) * sd(data) # silverman's Laurent
  }
  tmp <- stats::density(data, bw = h, adjust = 1, kernel = "gaussian", 
                        from = lower, to = upper, n = n_mesh)
  x <- as.numeric(tmp$x)
  y <- as.numeric(tmp$y)
  y[x < min(data) - h] <- 0.0
  y[x > max(data) + h] <- 0.0
  
  # ensure sum to 1
  try(y <- f_norm_pdf(x, y, type = type))
  
  out <- list("x"   = x, 
              "pdf" = y, 
              "h"   = h)
  out
}

#' Function which normalizes a PDF
f_norm_pdf <- function(x, pdf, type = 1) {
      pdf[pdf < 1e-8] <- 0
      dens <- f_integrate(x, pdf, type = type)
      pdf  <- pdf / dens
      pdf
}

#' Function for numerical integration 
#' Used as a wrapper for better flexibility
f_integrate <- function(x, y, type = 1) {
      
      if (type == 1) {
            out <- pracma::trapz(x, y)
      }
      if (type == 2) {
            f.int <- stats::approxfun(x = x, y = y, method = "linear", yleft = 0, yright = 0, rule = 2, f = 0)
            out   <- stats::integrate(f = f.int, lower = min(x), upper = max(x), subdivisions = 1000L)$value
      }
      out <- as.numeric(out)
      out
}



