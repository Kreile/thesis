################################################################################
##### FUNCTION FOR USING THE TEST OF EXCESS SIGNIFICANCE AS DESCRIBED IN   #####
##### IOANNIDIS AND TRIKALINOS (2007)                                      #####
##### Author: Robbie C.M. van Aert                                         #####
################################################################################

tes <- function(yi, sei, alpha, side) 
{
  
  ### FE meta-analysis for statistical power analysis
  est.fe <- rma(yi = yi, sei = sei, method = "FE")$b[1] 
  
  ### Compute statistical power and determine the number of observed 
  # statistically significant results
  if (side == "right") { 
    pow <- pnorm(qnorm(alpha, lower.tail = FALSE, sd = sei), mean = est.fe,
                 sd = sei, lower.tail = FALSE)
    O <- sum(pnorm(yi/sei, lower.tail = FALSE) < alpha)
  } else if (side == "left") {
    pow <- pnorm(qnorm(alpha, sd = sei), mean = est.fe, sd = sei)
    O <- sum(pnorm(yi/sei) < alpha)
  }
  
  E <- sum(pow) # Expected number of statistically significant result  
  n <- length(yi) # Number of studies in meta-analysis
  A <- (O - E)^2/E + (O - E)^2/(n - E) # Compute chi-square statistic
  pval.tes <- pchisq(A, 1, lower.tail = FALSE) # Compute p-value
  pval.tes <- ifelse(pval.tes < 0.5, pval.tes*2, (1-pval.tes)*2)
  
  return(data.frame(A = A, pval.tes = pval.tes, O = O, E = E, n = n))
  
}
