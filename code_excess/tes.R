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
                 sd = sei, lower.tail = FALSE) #Probability to 
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



data <- data.ext2 %>% filter(meta.id == 168988)
168988 171006 171029 171032 171592 178208 178210
tes(ts)

tes <- function(data){
  yi = data$z
  sei = sqrt(data$var.z)
  alpha = 0.05
  
  ### FE meta-analysis for statistical power analysis
  est.fe <- rma(yi = yi, sei = sei, method = "FE")$b[1] 
  
  side <- ifelse(sign(est.fe) == -1, "left", "right")
  
  ### Compute statistical power and determine the number of observed 
  # statistically significant results
  if (side == "right") { 
    pow <- pnorm(qnorm(alpha, lower.tail = FALSE, sd = sei), mean = est.fe,
                 sd = sei, lower.tail = FALSE) #Probability to 
    O <- sum(pnorm(yi/sei, lower.tail = FALSE) < alpha)
  } else if (side == "left") {
    pow <- pnorm(qnorm(alpha, sd = sei), mean = est.fe, sd = sei)
    O <- sum(pnorm(yi/sei) < alpha)
  }
  
  E <- sum(pow) # Expected number of statistically significant result  
  n <- length(yi) # Number of studies in meta-analysis
  A <- (O - E)^2/E + (O - E)^2/(n - E) # Compute chi-square statistic
  pval.chi <- pchisq(A, 1, lower.tail = FALSE) # Compute p-value
  pval.chi <- ifelse(pval.chi < 0.5, pval.chi*2, (1-pval.chi)*2)
  
  pval.bin <- 1 - pbinom(q = O, n = n, prob = E/n)
  
  return(c(A = A, A.wald = A.wald, pval.chi = pval.chi, pval.bin = pval.bin, O = O, E = E, n = n))
  
}

O <- 3
E <- 1.06 # Expected number of statistically significant result  
n <- 20 # Number of studies in meta-analysis
A <- (O - E)^2/E + (O - E)^2/(n - E) # Compute chi-square statistic
pval.chi <- pchisq(A, 1, lower.tail = FALSE) # Compute p-value
pval.chi <- ifelse(pval.chi < 0.5, pval.chi*2, (1-pval.chi)*2)

est.p <- O/n

# pval.chi <- (1- pchisq(A, 1) )
pval.bin <- pbinom(q = O-1, size = n, prob = E/n, lower.tail = F)

c(A = A, A.wald = A.wald, pval.chi = pval.chi, pval.bin = pval.bin, O = O, E = E, n = n)






