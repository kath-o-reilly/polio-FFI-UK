#elimination functions

# based on AR_2/AR_1 = R_2/R_1 ... etc and (AR_1 x PrP_1) + (AR_2 x PrP_2) + ... + (AR_n x PrP_n) = 1
# we can solve for AR_1 ... AR_n

# inputs
# n = number of RR (length(RR))
# RR = values of RR
# x is AR

# ARfunct <- function(x,data) {
#   RR <- data$RR          # estimated risk
#   PrP <- data$PrP
#   n <- length(RR)                          # number of equations
#   f <- numeric(n) 					               # read as:
#   tmp <- numeric(n)                        # tmp for final equation
#   for(i in 1:(n-1)){
#     f[i] <-  RR[i+1]/RR[1] - x[i+1]/x[1]   # AR_2/AR_1 - R_2/R_1 = 0 etc...
#     tmp[i] <- x[i]*PrP[i]
#   }
#   tmp[n] <- x[n]*PrP[n]
#   f[n] <- sum(tmp) - 1                     # (AR_1 x PrP_1) + (AR_2 x PrP_2) + ... + (AR_n x PrP_n) = 1
#   #f
#   #browser()
#   return(abs(sum(f)))
# } 

ARfunct2 <- function(data) {
  # add in some check
  if(sum(data$PrP)!=1){
    warning("check sum of PrP")
  }
  if(sum(data$RR==1)<1){
    warning("check that 1 RR == 1.00")
  }
  RR <- data$RR          # estimated risk
  PrP <- data$PrP
  n <- length(RR)                          # number of equations
  # which RR==1
  small <- which(RR==1)
  AR <- rep(NA,n)
  # What is AR_small? ie. eqn 1
  AR_tmp <- rep(NA,n-length(small))
  tmp1 <- sum(PrP[small])
  tmp2 <- RR[-small]*PrP[-small]
  AR1 <- 1/(tmp1 + sum(tmp2))
  # eqn 2 - use AR1 in eqns
  vals <- c(1:n)
  vals2 <- vals[-small]
  AR[small] <- AR1  # those with RR=1
  for(i in vals2){
    AR[i] <- RR[i]*AR1  # those left
  }
  #browser()
  return(AR)
} 

# check using example
# n <- 3
# RR <- rep(1,n) #c(3,2,1)
# RR[2:3] <- c(1.5,1.6)
# PrP <- c(0.8,0.1,0.1) #rep(1/n,n)   #c(0.1,0.1,0.8)
# data <- data.frame(RR=RR,PrP=PrP)
# 
# ARfunct2(data)
# 
# tmp <- RR
# tmp[tmp==1] <- 0.5
# out <- optim(par = tmp, fn = ARfunct, data = data,
#                         method="L-BFGS-B",lower=rep(0.0001,n),upper=rep(1e+10,n))#,
# # #control = list(maxit = 200000))  # answers["x"] = x,y,z are the solutions closest to startx if there are multiple solutions
# sum(out$par*PrP)
# 

# end
