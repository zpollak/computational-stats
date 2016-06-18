### PROBLEM 1 ###
# The logarithmic series random variable has pmf f_{alpha}(x) = (-(1 - alpha)^x)/(x * log(alpha)), 
# for x = [1:n] and alpha in (0,1). Write R function that takes two arguments, alpha & n,
# and returns n random variables each with density f_{alpha}. For alpha=0.3, use funciton
# to estimate probablity mass funciton f_{alpha} and fill table:
#   x_{i} =     |  1  |  2  |  3  |  4  |  5  |  6  |
#   true pmf =  |  ?  |  ?  |  ?  |  ?  |  ?  |  ?  |
#   est pmf  =  |  ?  |  ?  |  ?  |  ?  |  ?  |  ?  |
# ---------------------------------------------------------------------------------
n = 1000
alpha = .3
# True pmf values
Xtrue = rep(0, 6)
for (j in 1:6) {
  Xtrue[j] = (-(1 - alpha)^j) / (j * log(alpha))
}
Xtrue # true pmf

# Inversion function
set.seed(222)
logseries_inv = function(alpha, n) {
  X1a = rep(0,n)
  Xest = rep(0,n)
  for (j in 1:n) {
    U = runif(1);
    k = 1;
    S = ((-(1 - alpha)^k) / (k * log(alpha)));
    while (U>S) {
      k = k + 1;
      S = S + (-((1 - alpha) ^ k) / (k*log(alpha)));
      
    }
    X1a[j] = k
  }
  return(X1a)
}
# Compute frequencies of x's appearing and average --> f estimates
Xest = logseries_inv(alpha, n)
pmf.table = as.data.frame(table(Xest))[1:6,]
est.vals = pmf.table[,'Freq'] / n
pmf.table$f.true = Xtrue
pmf.table$f.est= est.vals
print(pmf.table[,-2])
#   Xest     f.true f.est
# 1    1 0.58140848 0.586
# 2    2 0.20349297 0.213
# 3    3 0.09496339 0.089
# 4    4 0.04985578 0.042
# 5    5 0.02791924 0.032
# 6    6 0.01628622 0.015


### PROBLEM 2 ###
# We consider the following density f(x) = {2/(pi(1+x)^2) if x>=0, else 0}. Describe an
# inversion method to sample from f. Implement the method in R. Write a funciton that
# takes an argument, n, and returns n RV's with density f. Use  this function to approximate
# the integral $f(x)dx from [0,2]. 
# ---------------------------------------------------------------------------------
# First, find the CDF: tan(pi*runif()/2)
# U ~ uniform(0,1)
# Define X = F^-1 ( U) 
f_inv2 = function(n) {
  U2 = runif(n)
  X2_1 = tan((pi / 2) * U2)
  return(X2_1)
}
X2_1 = f_inv2(n)
head(X2_1)
# 0.05602102 7.00491590 3.79598456 0.50100232 2.40318154 1.20281254

# Approximate the integral
mean(X2_1 <= 2) # est integral
# 0.694
# Find actual integral
integrand2 = function(x) {(2 / (pi * (1 + x^2)))}
integrate(integrand2, lower = 0, upper = 2) # true integral
# 0.7048328 with absolute error < 6.4e-12
atan(2) * 2 / pi
# 0.7048328 confirming the integrate function was properly used


### PROBLEM 3 ###
# Describe a rejection algorithm to sample from the density f given below. Using this
# algorithm, write an R function that returns n independent samples from f. Use function
# to approximate the integral $f(z)dz from [0.5,1].
#   f(x) = {(2/(1-e^(-1)))*x*e^(-x^2) if x in (0,1), else 0}
# ---------------------------------------------------------------------------------
# Pre-work:
#   f'(x) = (2/(1-e^(-1))) * [-2*x^2*e^(-x^2) + e^(-x^2)]
#   f' = (1-e^(-1))^(-1) * [2*(1-2*x^2)*e^(-x^2)]
#   f' = 0 = 1-2*x^2    -->    x = sqrt(1/2)
#   F(sqrt(1/2)) = M = (2/(1-e^(-1))) * sqrt(1/2) * e^(-1/2) ~= 1.357
#   U*M*g(Y) <= f(Y) --> U <= f(Y) / M = Y*exp(-Y^2)) / (sqrt(1/2) * exp(-1/2))
#     where U ~ uniform(0,1), Y ~ uniform(0,1)
# ---------------------------------------------------------------------------------
# Generate Y~U(0,1), U~U(0,1);
# If U <= (Y*exp(-Y^2)) / (sqrt(1/2) * exp(-1/2)) --> accept
# Else, reject and go back to step 1
f_reject3 = function(n) {
  Vy = numeric(n); 
  Vk = integer(n);
  i = 1; k = 0;
  while (i <= n){
    Y3 = runif(n);
    U3 = runif(n);
    k = k + 1;
    if (U3 <= ((Y3*exp(-(Y3)^2)) / (sqrt(1/2)*exp(-sqrt(1/2)^2)))) {
      Vy[i] = Y3; Vk[i] = k;
      i = i + 1; k = 0;
    }
  }
  return(list(Vy,Vk))
}
accept = f_reject3(n)[[1]]; # freject[[2]] = number of trials before each acceptance
# Plot histogram, add line
hist(accept, xlim = range(0, 1), prob = TRUE)
y = seq(0,1, length = n);
lines(y, (2/(1-exp(-1))) * y * exp(-y^2))
# Approximate integral
mean(accept >= 0.5)
# 0.654
integrand3 = function(x) {(2 / (1 - exp(-1))) * x * exp(-x^2)}
integrate(integrand3, lower = 1/2, upper = 1)
# 0.650068 with absolute error < 7.2e-15