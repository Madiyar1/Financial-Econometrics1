###Financial Econometrics Project###

# Downloading dataset:
EUR.CAD <- Quandl("ECB/EURCAD",type = "xts")
# Convert the series into returns:
names(EUR.CAD)[1] <- "rates"
EUR_CAD <- diff(log(EUR.CAD$rates), lag = 1)
EUR_CAD <- na.omit(EUR_CAD)
# Plotting the log returns.
chartSeries(EUR_CAD, theme = "white")
# Basic Stats:
stats <- basicStats(EUR_CAD)
# Testing mean expected return.
stat_test <- t.test(EUR_CAD)
# Empirical Density Function of log returns
d1 = density(EUR_CAD)
plot(d1$x, d1$y, type="l") 
# ACF & PACF: for ARCH effect testing
Box.test(EUR_CAD , lag = 20, type = "Ljung")
t.test(EUR_CAD)
Box.test((EUR_CAD )**2, lag = 20, type = "Ljung")
ArchTest(EUR_CAD , 20)
# Building ARCH model:
acf(EUR_CAD , lag = 20)
pacf(EUR_CAD , lag = 20)
acf(EUR_CAD**2, lag = 20)
pacf(EUR_CAD**2, lag = 20)
# GARCH (1,1)
m1 = garchFit(~1+garch(1,1), data = EUR_CAD , trace = F)
summary(m1)
# Obtain ACF and PACF
v1 = volatility(m1)
resi = residuals(m1, standardize=T) # obtain residuals
acf(resi, lag = 20) 
pacf(resi, lag = 20)
acf(resi**2, lag = 20)
pacf(resi**2, lag = 20)
# Using t dstrb:
m2=garchFit(~1+garch(1,1), data = EUR_CAD , trace = F, cond.dist = "std")
summary(m2)
plot(m1) # QQ plot validity of the dstb asssumption
plot(m2)# QQ plot
resi = residuals(m2, standardize=T) # obtain residuals
acf(resi_0, lag = 20) 
pacf(resi_0, lag = 20)
acf(resi_0**2, lag = 20)
pacf(resi_0**2, lag = 20)
# Estimating GARCH(1,2);(2,1);(2,2)
m3=garchFit(~1+garch(2,2), data = EUR_CAD, trace = F, cond.dist = "std")
summary(m3)
plot(m3)
resi_1 = residuals(m3, standardize=T) # obtain residuals
acf(resi_1, lag = 20) 
pacf(resi_1, lag = 20)
acf(resi_1**2, lag = 20)
pacf(resi_1**2, lag = 20)
m4=garchFit(~1+garch(1,2), data = EUR_CAD, trace = F, cond.dist = "std")
summary(m4)
plot(m4)
resi_2 = residuals(m4, standardize=T) # obtain residuals
acf(resi_2, lag = 20) 
pacf(resi_2, lag = 20)
acf(resi_2**2, lag = 20)
pacf(resi_2**2, lag = 20)
m5=garchFit(~1+garch(2,1), data = EUR_CAD, trace = F, cond.dist = "std")
summary(m5)
plot(m5)
resi = residuals(m5, standardize=T) # obtain residuals
acf(resi_3, lag = 20) 
pacf(resi_3, lag = 20)
acf(resi_3**2, lag = 20)
pacf(resi_3**2, lag = 20)
# Estimating Aparch
m7=garchFit(~1+aparch(1,1,1), data = EUR_CAD, trace = F, cond.dist = "std")
summary(m7)
plot(m7)
resi_4 = residuals(m7, standardize=T) # obtain residuals
acf(resi_4, lag = 20) 
pacf(resi_4, lag = 20)
acf(resi_4**2, lag = 20)
pacf(resi_4**2, lag = 20)
m8=garchFit(~1+aparch(1,2), data = EUR_CAD, trace = F, cond.dist = "std")
summary(m8)
plot(m8)

#####


# IGARCH:
"Igarch" <- function(EUR_CAD,include.mean=F,volcnt=F){
  # Estimation of a Gaussian IGARCH(1,1) model.
  # EUR_CAD: return series 
  # include.mean: flag for the constant in the mean equation.
  # volcnt: flag for the constant term of the volatility equation.
  #### default is the RiskMetrics model
  #
  Idata <<- EUR_CAD
  Flag <<- c(include.mean,volcnt)
  #
  Mean=mean(Idata); Var = var(Idata); S = 1e-6
  if((volcnt)&&(include.mean)){
    params=c(mu = Mean,omega=0.1*Var,beta=0.85)
    lowerBounds = c(mu = -10*abs(Mean), omega= S**2, beta= S)
    upperBounds = c(mu = 10*abs(Mean), omega = 100*Var, beta = 1-S)
  }
  if((volcnt)&&(!include.mean)){
    params=c(omega=0.1*Var, beta=0.85)
    lowerBounds=c(omega=S**2,beta=S)
    upperBounds=c(omega=100*Var,beta=1-S)
  }
  #
  if((!volcnt)&&(include.mean)){
    params=c(mu = Mean, beta= 0.8)
    lowerBounds = c(mu = -10*abs(Mean), beta= S)
    upperBounds = c(mu = 10*abs(Mean), beta = 1-S)
  }
  if((!volcnt)&&(!include.mean)){
    params=c(beta=0.85)
    lowerBounds=c(beta=S)
    upperBounds=c(beta=1-S)
  }
  # Step 3: set conditional distribution function:
  igarchDist = function(z,hh){dnorm(x = z/hh)/hh}
  # Step 4: Compose log-likelihood function:
  igarchLLH = function(parm){
    include.mean=Flag[1]
    volcnt=Flag[2]
    mu=0; omega = 0
    if((include.mean)&&(volcnt)){
      my=parm[1]; omega=parm[2]; beta=parm[3]}
    if((!include.mean)&&(volcnt)){
      omega=parm[1];beta=parm[2]}
    if((!include.mean)&&(!volcnt))beta=parm[1]
    if((include.mean)&&(!volcnt)){mu=parm[1]; beta=parm[2]}
    #
    z = (Idata - mu); Meanz = mean(z**2)
    e= omega + (1-beta)* c(Meanz, z[-length(Idata)]**2)
    h = filter(e, beta, "r", init=Meanz)
    hh = sqrt(abs(h))
    llh = -sum(log(igarchDist(z, hh)))
    llh
  }
  # Step 5: Estimate Parameters and Compute Numerically Hessian:
  fit = nlminb(start = params, objective = igarchLLH,
               lower = lowerBounds, upper = upperBounds)
  ##lower = lowerBounds, upper = upperBounds, control = list(trace=3))
  epsilon = 0.0001 * fit$par
  cat("Estimates: ",fit$par,"\n")
  npar=length(params)
  Hessian = matrix(0, ncol = npar, nrow = npar)
  for (i in 1:npar) {
    for (j in 1:npar) {
      x1 = x2 = x3 = x4  = fit$par
      x1[i] = x1[i] + epsilon[i]; x1[j] = x1[j] + epsilon[j]
      x2[i] = x2[i] + epsilon[i]; x2[j] = x2[j] - epsilon[j]
      x3[i] = x3[i] - epsilon[i]; x3[j] = x3[j] + epsilon[j]
      x4[i] = x4[i] - epsilon[i]; x4[j] = x4[j] - epsilon[j]
      Hessian[i, j] = (igarchLLH(x1)-igarchLLH(x2)-igarchLLH(x3)+igarchLLH(x4))/
        (4*epsilon[i]*epsilon[j])
    }
  }
  cat("Maximized log-likehood: ",igarchLLH(fit$par),"\n")
  # Step 6: Create and Print Summary Report:
  se.coef = sqrt(diag(solve(Hessian)))
  tval = fit$par/se.coef
  matcoef = cbind(fit$par, se.coef, tval, 2*(1-pnorm(abs(tval))))
  dimnames(matcoef) = list(names(tval), c(" Estimate",
                                          " Std. Error", " t value", "Pr(>|t|)"))
  cat("\nCoefficient(s):\n")
  printCoefmat(matcoef, digits = 6, signif.stars = TRUE)
  
  if((include.mean)&&(volcnt)){
    mu=fit$par[1]; omega=fit$par[2]; beta = fit$par[3]
  }
  if((include.mean)&&(!volcnt)){
    mu = fit$par[1]; beta = fit$par[2]; omega = 0
  }
  if((!include.mean)&&(volcnt)){
    mu=0; omega=fit$par[1]; beta=fit$par[2]
  }
  if((!include.mean)&&(!volcnt)){
    mu=0; omega=0; beta=fit$par[1]
  }
  z=Idata-mu; Mz = mean(z**2)
  e= omega + (1-beta)*c(Mz,z[-length(z)]**2)
  h = filter(e,beta,"r",init=Mz)
  vol = sqrt(abs(h))
  
  Igarch <- list(par=fit$par,volatility = vol)
}
m9=Igarch(EUR_CAD)
# NGARCH:
"Ngarch" <- function(EUR_CAD){
  # Estimation of a non-symmertic GARCH, NGARCH(1,1), model. 
  # Assume normal innovations
  # EUR_CAD: return series 
  #
  # The likelihood function "glkn" can be modified to fit more general NGARCH 
  #  models.
  write(EUR_CAD,file='tmp.txt',ncol=1)
  # obtain initial estimates
  mu=mean(EUR_CAD)
  par=c(mu,0.01,0.8,0.01,0.7)
  #
  #
  mm=optim(par,glkn,method="Nelder-Mead",hessian=T)
  low=c(-10,0,0,0,0)
  upp=c(10,1,1,0.4,2)
  #mm=optim(par,glkn,method="L-BFGS-B",hessian=T,lower=low,upper=upp)
  ## Print the results
  par=mm$par
  H=mm$hessian
  Hi = solve(H)
  cat(" ","\n")
  cat("Estimation results of NGARCH(1,1) model:","\n")
  cat("estimates: ",par,"\n")
  se=sqrt(diag(Hi))
  cat("std.errors: ",se,"\n")
  tra=par/se
  cat("t-ratio: ",tra,"\n")
  # compute the volatility series and residuals
  ht=var(EUR_CAD)
  T=length(EUR_CAD)
  if(T > 40)ht=var(EUR_CAD[1:40])
  at=EUR_CAD-par[1]
  for (i in 2:T){
    sig2t=par[2]+par[3]*ht[i-1]+par[4]*(at[i-1]-par[5]*sqrt(ht[i-1]))^2
    ht=c(ht,sig2t)
  }
  sigma.t=sqrt(ht)
  Ngarch <- list(residuals=at,volatility=sigma.t)
}

glkn <- function(par){
  EUR_CAD=read.table("tmp.txt")[,1]
  glkn=0
  ht=var(EUR_CAD)
  T=length(EUR_CAD)
  if(T > 40)ht=var(EUR_CAD[1:40])
  at=EUR_CAD[1]-par[1]
  for (i in 2:T){
    ept=EUR_CAD[i]-par[1]
    at=c(at,ept)
    sig2t=par[2]+par[3]*ht[i-1]+par[4]*ht[i-1]*(at[i-1]/sqrt(ht[i-1])-par[5])**2
    ht=c(ht,sig2t)
    glkn=glkn + 0.5*(log(sig2t) + ept**2/sig2t)
  }
  glkn
}
m10=Ngarch(EUR_CAD)
### Data Wrangling




#######

#Importing Intra Daily Data

####### Data Wrangling:
EUR240 <- slice(EURCAD240.2, 1774:15657)
names(EUR240)[1] <- "Date"
names(EUR240)[2] <- "Time"
names(EUR240)[3] <- "Open"
names(EUR240)[4] <- "High"
names(EUR240)[5] <- "Low"
names(EUR240)[6] <- "Close"
names(EUR240)[7] <- "Volume"

###
OHLC <- EUR240 %>% 
  select(3:6)
OHLC <- as.matrix(OHLC)
my.array <- array(OHLC, dim = c(6,4,2314))

### Garman Klass realised volatility:
GarmanKlass = function(x){   
  # x is an array of prices with dimension:
  # (num.intervals.each.day)*( (5 or 6), "open", "high" etc)*(number.days) 
  n <- dim(x)[1] # number of intervals each day. 
  l <- dim(x)[3] # number of days 
  gk = NULL
  for (i in 1:l) {
    log.hi.lo <- log( x[1:n,2,i]/x[1:n,3,i] )
    log.cl.to.cl <- log( x[2:n,4,i]/x[1:(n-1),4,i] )
    gk[i] = ( sum(.5*log.hi.lo**2) - sum( (2*log(2) - 1)*(log.cl.to.cl**2) ) ) /n
  }
  return(gk)
}
###realised variance:
volat <- GarmanKlass(my.array)
###### Computing MSE MASE RMSE :




###### VaR:
p1=c(0.95,0.99)
quantile(EUR_CAD,p1)
qe=garchFit(~1+garch(1,1), data = EUR_CAD , trace = F, cond.dist = "std")
summary(qe)
vol=volatility(qe)
###Creating y:
y=EUR_CAD[2:2704]
y=cbind(y,vol[1:2703],abs(EUR_CAD[1:2703]))
colnames(y) <- c("EUR_CAD","vol1","abs1")
require(quantreg)
y=data.frame(y)
y=na.omit(y)
qe_1=rq(EUR_CAD~vol1+abs1,data=y,tau=0.95)
summary(qe_1)

####Change numbers here:
fit=0.00133+1.48123*vol[2704]-0.05159*abs(EUR_CAD[2704])
fit
m=rq(EUR_CAD~vol1+abs1,data=y,tau=0.99)
summary(mm)
fit=0.00832+1.87035*vol[2704]+0.25241*abs(EUR_CAD[2704])
fit





