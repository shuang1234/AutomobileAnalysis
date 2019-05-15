# Read data
data <- read.csv("/auto-data.csv")

# Clean data
# Check missing values
summary(data)# six missing values in horsepower
# Delete rows with missings
data <- na.omit(data)
# Select variables needed
data <- data[, c("mpg", "displacement", "horsepower", "weight", "acceleration")]
# Attach data
attach(data)


# DEA
cor(data)
pairs(mpg~displacement + horsepower + weight + acceleration, 
      col="darkblue", pch=20, main="Scatterplot Matrices")




################################################################################################################
################################################################################################################
################################################################################################################
# Polynominal, Bspline, BinSoomth and highestAdjR2 functions from "P01-code-stats-functions" file
#-------------------------------------
#Polynomial Regression
#-------------------------------------
#The function polreg performs a polynomial regression with p order polynomials 
#Input Arguments:
#       x - vector containing the explonatory variable
#       y - vector containing the dependent variable
#       p - nth order polynomial 
#       ouptut - 1: delivers some output, 0: no output
#       opt   - 1: returns adj R-squared, 0: returns nothing
#       ploto - 1: Create new plot, 0: no new plot
polreg <- function(x,y,p, output=1, ploto=1, opt=0)
{
  #Scale data to [0,1]
  x<-x-min(x)
  x<-x/max(x)
  #Create a Design Matrix DM
  n <- length(x)
  q <- p + 1
  DM = matrix(1,n,q)
  DM[,2] <- x
  if(q>2) for(i in 3:q) { DM[,i] <- x**(i-1)}
  #Perform regression 
  reg <- lm(y~0+DM)
  
  #Calculate goodness of fit measures
  #Residual sum of squares
  rss <-  sum(sapply(residuals(reg), function(x) { x^2 }))
  #Coefficient of determination: R^2
  R2 <- 1 - (rss/ (t(y)%*%y-(mean(y)**2*n)))
  #Adjusted Coefficient of determination: R^2
  R2adj <- 1 - ( (n-1)/(n-q) ) * (1-R2)
  #AIC
  aic <- AIC(reg)
  
  if(output==1)
  {
    #Summary output  
    cat("Number of polyn.: ", p, "\n")
    cat("RSS: ", rss, "\n")
    cat("TSS: ", t(y)%*%y-(mean(y)**2/n), "\n")
    cat("R-squared: ", R2, "\n")
    cat("Adjusted R-squared: ", R2adj, "\n")
    cat("AIC: ", aic, "\n")
    #cat("Coefficents: \n")
    #print(coef(reg))
    #print(summary(reg))
    #print(anova(reg))
    
    #Graphic
    xp <- 0:100/100
    n <- length(xp)
    DMp = matrix(1,n,q)
    DMp[,2] <- xp
    if(q>2) for(i in 3:q) { DMp[,i] <- xp**(i-1)}
    
    if(ploto==1) par(mfrow=c(1,2))
    if(ploto==1) matplot(xp, DMp, type="l", lwd=2, main="Individual polynomial functions")
    if(ploto==1) plot(x,y, main="Polynomial regression", pch=20, col="darkblue")
    lines(xp, DMp%*%coef(reg), col="darkgreen", type="l", lwd=2)
  }
  
  if(opt==1) {return(c(R2adj,aic))}   
}

#-------------------------------------
#Bspline regression 
#The following code is heavily based on code of
#Dr. Samiran Sinha (Department of Statistics, Texas A&M)
#http://www.stat.tamu.edu/~sinha/research/note1.pdf 
#And Jeffrey S. Racine: A PRIMER ON REGRESSION SPLINES
#-------------------------------------
#The function bsplinereg performs a bspline regression with user defined knots
#Input Arguments:
#Input Arguments bsplinereg(...):
#       x - vector containing the explonatory variable
#       y - vector containing the dependent variable
#       knots - number of knots in [0,1]
#       ouptut - 1: delivers some output, 0: no output
#       opt   - 1: returns adj R-squared, 0: returns nothing
#       ploto - 1: Create new plot, 0: no new plot

#Calculate basis (rekursiv)
basis <- function(x, degree, i, knots) 
{ 
  
  if(degree == 0)
  { B <- ifelse((x >= knots[i]) & (x < knots[i+1]), 1, 0) 
  } else { 
    if((knots[degree+i] - knots[i]) == 0) 
    { alpha1 <- 0 
    } else { 
      alpha1 <- (x - knots[i])/(knots[degree+i] - knots[i]) } 
    
    if((knots[i+degree+1] - knots[i+1]) == 0) 
    { alpha2 <- 0 
    } else { alpha2 <- (knots[i+degree+1] - x)/(knots[i+degree+1] - knots[i+1]) } 
    B <- alpha1*basis(x, (degree-1), i, knots) + alpha2*basis(x, (degree-1), (i+1), knots) 
  } 
  
  return(B) 
}

#Create bspline Desin Matrix
bspline <- function(x, degree, knotpos) 
{ 
  #Number of basis
  K <- length(knotpos) + degree + 1 
  #Number of observations
  n <- length(x)
  #Set Boundary knots
  Boundary.knots = c(0,1)
  #create new vector with knot positons 
  knotpos <- c(rep(Boundary.knots[1], (degree+1)), knotpos, rep(Boundary.knots[2], (degree+1))) 
  
  
  #Create design matrix
  DM <- matrix(0,n,K) 
  for(j in 1:K) DM[,j] <- basis(x, degree, j, knotpos) 
  if(any(x == Boundary.knots[2])) DM[x == Boundary.knots[2], K] <- 1 
  #Return DM  
  return(DM) 
}


bsplinereg <- function(x, y, knots=0, knotsdef=NULL,  degree, output=1, ploto=1, opt=0)
{
  #Scale data to [0,1]
  x<-x-min(x)
  x<-x/max(x)
  #Sort x values in ascending order
  y <- y[order(x)]
  x <- sort(x)
  n <- length(x)
  
  #Calculate knot postions
  if(knots == 0) knotpos <- NULL
  if(knots != 0) knotpos <- 1:knots / (knots+1)
  if(length(knotsdef)>0) knotpos <- knotsdef
  #Create Design Matrix 
  DM <- bspline(x, degree, knotpos) 
  
  #Perform penalized regression
  reg <- lm(y ~ 0 + DM)
  print(summary(reg))
  
  #Calculate goodness of fit measures
  q <- length(knotpos) + degree + 1
  #Residual sum of squares
  rss <-  sum(sapply(residuals(reg), function(x) { x^2 }))
  #Coefficient of determination: R^2
  R2 <- 1 - (rss/ (t(y)%*%y-(mean(y)**2*n)))
  #Adjusted Coefficient of determination: R^2
  R2adj <- 1 - ( (n-1)/(n-q) ) * (1-R2)   
  #AIC
  aic <- AIC(reg)
  
  if(output==1)
  {
    #Summary output  
    cat("Number of knots = ", knots, "\n")
    cat("Knot positions = ", knotpos, "\n")
    cat("RSS: ", rss, "\n")
    cat("TSS: ", t(y)%*%y-(mean(y)**2/n), "\n")
    cat("R-squared: ", R2, "\n")
    cat("Adjusted R-squared: ", R2adj, "\n")
    cat("AIC: ", aic , "\n")
    #cat("Coefficents: \n")
    #print(coef(reg))
    #print(summary(reg))
    #print(anova(reg))
    
    #Graphics
    
    #Values for prediction
    xp <- 0:100/100
    DMp <- bspline(xp, degree, knotpos)
    
    if(ploto==1)par(mfrow=c(1,2))
    if(ploto==1) matplot(xp, (DMp), type="l", lwd=2, main="Individual spline functions")
    if(ploto==1) for(i in 1:length(knotpos)) abline(v=knotpos[i], col="red", lty=2) 
    if(ploto==1) plot(x,y, main="BSpline Regression", pch=20, col="darkblue")
    
    points(xp,DMp%*%coef(reg), type="l", lwd=2, col="brown")
    if(ploto==1) for(i in 1:length(knotpos)) abline(v=knotpos[i], col="red", lty=2) 
  }
  if(opt==1) return(c(R2adj, aic))       
}

#-------------------------------------
#Binsmooth 
#-------------------------------------
#The function binsmoothREG performs a binsmooth regression with a user defined binlength
#Input Arguments:
#       x - vector containing the explonatory variable
#       y - vector containing the dependent variable
#       binlength - amount of x values per bin
#       ouptut - 1: delivers some output, 0: no output
#       opt   - 1: returns adj R-squared, 0: returns nothing
#       ploto - 1: Create new plot, 0: no new plot



binsmoothREG <- function(x, y, binlength=0, knotsdef=NULL, output=1, ploto=1, opt=0)
{
  #Scale data to [0,1]
  x<-x-min(x)
  x<-x/max(x)
  #Sort x values in ascending order
  y <- y[order(x)]
  x <- sort(x)
  n <- length(x)
  #Devide data into bins
  if(is.vector(knotsdef)) bins = knotsdef
  else bins = ceiling(length(x) / binlength)
  #Create Design Matrix without intercept
  DM <- matrix(1,length(x),bins)
  #Set all elements not corresponding to region j equal 0
  for(i in 1:bins)
  {
    if(i==1) { xstart = 1 }
    if(i>1) { xstart = (i-1)*binlength+1 }
    xend = min(xstart + binlength-1, length(x))
    binelements <- xstart:xend
    elements <- 1:length(x)
    elements[binelements] <- 0
    DM[elements,i] <- 0
  }
  
  #Perform Linear Regreesion
  reg <- lm(y~0+DM)
  
  #Calculate goodness of fit measures
  q <- bins
  #Residual sum of squares
  rss <-  sum(sapply(residuals(reg), function(x) { x^2 }))
  #Coefficient of determination: R^2
  R2 <- 1 - (rss/ (t(y)%*%y-(mean(y)**2*n)))
  #Adjusted Coefficient of determination: R^2
  R2adj <- 1 - ( (n-1)/(n-q) ) * (1-R2)   
  #AIC
  aic <- AIC(reg)
  
  if(output==1)
  {
    #Summary output  
    cat("Elements per bin: ", binlength, "\n")
    cat("Number of bins: ", bins, "\n")
    cat("RSS: ", rss, "\n")
    cat("TSS: ", t(y)%*%y-(mean(y)**2/n), "\n")
    cat("R-squared: ", R2, "\n")
    cat("Adjusted R-squared: ", R2adj, "\n")
    cat("AIC: ", aic, "\n")
    #cat("Coefficents: \n")
    #print(coef(reg))
    #print(summary(reg))
    #print(anova(reg))
    
    #Graphic 
    if(ploto==1) plot(x,y, main="Binsmooth regression", pch=20, col="darkblue")
    j<-1
    for(i in 1:length(coef(reg)))
    {
      if(i>1) lines(c(x[xend],x[xend]), c(as.numeric(coef(reg)[i-1]), as.numeric(coef(reg)[i])), col="red", lwd=2)
      xstart = j
      if(i>1) lines(c(x[xend],x[xstart]), c(as.numeric(coef(reg)[i]), as.numeric(coef(reg)[i])), col="red", lwd=2)
      xend = min(j+binlength-1, length(x))
      lines(c(x[xstart],x[xend]), rep(as.numeric(coef(reg)[i]), 2), col="red", lwd=2)
      j<-j+binlength
      
    }
  }
  
  if(opt==1) return(c(R2adj,aic))    
}


#-------------------------------------
#Lowest AIC
#-------------------------------------
#The function lowestAIC is determining the highest adjusted R-squared and lowest AIC
#Input Arguments:
#Input Arguments bsplinereg(...):
#       x - vector containing the explonatory variable
#       y - vector containing the dependent variable
#       iter - number of iterations
#       FUN - Function which should be analyzed 
#       ...  optional FUN specif arguments
lowestAIC <- function(x, y, iter, FUN, ...)
{
  
  R2adj <- numeric(iter)
  aic <- numeric(iter)
  for(i in 1:iter) 
  {
    p<-i
    back <- FUN(x,y,p, output=0, opt=1, ...)
    
    if(length(back)==2) {R2adj[i]<-back[1]; aic[i]<-back[2]}
    if(length(back)==1) {R2adj[i]<-back}
  }
  #!!!!!!!!!!!!!!!The next two lines are not straight forward, but easy to avoid the problem of NaN 
  #One could simply delete NaN, but then the order of the array would be chaos
  R2adj[!complete.cases(R2adj)] <- max(complete.cases(R2adj))-0.1  #!!!!!!!!!!!!!!!
  if(length(back)==2) aic[!complete.cases(aic)] <- max(aic)+10 #!!!!!!!!!!!!!!!
  #!!!!!!!!!!!!!!!
  if(length(back)==2) par(mfrow=c(1,2))
  plot(1:iter,R2adj, type="l", main="Maximized adj. R-squared", col="darkblue", pch=20, lwd=2, xlab="i", ylab="adj. R-squared")
  if(length(back)==2)  plot(1:iter,aic, type="l", main="Minimized AIC", col="darkblue", pch=20, lwd=2, xlab="i", ylab="AIC")
  cat("Lowest AIC: ", min(aic), "\n")
  cat("i: ", which(aic == min(aic)), "\n")
  
  i <-  ((1:iter)[aic==min(aic)])
  minaic <- min(aic)
  return(minaic)
}


# Function I write to choose degree for B-splines
aic_deg <- function(n, exp_variable){
  deg <- numeric()
for(i in 1:n){
  deg[i] <- lowestAIC(x = exp_variable,y = mpg,iter = 30, bsplinereg, degree = i)
}
  plot(x = 1:n, y = deg, type="l", main="Minimized AIC", col="darkblue", pch=20, lwd=2, xlab="Degree", ylab="Minimised AIC")
}
  
#########################################################################################################################################
#########################################################################################################################################
#########################################################################################################################################

# Fit polynominal linear Model
# displacement 
lowestAIC(x = displacement,y = mpg, iter=30, polreg)
polreg(x = displacement,y = mpg, p = 10, ploto=1)
# Horsepower
lowestAIC(x = horsepower, y = mpg, iter = 30, polreg)
polreg(x = horsepower,y = mpg, p = 25, ploto=1)
# Weight
lowestAIC(x = weight,y = mpg, iter = 30, polreg)
polreg(x = weight,y = mpg, p = 2, ploto=1)
# Acceleration
lowestAIC(x = acceleration,y = mpg, iter = 30, polreg)
polreg(x = acceleration,y = mpg, p = 4, ploto=1)


# Fit B-splines
# displacement 
aic_deg(n = 5, exp_variable = displacement)
lowestAIC(x = displacement,y = mpg, iter=30, bsplinereg, degree  = 5)
bsplinereg(x = displacement,y = mpg, knots=5, knotsdef=c(), degree=5, ploto=1)

# Horsepower
aic_deg(n = 5, exp_variable = horsepower)
lowestAIC(x = horsepower, y = mpg, iter = 30, bsplinereg, degree  = 5)
bsplinereg(x = horsepower,y = mpg, knots=15, knotsdef=c(), degree=5, ploto=1)
# Weight
aic_deg(n = 5, exp_variable = weight)
lowestAIC(x = weight,y = mpg, iter = 30, bsplinereg, degree  = 1)
bsplinereg(x = weight,y = mpg, knots=1, knotsdef=c(), degree=1, ploto=1)
# Acceleration
aic_deg(n = 5, exp_variable = acceleration)
lowestAIC(x = acceleration,y = mpg, iter = 30, bsplinereg, degree  = 1)
bsplinereg(x = acceleration,y = mpg, knots= 4, knotsdef=c(), degree = 1, ploto=1)


# Fit BinSmooth
# Displacement
lowestAIC(x = displacement, y = mpg, iter=100, binsmoothREG)
lowestAIC(x = displacement, y = mpg, iter=20, binsmoothREG)
binsmoothREG(x =displacement, y = mpg, binlength = 5, ploto=1)

# Horesepower
lowestAIC(x = horsepower, y = mpg, iter=100, binsmoothREG)
lowestAIC(x = horsepower, y = mpg, iter=50, binsmoothREG)
lowestAIC(x = horsepower, y = mpg, iter=20, binsmoothREG)
lowestAIC(x = horsepower, y = mpg, iter=5, binsmoothREG)
binsmoothREG(x = horsepower, y = mpg, binlength = 4, ploto=1)

# Weight
lowestAIC(x = weight, y = mpg, iter=100, binsmoothREG)
lowestAIC(x = weight, y = mpg, iter=47, binsmoothREG)
binsmoothREG(x = weight, y = mpg, binlength = 47, ploto=1)

# Acceleration
lowestAIC(x = acceleration, y = mpg, iter=100, binsmoothREG)
lowestAIC(x = acceleration, y = mpg, iter=10, binsmoothREG)
binsmoothREG(x = acceleration, y = mpg, binlength = 3, ploto=1)













