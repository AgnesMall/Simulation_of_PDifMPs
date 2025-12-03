#############################################################################
#
# Simulation of paths of the four test problems
# + Test problem 1 with process dependent jump rate (sigmoid function)
#
#############################################################################

# Load libraries
library(Rcpp)
library(RcppNumerical)
library(devtools)



h <- 10^(-2) # step size
t_fin <- 100 # time horizon



#############################################################################

# tp 1

sourceCpp(file="TP1_Cpp.cpp")


# Set parameters
sigma <- 1
b <- 2
lambda <- 0.1
X0 <- 0
z0 <- b
eta <- 0.5

# Set seed
set.seed(5)

# Simulate PDifMP
TP1 <- tp1_cpp_t(t_fin,h,X0,z0,sigma,b,lambda,eta)
TP1_path <- TP1$X
TP1_jumps <- TP1$jumps

# Compute grid
TP1_grid <- c()
for(i in 1: (length(TP1_jumps)-1)){
  TP1_grid <- append(TP1_grid, seq(from=TP1_jumps[i],to=TP1_jumps[i+1],by=h))
}
if(grid[length(TP1_jumps)] != t_fin){
  TP1_grid <- append(TP1_grid,t_fin)
}




#############################################################################

# tp 2

sourceCpp(file="TP2_Cpp.cpp")


# Set parameters
sigma <- 1
b <- 20
lambda <- 0.1
eta <- 1
X0 <- c(1,1)
z0 <- b

# Set seed
set.seed(1)

# Simulate PDifMP
TP2 <- tp2_cpp_t(t_fin,h,X0,z0,sigma,b,lambda,eta)
TP2_1_path <- TP2$X[1,]
TP2_jumps <- TP2$jumps

TP2_grid <- seq(0,t_fin,h)




#############################################################################

# tp 3

sourceCpp(file="TP3_Cpp.cpp")

# Set parameters
sigma <- 2
b <- 2
lambda <- 0.1
X0 <- 0
z0 <- 1

# Set seed
set.seed(1)

# Simulate PDifMP
TP3 <- tp3_cpp_t(t_fin,h,X0,z0,sigma,b,lambda)
TP3_path <- TP3$X
TP3_jumps <- TP3$jumps

# Compute grid
TP3_grid <- c()
for(i in 1: (length(TP3_jumps)-1)){
  TP3_grid <- append(TP3_grid, seq(from=TP3_jumps[i],to=TP3_jumps[i+1],by=h))
}
if(TP3_grid[length(TP3_grid)] != t_fin){
  TP3_grid <- append(TP3_grid,t_fin)
}




#############################################################################

# tp 4

sourceCpp(file="TP4_Cpp.cpp")

# Set parameters
sigma <- 1
b <- 0.9
lambda <- 0.1
eta <- 2 
X0 <- c(1,1)
z0 <- b

# Set seed
set.seed(11)

# Simulate PDifMP
TP4 <- tp4_cpp_t(t_fin,h,X0,z0,sigma,b,lambda,eta)
TP4_1_path <- TP4$X[1,]
TP4_jumps <- TP4$jumps

TP4_grid <- seq(0,t_fin,h)








#############################################################################

# tp 5

sourceCpp(file="TP5_Cpp.cpp")


# Set parameters
sigma <- 1
b <- 2
lambda <- 0.1
X0 <- 0
z0 <- b
eta <- 0.5

# Set seed
set.seed(3)

# Simulate PDifMP
TP5 <- tp5_cpp_t(t_fin,h,X0,z0,sigma,b,lambda,eta)
TP5_path <- TP5$X
TP5_jumps <- TP5$jumps
TP5_gridPoints <- TP5$gridPoints

# Compute grid
TP5_grid <- c()
for(i in 1: (length(TP5_gridPoints)-1)){
  TP5_grid <- append(TP5_grid, seq(from=TP5_gridPoints[i],to=TP5_gridPoints[i+1],by=h))
}
if(grid[length(TP5_gridPoints)] != t_fin){
  TP5_grid <- append(TP5_grid,t_fin)
}


###########################################################################

# plot paths

par(mfrow=c(5,1),mai=c(0.25,0.23,0.1,0.1)) 



#-----------------------------
# Path of TP1: Ornstein Uhlenbeck 
#-----------------------------

plot(TP1_grid,TP1_path,type="l",col="blue",ylim=c(-5,5),xlab="",ylab="")

abline(v=TP1_jumps,lty=2)
abline(h=0)

mtext(side=2,expression(paste("TP 1 (OU)")),line=2.0,cex=0.9)
mtext(expression(X[t]),side=3,line=-2,at=-2,cex=0.6)


#-----------------------------
# Path of TP2: Weakly damped stochastic harmonic oscillator
#-----------------------------

plot(TP2_grid,TP2_1_path,type="l",col="blue",ylim=c(-1,1),xlab="",ylab="")


abline(v=TP2_jumps,lty=2)
abline(h=0)

mtext(side=2,expression(paste("TP 2 (WDSHO)")),line=2.0,cex=0.9)
mtext(expression(X[t]^{(1)}),side=3,line=-2,at=-2,cex=0.6)


#-----------------------------
# Path of TP3: Wiener process with drift
#-----------------------------

plot(TP3_grid,TP3_path,type="l",col="blue",ylim=c(-50,30),xlab="",ylab="")

abline(v=TP3_jumps,lty=2)
abline(h=0)

mtext(side=2,expression(paste("TP 3 (WPWD)")),line=2.0,cex=0.9)
mtext(expression(X[t]),side=3,line=-2,at=-2,cex=0.6)


#-----------------------------
# Path of TP4: Switched stochastic harmonic oscillator
#-----------------------------

plot(TP4_grid,TP4_1_path,type="l",col="blue",xlab="",ylab="")

abline(v=TP4_jumps,lty=2)
abline(h=0)

mtext(side=2,expression(paste("TP 4 (Switched SHO)")),line=2.0,cex=0.9)
mtext(expression(X[t]^{(1)}),side=3,line=-2,at=-2,cex=0.6)




#-----------------------------
# Path of TP5: Ornstein Uhlenbeck with sigmoid jump rate function
#-----------------------------

plot(TP5_grid,TP5_path,type="l",col="blue",ylim=c(-5,5),xlab="",ylab="")

abline(v=TP5_jumps,lty=2)
abline(h=0)

mtext(side=2,expression(paste("TP 5 (OU)")),line=2.0,cex=0.9)
mtext(expression(X[t]),side=3,line=-2,at=-2,cex=0.6)





