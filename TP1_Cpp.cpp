#include <Rcpp.h>
using namespace Rcpp;
using namespace std;


//------------------------------------------------------------------------------
// Solution of the Ornstein Uhlenbeck process (OU)
//------------------------------------------------------------------------------

// Input:
// z      mean of OU
// sig    diffusion coefficient
// X0     initial value
// t_fin  end point of one interval
// t0     starting point of one interval
// h      step size
// eta    mean reversion rate

// Output: path of the OU from t_0 to t_fin


// [[Rcpp::export]]
NumericVector sde_cpp_t(double z,double sig,double X0,double t_fin,double t0,double h, double eta){
  int N = floor((t_fin - t0)/h);             // number of steps of size h
  double lastStepSize = t_fin - t0 - N*h;    // step size from last equidistant grid point to interval endpoint  
  
  
  NumericVector X(N+1); // create a Vector object for the path of size N + 1
  X[0] = X0;            // set first value to X0
  
  NumericVector phi = rnorm(N,0,1); // generate N standard normal distributed random values 
  
  
  // path of the OU
  
  if(N > 0){
    for(int i=0; i<N; i++){
      X[i+1] = exp(-eta*h)*X[i] + z*(1-exp(-eta*h)) + phi[i]*sig*sqrt((1-exp(-2.0*eta*h))/(2.0*eta)); 
    }
  }
  
  if(lastStepSize != 0){
    phi = {rnorm(1)};        // generate 1 standard normal distributed random value
    double temp = exp(-eta*h)*X[X.length()-1] + z*(1-exp(-eta*h)) + phi[0]*sig*sqrt((1-exp(-2.0*eta*h))/(2.0*eta));
    X.push_back(temp);
    return X[Range(1, N+1)];
  }
  return X[Range(1, N)]; // return path without the initial value
  
}


//------------------------------------------------------------------------------
// Append function for two vectors
//------------------------------------------------------------------------------

// Input:
// x      NumericVector
// y      NumericVector


// Output: vector containing x and y


// [[Rcpp::export]]
NumericVector vectorAppend_cpp_t(NumericVector x, NumericVector y) {
  NumericVector result(x.size() + y.size());
  std::copy(x.begin(), x.end(), result.begin()); 
  std::copy(y.begin(), y.end(), result.begin() + x.size()); 
  return result;
}


//------------------------------------------------------------------------------
// Solution of the OU PDifMP
//------------------------------------------------------------------------------

// Input:
// t_fin    time horizon
// h        step size
// X0       initial value of the continuous component
// z0       initial value of the jump component
// sig      diffusion coefficient
// b        parameter of the jump component
// lambda0  constant jump rate
// eta      mean reversion rate

// Output: path of the OU PDifMP from 0 to t_fin 
//            (
//              X           vector of continuous componen t
//              zVec        vector of jump component
//              jumps       vector of jump times
//              trueJumps   vector of jumps without no-move jumps
//              jumpIndices vector of indeces of vector X, where a true jump occurs
//            )

// [[Rcpp::export]]
List tp1_cpp_t(int t_fin,double h,double X0,double z0,double sig,double b,double lambda0,double eta){
  
  // set initial values
  
  double z = z0;
  NumericVector zVec = {z0};
  
  NumericVector X = {X0};
  int lengthX = X.length();
  
  double t0 = 0;
  NumericVector jumps = {t0};
  NumericVector trueJumps = {t0};  // jumps without no-move jumps
  NumericVector jumpIndices = {0}; // index of vector X, where a true jump occurs
  double jumpIndicesCounter = 0;
  
  NumericVector u = {runif(1)};             // generate standard uniform distributed random variable
  double t = t0 - (1.0/lambda0)*log(u[0]);  // generate exponential distributed random variable with rate lambda0
  
  while(t < t_fin){
    jumps.push_back(t);
    
    // Simulate solution of OU until jump time t
    
    NumericVector tempX = sde_cpp_t(z,sig,X[lengthX-1],t,t0,h,eta);
    X = vectorAppend_cpp_t(X,tempX);
    lengthX = X.length();
    
    jumpIndicesCounter += tempX.length();
    
    
    // new value for z according to transition kernel
    
    if(X[lengthX-1] > 0){
      z = -b;
    }
    else{
      z = b;
    }
    
    // check if z did change, i.e. a move happend
      
    zVec.push_back(z);
    if(z *zVec[zVec.length()-2] < 0){
      trueJumps.push_back(t);
      jumpIndices.push_back(jumpIndicesCounter);
    }
    
    // draw new jump time
    
    t0 = t;
    u = {runif(1)};
    t = t - (1/lambda0)*log(u[0]); 
  };
  
  // Simulate solution of OU until time horizon t_fin
  
  NumericVector tempX = sde_cpp_t(z,sig, X[lengthX-1], t_fin, t0, h,eta);
  X = vectorAppend_cpp_t(X,tempX);
  lengthX = X.length();
  
  jumpIndicesCounter += tempX.length();
  trueJumps.push_back(t_fin);
  jumpIndices.push_back(jumpIndicesCounter);
  
  jumps.push_back(t_fin);
  
  if(X[lengthX-1] > 0){
    z = -b;
  }
  else{
    z = b;
  }
  
  zVec.push_back(z);
  
  return Rcpp::List::create(Rcpp::Named("X") = X,
                            Rcpp::Named("zVec") = zVec,
                            Rcpp::Named("jumps") = jumps,
                            Rcpp::Named("trueJumps") = trueJumps,
                            Rcpp::Named("jumpIndices") = jumpIndices);
}
