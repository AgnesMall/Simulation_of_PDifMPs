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
  int N = floor((t_fin - t0)/h);            // number of steps of size h
  double lastStepSize = t_fin - t0 - N*h;   // step size from last equidistant grid point to interval endpoint  
  
  NumericVector X(N+1);                     // create a Vector object for the path of size N + 1
  X[0] = X0;                                // set first value to X0
  
  NumericVector phi = rnorm(N,0,1);         // generate N standard normal distributed random values 
  
  
  // path of the OU
  
  if(N > 0){
    for(int i=0; i<N; i++){
      X[i+1] = exp(-eta*h)*X[i] + z*(1-exp(-eta*h)) + phi[i]*sig*sqrt((1-exp(-2.0*eta*h))/(2.0*eta)); 
    }
  }
  
  if(lastStepSize != 0){
    phi = {rnorm(1)};                       // generate 1 standard normal distributed random value
    double temp = exp(-eta*h)*X[X.length()-1] + z*(1-exp(-eta*h)) + phi[0]*sig*sqrt((1-exp(-2.0*eta*h))/(2.0*eta));
    X.push_back(temp);
    return X[Range(1, N+1)];
  }
  return X[Range(1, N)];                  // return path without the initial value
  
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
// Transition kernel
//------------------------------------------------------------------------------

// Input:
// loc    current value of the continuous component
// b      absolute value of z


// Output: new value of z


// [[Rcpp::export]]
double drawz_cpp_t(double loc, double b){
  if(loc > 0){
    return -b;
  }
  return b;
}



//------------------------------------------------------------------------------
// Jump rate function
//------------------------------------------------------------------------------

// Input:
// x           current value of the continuous component
// lambda0     sclaing of the jump rate function


// Output: jump rate


// [[Rcpp::export]]
double jumpRateFunction(double x, double lambda0){
  return lambda0/(1+exp(-x));
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
List tp5_cpp_t(int t_fin,double h,double X0,double z0,double sig,double b,double lambda0,double eta){
  
  // set initial values
  
  double z = z0;
  NumericVector zVec = {z0};
  
  NumericVector X = {X0};
  int lengthX = X.length();
  
  double t0 = 0;
  NumericVector jumps = {t0};
  NumericVector gridPoints = {t0};
  
  NumericVector u = {runif(1)};             // generate standard uniform distributed random variable
  double t = t0 - (1.0/lambda0)*log(u[0]);  // generate exponential distributed random variable with rate lambda0
  
  while(t < t_fin){
    
    gridPoints.push_back(t);
    
    // Simulate solution of OU until jump time t
    
    NumericVector tempX = sde_cpp_t(z,sig,X[lengthX-1],t,t0,h,eta);
    X = vectorAppend_cpp_t(X,tempX);
    lengthX = X.length();
    
    NumericVector v = {runif(1)};       // generate standard unifrom distributed random variabel
    if(v[0] <= jumpRateFunction(X[lengthX-1],lambda0)/lambda0){ // accept jump with probability jumpRateFunction(x)/lambda0
      jumps.push_back(t);
      
      // new value for z according to transition kernel
      
      z = drawz_cpp_t(X[lengthX-1],b); 
      zVec.push_back(z); 
      
    }
    
    // draw new jump time candidate
    
    t0 = t;
    u = {runif(1)};
    t = t - (1/lambda0)*log(u[0]); 
  };
  
  // Simulate solution of OU until time horizon t_fin
  
  X = vectorAppend_cpp_t(X,sde_cpp_t(z,sig, X[lengthX-1], t_fin, t0, h,eta));
  lengthX = X.length();
  
  jumps.push_back(t_fin);
  gridPoints.push_back(t_fin);
  zVec.push_back(drawz_cpp_t(X[lengthX-1],b));
  
  return Rcpp::List::create(Rcpp::Named("X") = X,
                            Rcpp::Named("zVec") = zVec,
                            Rcpp::Named("jumps") = jumps,
                            Rcpp::Named("gridPoints") = gridPoints);
}
