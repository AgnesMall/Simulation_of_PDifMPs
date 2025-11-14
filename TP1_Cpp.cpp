#include <Rcpp.h>
using namespace Rcpp;
using namespace std;


// [[Rcpp::export]]
NumericVector sde_cpp_t(double z,double sig,double X0,double t_fin,double t0,double h, double eta){
  int N = floor((t_fin - t0)/h);
  double lastStepSize = t_fin - t0 - N*h;
  NumericVector X(N+1);
  X[0] = X0;
  
  NumericVector phi = rnorm(N,0,1);
  
  if(N > 0){
    for(int i=0; i<N; i++){
      X[i+1] = exp(-eta*h)*X[i] + z*(1-exp(-eta*h)) + phi[i]*sig*sqrt((1-exp(-2.0*eta*h))/(2.0*eta)); 
    }
  }
  
  if(lastStepSize != 0){
    phi = {rnorm(1)};
    double temp = exp(-eta*h)*X[X.length()-1] + z*(1-exp(-eta*h)) + phi[0]*sig*sqrt((1-exp(-2.0*eta*h))/(2.0*eta));
    X.push_back(temp);
    return X[Range(1, N+1)];
  }
  return X[Range(1, N)];
  
}


// [[Rcpp::export]]
NumericVector vectorAppend_cpp_t(NumericVector x, NumericVector y) {
  NumericVector result(x.size() + y.size());
  std::copy(x.begin(), x.end(), result.begin()); 
  std::copy(y.begin(), y.end(), result.begin() + x.size()); 
  return result;
}




// [[Rcpp::export]]
List tp1_cpp_t(int t_fin,double h,double X0,double z0,double sig,double b,double lambda0,double eta){
  double z = z0;
  NumericVector zVec = {z0};
  
  NumericVector X = {X0};
  int lengthX = X.length();
  
  double t0 = 0;
  NumericVector jumps = {t0};
  NumericVector trueJumps = {t0};
  NumericVector jumpIndices = {0};
  double jumpIndicesCounter = 0;
  
  NumericVector u = {runif(1)};
  double t = t0 - (1.0/lambda0)*log(u[0]);
  
  while(t < t_fin){
    jumps.push_back(t);
    
    NumericVector tempX = sde_cpp_t(z,sig,X[lengthX-1],t,t0,h,eta);
    X = vectorAppend_cpp_t(X,tempX);
    lengthX = X.length();
    
    jumpIndicesCounter += tempX.length();
    
    if(X[lengthX-1] > 0){
      z = -b;
    }
    else{
      z = b;
    }
      
    zVec.push_back(z);
    if(z *zVec[zVec.length()-2] < 0){
      trueJumps.push_back(t);
      jumpIndices.push_back(jumpIndicesCounter);
    }
    
    
    t0 = t;
    u = {runif(1)};
    t = t - (1/lambda0)*log(u[0]); 
  };
  
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
