// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace std;






//------------------------------------------------------------------------------
// Exponential matrix (Simple Oscillator)
//------------------------------------------------------------------------------

// Input:
// t      time point
// z      frequency of oscillator 

// Output: Exponential matrix

// [[Rcpp::export]]
NumericMatrix exp_mat_simple_Cpp(double t, double z){
  NumericMatrix ret(2,2);

  ret(0,_) = NumericVector::create(cos(z*t), sin(z*t)/z);
  ret(1,_) = NumericVector::create(-z*sin(z*t), cos(z*t));
  return ret;
}


//------------------------------------------------------------------------------
// Covariance matrix (Simple Oscillator)
//------------------------------------------------------------------------------

// Input:
// t      time point
// sig    diffusion coefficient
// z      frequency

// Output: Covariance matrix

// [[Rcpp::export]]
NumericMatrix cov_mat_simple_Cpp(double t, double sig, double z){
  NumericMatrix ret(2,2);

  ret(0,_) = NumericVector::create((2*z*t-sin(2*z*t))/(2*pow(z,3)), (pow(sin(z*t),2))/(pow(z,2)));
  ret(1,_) = NumericVector::create((pow(sin(z*t),2))/(pow(z,2)), (2*z*t+sin(2*z*t))/(2*z));
  return (pow(sig,2)/2)*ret;
}


//------------------------------------------------------------------------------
// Exponential matrix (Damped Oscillator)
//------------------------------------------------------------------------------

// Input:
// t      time point
// a      damping parameter
// b      frequency of oscillator  

// Output: Exponential matrix

// [[Rcpp::export]]
NumericMatrix exp_mat_Cpp(double t, double a, double b){
  NumericMatrix ret(2,2);
  double d = sqrt(pow(b,2)-pow(a,2));
  double cdt = cos(d*t);
  double sdt = sin(d*t);
  ret(0,_) = NumericVector::create(cdt+(a/d)*sdt,sdt/d);
  ret(1,_) = NumericVector::create((-pow(b,2)/d)*sdt,cdt-(a/d)*sdt);
  return exp(-a*t)*ret;
}


//------------------------------------------------------------------------------
// Covariance matrix (Damped Oscillator)
//------------------------------------------------------------------------------

// Input:
// t      time point
// sig    diffusion coefficient
// a      damping parameter
// b      frequency of oscillator  

// Output: Covariance matrix

// [[Rcpp::export]]
NumericMatrix cov_mat_Cpp(double t, double sig, double a, double b){
  NumericMatrix ret(2,2);
  double d1 = sqrt(pow(b,2)-pow(a,2));
  double d2 = pow(b,2)-pow(a,2);
  double em2at = exp(-2*a*t);
  double c2dt = cos(2*d1*t);
  double s2dt = sin(2*d1*t);
  double sdt2 = pow(sin(d1*t),2);
  ret(0,_) = NumericVector::create(((pow(sig,2))/(-4*a*pow(b,2)*d2))*(em2at*pow(b,2)+pow(a,2)-pow(b,2)-em2at*pow(a,2)*c2dt+em2at*a*d1*s2dt),((pow(sig,2))/(2*d2))*em2at*sdt2);
  ret(1,_) = NumericVector::create(((pow(sig,2))/(2*d2))*em2at*sdt2,((pow(sig,2))/(-4*a*d2))*(em2at*pow(b,2)+pow(a,2)-pow(b,2)-em2at*pow(a,2)*c2dt-em2at*a*d1*s2dt));
  return ret;
}



//------------------------------------------------------------------------------
// Cholesky decomposition
//------------------------------------------------------------------------------

// Input:
// mat    NumericMatrix
// z      frequency of oscillator 
// t      time point


// Output: Matrix of Cholesky decomposition

// [[Rcpp::export]]
NumericMatrix cholesky_decomp(Rcpp::NumericMatrix mat, double z, double t) {
  arma::mat A = Rcpp::as<arma::mat>(mat);     
  arma::mat R;
  try {
    R = arma::chol(A);  // Cholesky decomposition
  } catch (const std::exception& e) {
    Rcpp::Rcout << "Error in the cholesky decomposition: " << e.what() << std::endl;
    Rcpp::Rcout << "A: " << A << std::endl;
    Rcpp::Rcout << "z: " << z << std::endl;
    Rcpp::Rcout << "t: " << t << std::endl;
    return Rcpp::wrap(NumericMatrix());  // Return empty matrix, if an error occurs
  }
  
  return Rcpp::wrap(R);                      
}

//------------------------------------------------------------------------------
// Matrix Vector multiplication
//------------------------------------------------------------------------------

// Input:
// M      Numericmatrix
// v      NumericVector

// Output: product of matrix and vector

// [[Rcpp::export]]
NumericVector mv_zlt_(NumericMatrix M, NumericVector v) {
  arma::mat A(M.begin(), M.nrow(), M.ncol(), false);  // 'false' = no copy
  arma::colvec x(v.begin(), v.size(), false);
  arma::colvec result = A * x;
  return NumericVector(result.begin(), result.end());
}


//------------------------------------------------------------------------------
// Transpose matrix
//------------------------------------------------------------------------------

// Input:
// A      NumericMatrix

// Output: Transposed matrix


// [[Rcpp::export]]
NumericMatrix transpose(NumericMatrix A) {

  arma::mat mat_A = Rcpp::as<arma::mat>(A);
  arma::mat transposed_A = mat_A.t();
  
  return Rcpp::wrap(transposed_A);
}



//------------------------------------------------------------------------------
// Solution of the Switched stochastic harmonic oscillator (SwitchedSHO)
//------------------------------------------------------------------------------

// Input:
// startv   initial value
// t_fin    end point of one interval
// t0       starting point of one interval
// h        step size
// sig      diffusion coefficient
// z        damping parameter 
// a        frequency of WDSHO

// Output: path of the Switched SHO PDifMP from t_0 to t_fin

// [[Rcpp::export]]
NumericMatrix sde_SHO_cpp_t(NumericVector startv,double t_fin,double t0,double h,double sig, double z, double a){
  
  double eps = 1e-8;                      // add eps to avoid eroors caused by floating-point rounding
  int N = floor((t_fin - t0) / h + eps);  // number of steps of size h
  
  NumericVector newv = startv;            // set first value to startv
  NumericVector randvec(2);
  
  int num_cols = N + 1;
  
  
  NumericMatrix sol(2,num_cols);          // create a Matrix object for the path of size N + 1
  NumericMatrix randarr(2,num_cols-1);    // generate 2N standard normal distributed random values 
  randarr(0,_) = rnorm(num_cols-1);
  randarr(1,_) = rnorm(num_cols-1);
  
  
  sol(_, 0) = startv;
  
  // path of the WDSHO / Simple SHO
  
  if(N > 0){
    for(int i=0;i<N;i++)
    {
      NumericMatrix dm;
      NumericMatrix cm;
      
      if(z==0){
        dm  = exp_mat_simple_Cpp(h,a);
        cm = transpose(cholesky_decomp(cov_mat_simple_Cpp(h,sig,a),z,h));
      }
      else{
        dm  = exp_mat_Cpp(h,z,a);
        cm = transpose(cholesky_decomp(cov_mat_Cpp(h,sig,z,a),z,h));
      }
      
      
      
      randvec = randarr(_,i);
      newv = mv_zlt_(dm,newv) + mv_zlt_(cm,randvec);
      sol(_,(i+1)) = newv;
    }
  }
  
  return sol(_, Range(1,N));      // return path without the initial value
}


//------------------------------------------------------------------------------
// Append function for two matrices
//------------------------------------------------------------------------------

// Input:
// A      NumericMatrix
// B      NumericMatrix


// Output: matrix containing A and B


// [[Rcpp::export]]
NumericMatrix matrixAppend(NumericMatrix A, NumericMatrix B) {
  int nrow_A = A.nrow();
  int ncol_A = A.ncol();
  int ncol_B = B.ncol();
  
  NumericMatrix result(nrow_A, ncol_A + ncol_B);
  std::copy(A.begin(), A.end(), result.begin());
  std::copy(B.begin(), B.end(), result.begin() + ncol_A * nrow_A);
  
  return result;
}


//------------------------------------------------------------------------------
// Transition kernel
//------------------------------------------------------------------------------

// Input:
// z      previous value of z
// b      possible new value for z


// Output: new value for z

// [[Rcpp::export]]
double drawZ_binary_cpp_t(double z, double b){
  if(z == b){
    return 0;
  }
  return b;
}


//------------------------------------------------------------------------------
// Solution of the OU PDifMP
//------------------------------------------------------------------------------

// Input:
// t_fin    time horizon
// h        step size
// startv   initial value of the continuous component
// z0       initial value of the jump component
// sig      diffusion coefficient
// b        parameter of the jump component
// lambda0  constant jump rate
// a        frequency

// Output: path of the Switched SHO PDifMP from 0 to t_fin 
//            (
//              X           vector of continuous componen t
//              zVec        vector of jump component
//              jumps       vector of jump times
//            )


// [[Rcpp::export]]
List tp4_cpp_t(int t_fin,double h,NumericVector startv,double z0,double sig,double b,double lambda0,double a){
  
  // set initial values
  
  double z = z0;
  NumericVector zVec = {z0};
  
  NumericMatrix X(2, 1);
  X(_, 0) = startv;
  
  
  double t0 = 0;
  NumericVector jumps = {t0};
  
  
  NumericVector u;
  NumericVector tVec;
  double t = t0; // enter while loop
  
  
  while(t == t0){                         // jump time should be different from t0
    u = {runif(1)};                       // generate standard uniform distributed random variable
    tVec = t0 - (1.0/lambda0)*log(u[0]);  // generate exponential distributed random variable with rate lambda0
    tVec = round(tVec,-log10(h));         // round jump time to avoid numerical instability
    t = tVec[0];
  }
  
  while(t < t_fin){
    jumps.push_back(t);
    
    // Simulate solution of the Switched SHO until jump time t
    
    NumericMatrix tempX =  sde_SHO_cpp_t(X(_, X.ncol() - 1),t,t0,h,sig,z,a);
    X = matrixAppend(X,tempX);
    
    // new value for z according to transition kernel
    
    z = drawZ_binary_cpp_t(z,b);
    zVec.push_back(z); 
    
    t0 = t;
    
    while(t == t0){                     // jump time should be different from previous jump time t0
      u = {runif(1)};                   // generate standard uniform distributed random variable
      tVec = t - (1/lambda0)*log(u[0]); // generate exponential distributed random variable with rate lambda0
      tVec = round(tVec,-log10(h));     // round jump time to avoid numerical instability
      t = tVec[0];
    }
  };
  
  // Simulate solution of the Switched SHO until time horizon t_fin
  
  X = matrixAppend(X, sde_SHO_cpp_t(X(_, X.ncol() - 1),t_fin,t0,h,sig,z,a));
  
  jumps.push_back(t_fin);
  z = drawZ_binary_cpp_t(z,b);
  zVec.push_back(z); 
  
  return Rcpp::List::create(Rcpp::Named("X") = X,
                            Rcpp::Named("zVec") = zVec,
                            Rcpp::Named("jumps") = jumps);
}




