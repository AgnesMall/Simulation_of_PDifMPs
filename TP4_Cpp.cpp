// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace std;


// Matrix vector multiplications
// [[Rcpp::export]]
NumericVector mv_zlt_(NumericMatrix M, NumericVector v) {
  arma::mat A(M.begin(), M.nrow(), M.ncol(), false);  // 'false' = no copy
  arma::colvec x(v.begin(), v.size(), false);
  arma::colvec result = A * x;
  return NumericVector(result.begin(), result.end());
}



// [[Rcpp::export]]
NumericMatrix exp_mat_simple_Cpp(double t, double z){
  NumericMatrix ret(2,2);

  ret(0,_) = NumericVector::create(cos(z*t), sin(z*t)/z);
  ret(1,_) = NumericVector::create(-z*sin(z*t), cos(z*t));
  return ret;
}


// [[Rcpp::export]]
NumericMatrix cov_mat_simple_Cpp(double t, double sig, double z){
  NumericMatrix ret(2,2);

  ret(0,_) = NumericVector::create((2*z*t-sin(2*z*t))/(2*pow(z,3)), (pow(sin(z*t),2))/(pow(z,2)));
  ret(1,_) = NumericVector::create((pow(sin(z*t),2))/(pow(z,2)), (2*z*t+sin(2*z*t))/(2*z));
  return (pow(sig,2)/2)*ret;
}




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




// Chelosky decomposition
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

// Transpose matrix
// [[Rcpp::export]]
NumericMatrix transpose(NumericMatrix A) {

  arma::mat mat_A = Rcpp::as<arma::mat>(A);
  arma::mat transposed_A = mat_A.t();
  
  return Rcpp::wrap(transposed_A);
}



// [[Rcpp::export]]
NumericMatrix sde_SHO_cpp_t(NumericVector startv,double t_fin,double t0,double h,double sig, double z, double a){
  
  double eps = 1e-8;
  int N = floor((t_fin - t0) / h + eps);
  
  NumericVector newv = startv;
  NumericVector randvec(2);
  
  int num_cols = N + 1;
  
  
  NumericMatrix sol(2,num_cols);
  NumericMatrix randarr(2,num_cols-1);
  randarr(0,_) = rnorm(num_cols-1);
  randarr(1,_) = rnorm(num_cols-1);
  
  
  sol(_, 0) = startv;
  
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
  
  return sol(_, Range(1,N));
}



// Append matrix to other matrix
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




// [[Rcpp::export]]
double drawZ_binary_cpp_t(double z, double b){
  if(z == b){
    return 0;
  }
  return b;
}




// [[Rcpp::export]]
List tp4_cpp_t(int t_fin,double h,NumericVector startv,double z0,double sig,double b,double lambda0,double a){
  double z = z0;
  NumericVector zVec = {z0};
  
  NumericMatrix X(2, 1);
  X(_, 0) = startv;
  
  
  double t0 = 0;
  NumericVector jumps = {t0};
  
  
  NumericVector u;
  NumericVector tVec;
  double t = t0; // enter while loop
  
  
  while(t == t0){
    u = {runif(1)};
    tVec = t0 - (1.0/lambda0)*log(u[0]);
    tVec = round(tVec,-log10(h));
    t = tVec[0];
  }
  
  while(t < t_fin){
    jumps.push_back(t);
    
    NumericMatrix tempX =  sde_SHO_cpp_t(X(_, X.ncol() - 1),t,t0,h,sig,z,a);
    X = matrixAppend(X,tempX);
    z = drawZ_binary_cpp_t(z,b);
    zVec.push_back(z); 
    
    
    t0 = t;
    
    while(t == t0){
      u = {runif(1)};
      tVec = t - (1/lambda0)*log(u[0]); 
      tVec = round(tVec,-log10(h));
      t = tVec[0];
    }
  };
  
  X = matrixAppend(X, sde_SHO_cpp_t(X(_, X.ncol() - 1),t_fin,t0,h,sig,z,a));
  
  jumps.push_back(t_fin);
  z = drawZ_binary_cpp_t(z,b);
  zVec.push_back(z); 
  
  return Rcpp::List::create(Rcpp::Named("X") = X,
                            Rcpp::Named("zVec") = zVec,
                            Rcpp::Named("jumps") = jumps);
}




