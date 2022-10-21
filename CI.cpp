// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace arma;
using namespace std;

// #define C0 1/(13/exp(3)-40/exp(2)+13/datum::e)
#define C0 1

// [[Rcpp::export]]
vec mvnormpdf(const mat& x) {
  return exp(-sum(pow(x, 2), 1) / 2.0) / pow(2 * datum::pi, x.n_cols / 2.0);
}

// [[Rcpp::export]]
vec Kernel(const mat& x, const double& h){
  if (x.n_cols == 1)  // normal Kernel N(0,1)
  {
    return normpdf(x / h);
  }
  else
  {
    return mvnormpdf(x / h);
  }
}

// [[Rcpp::export]]
double EstimateUV(const mat& z, const vec& y, const int& i, const double& h){
  return sum(Kernel(z.each_row() - z.row(i), h) % (y <= y(i))) / (sum(Kernel(z.each_row() - z.row(i), h)) + 0.0);
}

// [[Rcpp::export]]
mat EstimateW(const mat& z, const int& i, const double& h) {
  if (z.n_cols == 1)
  {
    mat temp = mat(1,1);
    temp = sum(z.col(0) <= z(i, 0)) / (z.n_rows + 0.0);
    return temp; // +0.0 to convert int to double; or else int/int gives int; double/int or int/double gives double
  }
  else
  {
    mat w = mat(1,z.n_cols);
    w(0,0) = sum(z.col(0) <= z(i,0)) / (z.n_rows + 0.0);
    for (size_t m = 1; m < z.n_cols; m++)
    {
      mat temp = z.cols(0, m - 1);
      w(0,m) = sum(Kernel(temp.each_row() - z.submat(i,0,i,m-1), h) % (z.col(m) <= z(i,m))) / sum(Kernel(temp.each_row() - z.submat(i,0,i,m-1), h));
    }
    return w;
  }
    
}

// [[Rcpp::export]]
double EstimateRho(const vec& u, const vec& v, const mat& w){
  double temp = 0;
  int n = u.size();
  for(int i=0;i<n;i++){
    for(int j=i+1;j<n;j++){
      temp += (exp(-abs(u(i)-u(j))) + exp(-u(i))+exp(u(i)-1)+exp(-u(j))+exp(u(j)-1)+2*exp(-1)-4) * (exp(-abs(v(i)-v(j))) + exp(-v(i))+exp(v(i)-1)+exp(-v(j))+exp(v(j)-1)+2 * exp(-1) -4) * exp(-sum(abs(w.row(i)-w.row(j))));
    }
  }
  temp *= 2;
  /*for(int i=0;i<n;i++){
    temp += (1 + 2*exp(-u(i))+2*exp(u(i)-1)+2 * exp(-1) -4) * (1 + 2*exp(-v(i))+2*exp(v(i)-1)+2 * exp(-1)-4);
  }*/
  //return temp*C0/pow(n,2);
  return temp * C0 / (n*(n-1));
}

// [[Rcpp::export]]
vec CIS(const mat& x, const vec& y, const mat& z, const double& h){
  int n = x.n_rows;
  int p = x.n_cols;
  
  vec rho = vec(p);
  vec v = vec(n);
  mat w = mat(n,z.n_cols);
  vec u = vec(n);
  
  for(int i=0; i<n; i++){
    v(i) = EstimateUV(z,y,i,h);
    w.row(i) = EstimateW(z,i,h);
  }
  
  for(int k=0;k<p;k++){
    for(int i=0; i<n; i++){
      u(i) = EstimateUV(z,x.col(k),i,h);
    }
    rho(k) = EstimateRho(u,v,w);
  }
  return rho;
}

