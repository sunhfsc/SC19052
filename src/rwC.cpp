#include <Rcpp.h>
using namespace Rcpp;

//' @title A random walk Metropolis sampler using Rcpp
//' @description A random walk Metropolis sampler for generating the standard Laplace distribution
//' @param x0 The initial point (double)
//' @param sigma The standard deviation in the normal distribution (double)
//' @param N The length of the chain (int)
//' @return The random numbers and the number of rejections (list)
//' @useDynLib SC19052
//' @import Rcpp
//' @examples
//' \dontrun{
//' N<-1000
//' sigma<-0.05
//' x0<-25
//' rwC(sigma,x0,N)
//' }
//' @export
// [[Rcpp::export]]
List rwC(double sigma,double x0,int N) {
  NumericVector x(N);
  x[0] = x0;
  NumericVector u(N);
  u = runif(N,0,1);
  int k=0;
  
  
  for (int i = 1; i < N; ++i) {
    double y = rnorm(1,x[i-1],sigma)[0];
    double dsk1 = exp(-abs(y))/2;
    double dsk2 = exp(-abs(x[i-1]))/2;
    if(u[i] <= dsk1/dsk2) {
      x[i] = y;
    }else {
      x[i] = x[i-1];
      k = k+1;
    }
  }
  return List::create(
    _["x"] = x,
    _["k"] = k
  );
}
