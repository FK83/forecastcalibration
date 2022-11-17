// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <cmath>
#define _USE_MATH_DEFINES
#include <math.h>
using namespace Rcpp;
using namespace R;
using namespace stats;

// [[Rcpp::export]]
double euclnormC(arma::colvec x){
  double out = sqrt(sum(square(x)));
  return(out);
}

// [[Rcpp::export]]
List copulapitC(arma::colvec y, arma::mat dat){
  double m = dat.n_cols; // size of MC sample
  double c1 = 0; // init count
  arma::colvec c2 = arma::zeros(m);
  for (int i = 1; i < (m+1); i++) {
    c1 += arma::all(dat.col(i-1) <= y);
    for (int j = (i+1); j < (m+1); j++) {
      c2(j-1) += arma::all(dat.col(i-1) <= dat.col(j-1));
      c2(i-1) += arma::all(dat.col(j-1) <= dat.col(i-1));
    }
  }
  double Hy = c1/m;
  arma::colvec KH = c2/(m-1);

  return List::create(Named("Hy") = Hy, Named("KH") = KH);
}

// [[Rcpp::export]]
List esC_sample(arma::colvec y, arma::mat dat){

  double s1 = 0;
  double m = dat.n_cols;
  for (int i = 1; i < (m+1); i++) {
    s1 += euclnormC(dat.col(i-1) - y);
  }

  arma::colvec dist = arma::zeros(m);
  double s2 = 0;
  for (int i = 1; i < (m+1); i++) {
    for (int j = i; j < (m+1); j++) {
      double d = euclnormC(dat.col(i-1) - dat.col(j-1));
      s2 += 2*d;
      dist[i-1] += d;
      dist[j-1] += d;
    }
  }

  double eyx = (s1 / m);
  double exx = (s2 / pow(m, 2));
  arma::colvec avgdist = dist/m;
  return List::create(Named("eyx") = eyx, Named("exx") = exx, Named("avgdist") = avgdist);

}

// [[Rcpp::export]]
List esC_split(arma::colvec y, arma::mat dat1, arma::mat dat2){

  double s1 = 0;
  double j1 = dat1.n_cols;
  for (int i = 1; i < (j1+1); i++) {
    s1 += euclnormC(dat1.col(i-1) - y);
  }

  double j2 = dat2.n_cols;
  arma::colvec dist = arma::zeros(j2);
  for (int i = 1; i < (j1+1); i++) {
    for (int j = 1; j < (j2+1); j++) {
      double d = euclnormC(dat1.col(i-1) - dat2.col(j-1));
      dist[j-1] += d;
    }
  }

  double eyx = (s1/j1);
  arma::colvec avgdist = (dist/j2);
  return List::create(Named("eyx") = eyx, Named("avgdist") = avgdist);

}
