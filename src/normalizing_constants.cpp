#include <Rcpp.h>
using namespace Rcpp;

//' @title Compute COM-Poisson normalizing constante in log scale
//' @description Computes the normalizing constant in the log scale
//'   usign the LogSumExp trick to avoid numerical issues (See details).
//' @details \code{logspace_add(a, b) = log(exp(a) + exp(b)) = b + log(1
//'   + exp(a - b))}, where \code{b > a}.
//' @param loglambda A vector of logarithm of the \eqn{\lambda}
//'   parameter.
//' @param mu A vector of the \eqn{\mu} parameter.
//' @param nu A vector of dispersion parameters \eqn{\nu}.
//' @return The normalizing constant.
//' @references Wikipedia, LogSumExp. \url{http://rstudio.com}.
//' @author Eduardo Jr <edujrrib@gmail.com>
//' @export
// [[Rcpp::export]]
NumericVector compute_logz(NumericVector loglambda,
                           NumericVector mu,
                           double nu) {
  // Control loop
  int maxiter = 1e4;
  double logepsilon = log(1e-8);
  // Output vector
  int n = loglambda.size();
  NumericVector out(n);
  // Compute logz
  for (int i = 0; i < n; ++i) {
    double logz = 0;
    int index = ceil(mu[i]);
    // Left side
    for (int j = 1; j < index; ++j) {
      double logz_ = j * loglambda[i] - nu * R::lgammafn(j + 1);
      logz = R::logspace_add(logz, logz_);
      // if (logz_ < logepsilon) break;
    }
    for (int j = index; j < maxiter; ++j) {
      double logz_ = j * loglambda[i] - nu * R::lgammafn(j + 1);
      logz = R::logspace_add(logz, logz_);
      if (logz_ < logepsilon) break;
    }
    out[i] = logz;
  }
  return out;
}

//' @title Compute Double Poisson normalizing constante in log scale
//' @description Computes the normalizing constant in the log scale
//'   usign the LogSumExp trick to avoid numerical issues (See details).
//' @details \code{logspace_add(a, b) = log(exp(a) + exp(b)) = b + log(1
//'   + exp(a - b))}, where \code{b > a}.
//' @param mu a vector of the \eqn{\mu} parameter.
//' @param lmu logarithm of \eqn{\mu}.
//' @param phi value of \eqn{\varphi} parameter.
//' @param lphi logarithm of \eqn{\varphi}.
//' @return The normalizing constant.
//' @references Wikipedia, LogSumExp. \url{http://rstudio.com}.
//' @author Eduardo Jr <edujrrib@gmail.com>
//' @export
// [[Rcpp::export]]
NumericVector compute_logk(NumericVector mu,
                           NumericVector lmu,
                           double phi,
                           double lphi) {
  // Control loop
  int maxiter = 1e4;
  double logepsilon = log(1e-8);
  // Output vector
  int n = mu.size();
  NumericVector out(n);
  // Compute logz
  for (int i = 0; i < n; ++i) {
    double logk = -lphi/2 - mu[i]/phi;
    int index = ceil(mu[i]);
    // Left side
    for (int j = 1; j < index; ++j) {
      double logk_ = -0.5 * lphi - mu[i]/phi - j + j * log(j) -
        R::lgammafn(j + 1) + (j/phi) * (1 + lmu[i] - log(j));
      logk = R::logspace_add(logk, logk_);
      // if (logk_ < logepsilon) break;
    }
    // Right side
    for (int j = index; j < maxiter; ++j) {
      double logk_ = -0.5 * lphi - mu[i]/phi - j + j * log(j) -
        R::lgammafn(j + 1) + (j/phi) * (1 + lmu[i] - log(j));
      logk = R::logspace_add(logk, logk_);
      if (logk_ < logepsilon) break;
    }
    out[i] = logk;
  }
  return out;
}
