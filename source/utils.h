#ifndef UTILS_H_
#define UTILS_H_

#include <cmath>
#include <iomanip>
#include <sstream>
#include "common.h"

inline double p2e(double p, double E0){ // convert momentum to energy
  return sqrt(p * p * gC * gC + E0 * E0) - E0; 
}

inline double e2p(double E, double E0){ // convert energy to momentum
  return sqrt(E * (E + 2 * E0)) / gC; 
}

inline double e2logMu(double E, double E0, double L){ // NOTICE: ignore alpha, m, BE terms, only for interpolation
  return log(pow(e2p(E, E0), 2) * pow(L, 3));
}

inline double dlogE_dp(double logE, double E0){
  double E = exp(logE);
  return e2p(E, gE0) * gC * gC / (E * (E + E0));
}

inline double a2Y(double alpha){
  double y = sin(alpha);
  double T0 = 1.3802;
  double T1 = 0.7405;

  return 2 * T0 * (1 - y) + (T0 - T1) * (y * log(y) + 2 * y - 2 * sqrt(y));
}

inline double a2T(double alpha){
  double y = sin(alpha);
  double T0 = 1.3802;
  double T1 = 0.7405;

  return T0 - (T0 - T1) * (y + sqrt(y)) / 2.0;
}

// for reading file to set date / hour number digit alignment
inline std::string intToStringWithPadding(int value, int totalWidth) {
    std::ostringstream oss;
    oss << std::setw(totalWidth) << std::setfill('0') << value;
    return oss.str();
}

// for common structured mesh cell interpolation
inline double interp3D(const Xtensor3d& raw, const Loc& loc){
  int i0,j0,k0;
  double wi,wj,wk; 

  i0 = loc.i0;
  j0 = loc.j0;
  k0 = loc.k0;
  wi = loc.wi;
  wj = loc.wj; 
  wk = loc.wk; 

  return raw(i0,j0,k0)*wi*wj*wk + raw(i0+1,j0,k0)*(1-wi)*wj*wk + raw(i0+1,j0+1,k0)*(1-wi)*(1-wj)*wk + raw(i0,j0+1,k0)*wi*(1-wj)*wk +
        raw(i0,j0,k0+1)*wi*wj*(1-wk) + raw(i0+1,j0,k0+1)*(1-wi)*wj*(1-wk) + raw(i0+1,j0+1,k0+1)*(1-wi)*(1-wj)*(1-wk) + raw(i0,j0+1,k0+1)*wi*(1-wj)*(1-wk);
}


// for common weight calculation including range judge
inline void calWeight(std::size_t &i, double &w, std::size_t n, double pos){
  if (i >= 0 && i < n) {
    w = 1.0 - (pos - i); 
  } 
  else if (i < 0) {
    i = 0;
    w = 1.0;
  }
  else if (i >= n) {
    i = n - 1; 
    w = 0.0; 
  }
}


#endif
