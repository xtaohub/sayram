/*
 * File:        common.h
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 *              Peng Peng <pp140594@mail.ustc.edu.cn>
 * Date:        05/12/2024 
 * 
 * Copyright (c) Xin Tao 
 *
 */

#ifndef COMMON_H_
#define COMMON_H_

#define EIGEN_NO_DEBUG

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <limits>
#include "Eigen/Dense"
#include "Eigen/Sparse" 
#include "Eigen/Core"
#include "xtensor/xarray.hpp"
#include "xtensor/xtensor.hpp"
#include "xtensor/xio.hpp"
#include "xtensor-io/xhighfive.hpp"

typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SpMat;
typedef Eigen::Triplet<double> T;
typedef Eigen::Vector3d Vector3;
typedef Eigen::Vector3d Point; 

typedef xt::xarray<double> Xarray1d;
typedef xt::xtensor<double, 2> Xtensor2d; 
typedef xt::xtensor<double, 3> Xtensor3d; 

struct Loc{
  int i0; 
  int j0; 
  int k0;
  double wi; 
  double wj; 
  double wk; 
}; 

// constants
const double gEPS = std::numeric_limits<double>::epsilon();
const double gPI = 3.141592653589793238462;
const double gD2R = gPI / 180.0; // convert degree to radian
const double gC = 1;
const double gE0 = 0.511875; // MeV
const double gME = gE0 / (gC * gC); 
const double gRE = 6371000;
const double gBE = 0.312; // G

#endif
