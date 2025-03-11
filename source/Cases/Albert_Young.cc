/*
 * File:        Albert_Young.cc
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 * Date:        11/27/2024 
 * 
 * Copyright (c) Xin Tao 
 *
 */

#include "Albert_Young.h"


Albert_Young::Albert_Young(const Parameters& paras_in, const Mesh& mesh_in):Equation(mesh_in),paras(paras_in), m(mesh_in), io(paras_in){ // initialize G, D, and etc.
  init();
  constructD(paras, io);
}


double Albert_Young::alpha_lc(double L) const{ // calculate alphaLC from L

  return asin(pow(pow(L,5)*(4*L-3), -0.25)); 
}

double Albert_Young::init_f(const Ind& ind) const { 
  double x, y, z;
  x = m.x(ind.i);
  y = m.y(ind.j);
  z = m.z(ind.k); 

  return calculate_init_f(x, y, z);
}

double Albert_Young::bounce_period(double a0, double logE, double L) const{

  return 4*L * gRE * ((gE0 + exp(logE)) / (gC * gC)) / e2p(exp(logE), gE0) * a2T(a0) / (3e8 * 3600 * 24); 
}


void Albert_Young::init(){

  double x, y, z;

  for (std::size_t i = 0; i < m.nx(); i++){
    x = m.x(i);
    for (std::size_t j = 0; j < m.ny(); j++){
      y = m.y(j);
      for (std::size_t k = 0; k < m.nz(); k++){
        z = m.z(k);

        G_(i,j,k) = calculate_G(x, y, z);

        if (x < alpha_lc(z)) {
          tau_(i,j,k) = bounce_period(x, y, z) / 4.0; 
        }
        else {
          tau_(i,j,k) = std::numeric_limits<double>::max();
        }
      }
    }
  }
}


void Albert_Young::update(double t) {}


void Albert_Young::apply_bcs(Xtensor3d* vertex_fp) const { 
  Xtensor3d& vertex_f = *vertex_fp;  

  double x, y, z;

  // i == 0 and m.nx() boundary condition case
  for (std::size_t j = 0; j<m.ny()+1; ++j){
    for (std::size_t k = 0; k<m.nz()+1; ++k){
      vertex_f(0, j, k) = vertex_f(1, j, k);
      vertex_f(m.nx(), j, k) = vertex_f(m.nx()-1, j, k);
    }
  }

  // j == 0  and j == m.ny() boundary
  for (std::size_t i = 0; i<m.nx()+1; ++i){
    x = m.xO() + i*m.dx(); 
    for (std::size_t k = 0; k<m.nz()+1; ++k){
      z = m.zO() + k*m.dz();
      vertex_f(i, 0, k) = ymin(x, z); 
      vertex_f(i, m.ny(), k) = ymax(x, z); 
    }
  }

  for (std::size_t i = 0; i<m.nx()+1; ++i) {
    for (std::size_t j = 0; j<m.ny()+1; ++j){
      x = m.xO() + i*m.dx(); 
      y = m.yO() + j*m.dy();
      vertex_f(i, j, 0) = zmin(x, y); 
      vertex_f(i, j, m.nz()) = zmax(x, y); 
    }
  }
}


void Albert_Young::locate(double alpha0, double logE, double L, Loc* locp){
    std::size_t i0, j0, k0;
    double wi, wj, wk;

    double pos_x = (alpha0 - io.xmin_D()) / (io.xmax_D() - io.xmin_D()) * (io.nx_D() - 1);
    double pos_y = (logE - log(io.ymin_D())) / (log(io.ymax_D()) - log(io.ymin_D())) * (io.ny_D() - 1);
    double pos_z = (L - io.zmin_D()) / (io.zmax_D() - io.zmin_D()) * (io.nz_D() - 1);

    i0 = floor(pos_x); 
    j0 = floor(pos_y); 
    k0 = floor(pos_z); 

    calWeight(i0, wi, io.nx_D()-1, pos_x);
    calWeight(j0, wj, io.ny_D()-1, pos_y);
    calWeight(k0, wk, io.nz_D()-1, pos_z);

    locp->i0 = i0;
    locp->j0 = j0; 
    locp->k0 = k0; 
    locp->wi = wi;
    locp->wj = wj; 
    locp->wk = wk; 
}

void Albert_Young::constructD(const Parameters& par, const Albert_Young_IO& io){

  const double denormalize_factor = gME * gME * gC * gC;
  const double second_to_day = 3600 * 24;
  double alpha0, logE, p, L;

  Loc loc; 

  for(std::size_t i = 0; i < m.nx(); i++){
      alpha0 = m.x(i);
      for(std::size_t j = 0; j < m.ny(); j++){
          logE = m.y(j);
          p = e2p(exp(logE), gE0);
          for(std::size_t k = 0; k < m.nz(); k++){
            L = m.z(k);

            locate(alpha0, logE, L, &loc);  

            Dxx_(i,j,k) = interp3D(io.Dxx_raw, loc) * denormalize_factor * second_to_day / (p*p); 
            Dxy_(i,j,k) = interp3D(io.Dxy_raw, loc) * denormalize_factor * second_to_day * dlogE_dp(logE, gE0) / p;
            Dyy_(i,j,k) = interp3D(io.Dyy_raw, loc) * denormalize_factor * second_to_day * pow(dlogE_dp(logE, gE0),2);

            Dxz_(i,j,k) = 0.0;
            Dyz_(i,j,k) = 0.0;
            Dzz_(i,j,k) = 0.0;
          }
      }
  }
}







