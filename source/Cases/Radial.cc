/*
 * File:        Radial.cc
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 *              Peng Peng <pp140594@mail.ustc.edu.cn>
 * Date:        12/02/2024 
 * 
 * Copyright (c) Xin Tao 
 *
 */

#include "Radial.h"


Radial::Radial(const Parameters& paras_in, const Mesh& mesh_in):Equation(mesh_in),paras(paras_in), m(mesh_in){ // initialize G, D, and etc.
  init();
  constructD(paras);
}

double Radial::init_f(const Ind& ind) const { 

  return calculate_init_f(m.z(ind.k)); 
} 

void Radial::init(){

  double z;

  for (std::size_t i = 0; i < m.nx(); i++){
    for (std::size_t j = 0; j < m.ny(); j++){
      for (std::size_t k = 0; k < m.nz(); k++){
        z = m.z(k);

        G_(i,j,k) = calculate_G(z);
        tau_(i,j,k) = std::numeric_limits<double>::max();
      }
    }
  }
}


void Radial::update(double t) {}


void Radial::apply_bcs(Xtensor3d* vertex_fp) const { 
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
      vertex_f(i, 0, k) = vertex_f(i, 1, k); 
      vertex_f(i, m.ny(), k) = vertex_f(i, m.ny()-1, k); 
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


void Radial::constructD(const Parameters& par){

  double alpha0, logE, p, L, papL, pppL;

  for(std::size_t i = 0; i < m.nx(); i++){
      alpha0 = m.x(i);
      for(std::size_t j = 0; j < m.ny(); j++){
          logE = m.y(j);
          p = e2p(exp(logE), gE0);
          for(std::size_t k = 0; k < m.nz(); k++){
            L = m.z(k);

            papL = -a2Y(alpha0) * tan(alpha0) / (4 * a2T(alpha0) * L);
            pppL = -p * (6 * a2T(alpha0) - a2Y(alpha0)) / (L * 4 * a2T(alpha0));

            Dzz_(i,j,k) = pow(L, 10) / pow(10, 9);

            Dxz_(i,j,k) = Dzz_(i,j,k) * papL;
            Dyz_(i,j,k) = Dzz_(i,j,k) * pppL * dlogE_dp(logE, gE0);
            Dxx_(i,j,k) = Dzz_(i,j,k) * pow(papL, 2);
            Dyy_(i,j,k) = Dzz_(i,j,k) * pow(pppL, 2) * pow(dlogE_dp(logE, gE0),2);
            Dxy_(i,j,k) = Dzz_(i,j,k) * papL * pppL * dlogE_dp(logE, gE0); 
          }
      }
  }
}







