/*
 * File:        Mesh.cc
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 *              Peng Peng <pp140594@mail.ustc.edu.cn>
 * Date:        05/12/2024 
 * 
 * Copyright (c) Xin Tao 
 *
 */

#include "Mesh.h"

Mesh::Mesh(const Parameters& p): x_(p.nalpha0()), y_(p.nE()), z_(p.nL()) {

  nx_ = p.nalpha0(); 
  ny_ = p.nE(); 
  nz_ = p.nL(); 

  dt_ = p.dt(); 

  xO_ = p.alpha0_min(); 
  yO_ = p.logEmin();
  zO_ = p.Lmin();

  dx_ = (p.alpha0_max() - p.alpha0_min()) / p.nalpha0(); 
  dy_ = (p.logEmax() - p.logEmin()) / p.nE(); 
  dz_ = (p.Lmax() - p.Lmin()) / p.nL();

  x_(0) = xO() + dx()/2.0; 
  y_(0) = yO() + dy()/2.0; 
  z_(0) = zO() + dz()/2.0;

  for (std::size_t i=1; i<nx(); ++i) x_(i) = x_(0) + i*dx(); 
  for (std::size_t j=1; j<ny(); ++j) y_(j) = y_(0) + j*dy(); 
  for (std::size_t k=1; k<nz(); ++k) z_(k) = z_(0) + k*dz(); 

  nbr_inds.resize({nx(), ny(), nz(), nnbrs()});
  faces.resize({nx(), ny(), nz(), nnbrs()});

  build_connectivity(); 

  // The reverse inbr number.
  // For example, if the current cell is K, its 0th neighbor is L.
  // Then, for cell L, K is its 2th neighbor.
  rinbr_(0) = 2; 
  rinbr_(1) = 3; 
  rinbr_(2) = 0; 
  rinbr_(3) = 1; 
  rinbr_(4) = 5; 
  rinbr_(5) = 4; 

}

void Mesh::build_connectivity() {

  int inbr; 

  Point A, B, C, D; 
  Vector3 dr; 

  for (std::size_t i=0; i<nx(); ++i) {
    for (std::size_t j=0; j<ny(); ++j) {
      for (std::size_t k=0; k<nz(); ++k) {
        // nbr 0
        inbr = 0; 
        nbr_inds(i,j,k,inbr).i = i-1;
        nbr_inds(i,j,k,inbr).j = j; 
        nbr_inds(i,j,k,inbr).k = k; 

        A = {xO() + i*dx(), yO() + j*dy(), zO() + k*dz()}; 
        dr = {0, 0, dz()};  

        B = A + dr; 

        dr = {0, dy(), 0}; 
        C = B + dr; 

        dr = {0, 0, -dz()}; 
        D = C + dr; 

        faces(i,j,k,inbr).set_vs_dir({A,B,C,D}, XNEG); 

        // nbr 1
        inbr = 1;
        nbr_inds(i,j,k,inbr).i = i;
        nbr_inds(i,j,k,inbr).j = j+1;
        nbr_inds(i,j,k,inbr).k = k;

        A = D;
        B = C;

        dr = {dx(), 0, 0};
        C = B + dr; 

        dr = {dx(), 0, 0};
        D = A + dr; 

        faces(i,j,k,inbr).set_vs_dir({A,B,C,D}, YPOS); 

        // nbr 2
        inbr = 2; 
        nbr_inds(i,j,k,inbr).i = i+1;
        nbr_inds(i,j,k,inbr).j = j;
        nbr_inds(i,j,k,inbr).k = k;

        A = D;
        B = C;

        dr = {0, -dy(), 0};
        C = B + dr; 

        D = A + dr; 

        faces(i,j,k,inbr).set_vs_dir({A,B,C,D}, XPOS); 

        // nbr 3
        inbr = 3;
        nbr_inds(i,j,k,inbr).i = i;
        nbr_inds(i,j,k,inbr).j = j-1;
        nbr_inds(i,j,k,inbr).k = k;

        A = D;
        B = C;

        dr = {-dx(), 0, 0};
        C = B + dr; 
        D = A + dr; 

        faces(i,j,k,inbr).set_vs_dir({A,B,C,D}, YNEG); 

        // nbr 4
        inbr = 4;
        nbr_inds(i,j,k,inbr).i = i;
        nbr_inds(i,j,k,inbr).j = j;
        nbr_inds(i,j,k,inbr).k = k-1;

        A = {xO() + i*dx(), yO() + j*dy(), zO() + k*dz()};

        dr = {0, dy(), 0}; 
        B = A + dr; 

        dr = {dx(), 0, 0};
        C = B + dr; 

        dr = {0, -dy(), 0}; 
        D = C + dr; 

        faces(i,j,k,inbr).set_vs_dir({A,B,C,D}, ZNEG); 

        // nbr 5
        inbr = 5;
        nbr_inds(i,j,k,inbr).i = i;
        nbr_inds(i,j,k,inbr).j = j;
        nbr_inds(i,j,k,inbr).k = k+1;

        A = {xO() + (i+1)*dx(), yO() + j*dy(), zO() + (k+1)*dz()};

        dr = {0, dy(), 0};
        B = A + dr; 

        dr = {-dx(), 0, 0};
        C = B + dr; 

        dr = {0, -dy(), 0};
        D = C + dr; 

        faces(i,j,k,inbr).set_vs_dir({A,B,C,D}, ZPOS); 
      }
    }
  }
}


