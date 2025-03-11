/*
 * File:        Equation.h
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 *              Peng Peng <pp140594@mail.ustc.edu.cn>
 * Date:        11/24/2024 
 * 
 * Copyright (c) Xin Tao 
 *
 */

#ifndef EQUATION_H_
#define EQUATION_H_

#include "common.h"
#include "utils.h"
#include "Mesh.h"
#include "Parameters.h"

//
// Definitions of the Equation are given here
// including D, tau, and boundary conditions.
// To handle different cases, we Equation as an abstract class.
//

class Equation{
  public:
    Equation(const Mesh& m) {
      std::size_t nx = m.nx();
      std::size_t ny = m.ny();
      std::size_t nz = m.nz();

      G_.resize({nx,ny,nz});

      Dxx_.resize({nx,ny,nz});
      Dyy_.resize({nx,ny,nz});
      Dzz_.resize({nx,ny,nz});
      Dxy_.resize({nx,ny,nz});
      Dxz_.resize({nx,ny,nz});
      Dyz_.resize({nx,ny,nz});

      tau_.resize({nx,ny,nz});
    }

    double G(const Ind& ind) const { return G_(ind.i, ind.j, ind.k); } 

    double Dxx(const Ind& ind) const { return Dxx_(ind.i, ind.j, ind.k); }
    double Dyy(const Ind& ind) const { return Dyy_(ind.i, ind.j, ind.k); } 
    double Dzz(const Ind& ind) const { return Dzz_(ind.i, ind.j, ind.k); } 

    double Dxy(const Ind& ind) const { return Dxy_(ind.i, ind.j, ind.k); } 
    double Dxz(const Ind& ind) const { return Dxz_(ind.i, ind.j, ind.k); } 
    double Dyz(const Ind& ind) const { return Dyz_(ind.i, ind.j, ind.k); } 

    double tau(const Ind& ind) const { return tau_(ind.i, ind.j, ind.k); }

    // pure virtual functions to be defined by the CASE of interest. 
    virtual double init_f(const Ind& ind) const = 0;
    virtual void apply_bcs(Xtensor3d* vertex_fp) const = 0; 
    virtual void update(double t) = 0;

  protected:
    Xtensor3d G_; 

    Xtensor3d Dxx_; 
    Xtensor3d Dyy_; 
    Xtensor3d Dzz_;
    Xtensor3d Dxy_;
    Xtensor3d Dxz_;
    Xtensor3d Dyz_;

    Xtensor3d tau_; 
};


#endif /* EQUATION_H */

