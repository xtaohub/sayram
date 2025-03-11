/*
 * File:        Radial.h
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 *              Peng Peng <pp140594@mail.ustc.edu.cn>
 * Date:        12/02/2024 
 * 
 * Copyright (c) Xin Tao 
 *
 */

#ifndef RADIAL_H_
#define RADIAL_H_

#include "Equation.h"


class Radial:public Equation{
  public:
    Radial(const Parameters& paras_in, const Mesh& m_in);
    void update(double t);  // where we update G, tau, D, and Bcs.

    double init_f(const Ind& ind) const; 

    void apply_bcs(Xtensor3d* vertex_fp) const; 

    double zmin(double a, double logE) const{
      // test
      return 0.0; 
    }

    double zmax(double a, double logE) const{
      // test
      return 1.0;
    }


  private:
    const Parameters& paras; 
    const Mesh& m;

    // Define your boundary condition functions here
    double calculate_init_f(double L) const{
      return (L - 3) / 5.0;
    }

    double calculate_G(double L){
      return 1.0 / pow(L, 2);
    }

    void init(); 
    void constructD(const Parameters& par);
}; 


#endif /* RADIAL_H_ */

