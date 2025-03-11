/*
 * File:        Albert_Young.h
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 *              Peng Peng 
 * Date:        11/27/2024 
 * 
 * Copyright (c) Xin Tao 
 *
 */

#ifndef ALBERT_YOUNG_H_
#define ALBERT_YOUNG_H_

#include "Equation.h"
#include "Albert_Young_IO.h"
#include "utils.h"


class Albert_Young:public Equation{
  public:
    Albert_Young(const Parameters& paras_in, const Mesh& m_in);

    double init_f(const Ind& ind) const; 
    void apply_bcs(Xtensor3d* vertex_fp) const; 
    void update(double t);  // where we update G, tau, D, and Bcs.
                            
    double zmin(double a, double logE) const{
      // test
      return 0.0; 
    }

    double zmax(double a, double logE) const{
      // test
      return 1.0;
    }

    double ymin(double a0, double L) const{
      return calculate_init_f(a0, paras.logEmin(), L); 
    }

    double ymax(double a0, double L) const{
      return 0.0;
    }

  private:
    const Parameters& paras; 
    const Mesh& m;
    Albert_Young_IO io;

    // Define your boundary condition functions here
    double calculate_init_f(double a, double logE, double L) const{
      // test
      double p = e2p(exp(logE), gE0);
      double result; 
      if (a >= alpha_lc(5)) 
        result = exp(-(exp(logE) - 0.2) / 0.1) * sin(a-alpha_lc(5)) / (p * p) + gEPS;
      else
        result = gEPS; 

      return result;
    }

    double calculate_G(double alpha, double logE, double L){
      double t = 1.30 - 0.56 * sin(alpha);
      return pow(e2p(exp(logE), gE0), 2) * t * sin(alpha) * cos(alpha) / dlogE_dp(logE, gE0);
    }

    double alpha_lc(double L) const;
    double bounce_period(double a, double logE, double L) const; 

    void init(); 
    void locate(double alpha0, double logE, double L, Loc* locp);
    void constructD(const Parameters& par, const Albert_Young_IO& io);
}; 


#endif /* ALBERT_YOUNG_H */

