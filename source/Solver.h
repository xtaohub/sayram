/*
 * File:        Solver.h
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 *              Peng Peng <pp140594@mail.ustc.edu.cn>
 * Date:        05/12/2024 
 * 
 * Copyright (c) Xin Tao 
 *
 */

#ifndef SOLVER_H
#define SOLVER_H

#include "common.h"
#include "Mesh.h"
#include "Equation.h"
#include "Parameters.h"

class Solver {
  public:
    Solver(const Parameters& paras_in, const Mesh& m_in, Equation* eqp);

    void update();
    double t() const { return istep_ * m.dt(); }
    const Xtensor3d& f() const { return f_; }
    double f(const Ind& ind) const { return f_(ind.i, ind.j, ind.k); }

  private:
    const Parameters& paras; 
    const Mesh& m;

    Equation& eq; 

    std::size_t istep_; // used to calcualte time = istep_ * dt() 

    // M f = R
    Eigen::BiCGSTAB<SpMat> iterSolver;

    SpMat M_;
    std::vector<T> M_coeffs_; 

    Xtensor3d f_;
    Eigen::VectorXd R_;
    Eigen::VectorXd ftmp_; // used to store results from the Eigen Solver.

    xt::xtensor<Eigen::Matrix3d, 3> Lambda_;  // Lambda_(i,j,k) would be the Lambda Matrix at cell (i,j,k) 
    const Eigen::Matrix3d& Lambda(int i, int j, int k) const { return Lambda_(i,j,k); }
    const Eigen::Matrix3d& Lambda(const Ind& ind) const { return Lambda_(ind.i,ind.j,ind.k); }
    void update_Lambda();                                            

    //
    // f at vertices to build a lookup table
    // 
    Xtensor3d vertex_f_; 

    double vertex_f(int i, int j, int k) const { return vertex_f_(i,j,k); }
    double vertex_f(const Ind& ind) const { return vertex_f_(ind.i, ind.j, ind.k); }
    void update_vertex_f(); 

    void assemble();

    // a_sigma_i = |sigma| * n_{sigma} * Lambda * beta_{sigma,i}
    void a_sigma_i_func(const Ind& ind, const Face& face, Eigen::Vector4d* a_sigma_i_p);  
                                                                                                                        
    // a_sigma = sum (a_sigma_i * fi)
    // a_sigma_i_sum = sum(a_sigma_i)
    void a_sigma_func(const Ind& ind, int inbr, double* a_sigmap, double* a_sigma_i_sump); 


    // update coefficients in M for cell Ind, and its neighbor. 
    // note that coefficients for both cell Ind and its neighbor are updated. 
    void update_coeff_inner_pair(const Ind& ind, int inbr);  

    // Here: the inbr neighbor is a Dirichlet boundary face.
    void update_coeff_dirbc(const Ind& ind, int inbr); 

    void calculate_mu(double a_sigma_K, double a_sigma_L, double* mu_Kp, double* mu_Lp){
      double denom = std::abs(a_sigma_K) + std::abs(a_sigma_L) + 2*gEPS; 
      (*mu_Kp) = (std::abs(a_sigma_L)+gEPS)/denom; 
      (*mu_Lp) = 1 - (*mu_Kp); 
    }

    void init();

};

#endif /* SOLVER_H */
