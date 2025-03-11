/*
 * File:        Solver.cc
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 *              Peng Peng <pp140594@mail.ustc.edu.cn>
 * Date:        05/12/2024 
 * 
 * Copyright (c) Xin Tao 
 *
 */

#include "Solver.h"

Solver::Solver(const Parameters& paras_in, const Mesh& m_in, Equation* eqp)
  : paras(paras_in), m(m_in), eq(*eqp){

    std::size_t nx = m.nx();
    std::size_t ny = m.ny();
    std::size_t nz = m.nz();

    M_.resize(nx*ny*nz, nx*ny*nz);

    f_.resize({nx,ny,nz});

    R_.resize(nx*ny*nz);
    ftmp_.resize(nx*ny*nz);

    Lambda_.resize({nx, ny, nz}); 
    vertex_f_.resize({nx+1,ny+1,nz+1}); 

    istep_ = 0;

    init(); 

  }

void Solver::init(){

  iterSolver.setTolerance(1e-8);

  for (std::size_t i = 0; i < m.nx(); i++){
    for (std::size_t j = 0; j < m.ny(); j++){
      for (std::size_t k = 0; k < m.nz(); k++){
        f_(i,j,k) = eq.init_f({i,j,k});
      }
    }
  }

  update_Lambda();
  update_vertex_f();
}

void Solver::update_Lambda(){

  for (std::size_t i = 0; i < m.nx(); i++){
    for (std::size_t j = 0; j < m.ny(); j++){
      for (std::size_t k = 0; k < m.nz(); k++){
        Lambda_(i,j,k) << eq.Dxx({i,j,k}) * eq.G({i,j,k}), eq.Dxy({i,j,k}) * eq.G({i,j,k}), eq.Dxz({i,j,k}) * eq.G({i,j,k}),
                   eq.Dxy({i,j,k}) * eq.G({i,j,k}), eq.Dyy({i,j,k}) * eq.G({i,j,k}), eq.Dyz({i,j,k}) * eq.G({i,j,k}),
                   eq.Dxz({i,j,k}) * eq.G({i,j,k}), eq.Dyz({i,j,k}) * eq.G({i,j,k}), eq.Dzz({i,j,k}) * eq.G({i,j,k});

      }
    }
  }
}

void Solver::a_sigma_i_func(const Ind& ind, const Face& face, Eigen::Vector4d* a_sigma_i_p){  
  
  const int NT = face.NT;  // number of tetrahedrons per face; In our case, each face is divided into four tetrahedrons
  std::array<Vector3,NT> alpha1, alpha2, alpha3; // alpha_{K,sigma,1,2,3} in notes

  Point K = {m.x(ind.i), m.y(ind.j), m.z(ind.k)};
  Point O = {0,0,0}; 

  for (auto i=0; i<face.v_size(); ++i) O += face.v(i); 
  O = O / 4.0; 

  int ip1; 

  double VT = m.cell_volume()/24.0; // the volume of the tetrahedrons

  for (int i=0; i<NT; ++i) {

    ip1 = (i+1) % NT; 

    alpha1[i] = (face.v(i) - K).cross(face.v(ip1)-K) / (6*VT); 
    alpha2[i] = (face.v(ip1) - K).cross(O - K) / (6*VT);
    alpha3[i] = (O - K).cross(face.v(i)-K) / (6*VT); 
  }

  Vector3 av_alpha1 = {0,0,0};
  for (int i=0; i<NT; ++i) av_alpha1 += alpha1[i]; 
  av_alpha1 /= NT; 

  std::array<Vector3,NT> beta;
  for (int i=0; i<NT; ++i){
    beta[i] = (alpha2[i] + alpha3[(i+NT-1)%(NT)] + av_alpha1); 
    beta[i] = beta[i] / NT; 
  }
  
  for (int i=0; i<NT; ++i) {
    (*a_sigma_i_p)(i) = face.area() * (Lambda(ind)*face.n()).transpose() * beta[i]; 
  }

}

void Solver::a_sigma_func(const Ind& ind, int inbr, double* a_sigmap, double* a_sigma_i_sump){ 
  // a_sigma = sum (a_sigma_i * fi)
  // a_sigma_i_sum = sum (a_sigma_i)

  Face face;
  m.get_nbr_face(ind, inbr, &face); 

  Eigen::Vector4d a_sigma_i;
  a_sigma_i_func(ind, face, &a_sigma_i);

  double a_sigma = 0, a_sigma_i_sum = 0;
  Ind indv;
  for (auto iv = 0; iv < face.v_size(); ++iv) {
    m.indO(face.v(iv), &indv); 
    a_sigma += a_sigma_i(iv)*vertex_f(indv);
    a_sigma_i_sum += a_sigma_i(iv); 
  }

  *a_sigmap = a_sigma;
  *a_sigma_i_sump = a_sigma_i_sum; 
}

void Solver::update_coeff_inner_pair(const Ind& ind, int inbr){ // add coefficient from a inner (no-a-boundary) face  

  Face face; 
  Ind nbr_ind;

  int rinb;

  double a_sigma_K, a_sigma_L, a_sigma_i_sum_K, a_sigma_i_sum_L; 
  double mu_K, mu_L; 

  double B_sigma, B_sigma_abs, B_sigma_plus, B_sigma_minus;
  double A_K, A_L; 

  // calculate A_K
  a_sigma_func(ind, inbr, &a_sigma_K, &a_sigma_i_sum_K);

  // calculate A_L
  m.get_nbr_ind(ind, inbr, &nbr_ind); 
  rinb = m.rinbr(inbr);
  a_sigma_func(nbr_ind, rinb, &a_sigma_L, &a_sigma_i_sum_L); 

  calculate_mu(a_sigma_K, a_sigma_L, &mu_K, &mu_L);  

  B_sigma = mu_L*a_sigma_L - mu_K*a_sigma_K; 

  B_sigma_abs = std::abs(B_sigma); 
  B_sigma_plus = (B_sigma_abs + B_sigma) / 2.0;
  B_sigma_minus = (B_sigma_abs - B_sigma) / 2.0;

  A_K = mu_K * a_sigma_i_sum_K + B_sigma_plus/(f(ind)+gEPS); 
  A_L = mu_L * a_sigma_i_sum_L + B_sigma_minus/(f(nbr_ind)+gEPS); 

  std::size_t K = m.flatten_index(ind), L = m.flatten_index(nbr_ind); 

  M_coeffs_.push_back(T(K, K, A_K));
  M_coeffs_.push_back(T(K, L, -A_L));
  M_coeffs_.push_back(T(L, L, A_L)); 
  M_coeffs_.push_back(T(L, K, -A_K)); 

} 

void Solver::update_coeff_dirbc(const Ind& ind, int inbr){ // add coefficient from a Dirichlet face  

  Face face; 

  double a_sigma_K, a_sigma_i_sum_K; 

  double B_sigma, B_sigma_abs, B_sigma_plus, B_sigma_minus;
  double A_K; 

  // calculate A_K
  a_sigma_func(ind, inbr, &a_sigma_K, &a_sigma_i_sum_K);

  B_sigma = -a_sigma_K; 
  B_sigma_abs = std::abs(B_sigma); 
  B_sigma_plus = (B_sigma_abs + B_sigma) / 2.0;
  B_sigma_minus = (B_sigma_abs - B_sigma) / 2.0;

  A_K = a_sigma_i_sum_K + B_sigma_plus/(f(ind)+gEPS); 

  std::size_t K = m.flatten_index(ind);
  M_coeffs_.push_back(T(K, K, A_K));
  R_(K) += B_sigma_minus;

}


void Solver::assemble(){ // obtain M and R 

  std::size_t i, j, k;

  // i direction
  // put the handling of the boundary conditions at i==0 and i==nx-1 here.
  // i==0; by default, df/dx = 0 at (a=0), so do nothing.
  //
  // if, on the other hand, for some reason, the b.c. at i=0 is is set to f(alpha_LC) = 0
  // the Dirichlet type, remember to uncomment the following lines.
  
  // i = 0;
  // for (j=0; j<m.ny(); ++j)
  //   for (k=0; k<m.nz(); ++k)
  //     update_coeff_dirbc({i,j,k}, m.inbr_im());

  // i==m.nx()-1;
  // it's neumann type (df/dx=0) in our case, so do nothing.
     
  for (i=1; i<m.nx(); ++i)
    for (j=0; j<m.ny(); ++j)
      for (k=0; k<m.nz(); ++k)
          update_coeff_inner_pair({i,j,k}, m.inbr_im());

  // j direction
  // j == 0: handle dirbc
  j = 0;
  for (i=0; i<m.nx(); ++i)
    for (k=0; k<m.nz(); ++k)
      update_coeff_dirbc({i,j,k}, m.inbr_jm());

  for (i=0; i<m.nx(); ++i)
    for (j=1; j<m.ny(); ++j)
      for (k=0; k<m.nz(); ++k)
        update_coeff_inner_pair({i,j,k}, m.inbr_jm());

  // j == ny-1: handle dirbc
  j = m.ny()-1;
  for (i=0; i<m.nx(); ++i)
    for (k=0; k<m.nz(); ++k)
      update_coeff_dirbc({i,j,k}, m.inbr_jp());

  // k direction
  k = 0;
  for (i=0; i<m.nx(); ++i)
    for (j=0; j<m.ny(); ++j)
      update_coeff_dirbc({i,j,k}, m.inbr_km());

  for (i=0; i<m.nx(); ++i)
    for (j=0; j<m.ny(); ++j)
      for (k=1; k<m.nz(); ++k)
        update_coeff_inner_pair({i,j,k}, m.inbr_km());

  k = m.nz()-1; 
  for (i=0; i<m.nx(); ++i)
    for (j=0; j<m.ny(); ++j)
      update_coeff_dirbc({i,j,k}, m.inbr_kp());

  long K;
  double UKK;

  for (std::size_t i=0; i<m.nx(); ++i) {
    for (std::size_t j=0; j<m.ny(); ++j) {
      for (std::size_t k=0; k<m.nz(); ++k){
        K = m.flatten_index({i,j,k});
        UKK = eq.G({i,j,k}) * m.cell_volume_dt();
        M_coeffs_.push_back(T(K, K, UKK));
        R_(K) += UKK * f_(i,j,k);
      }
    }
  }

  M_.setFromTriplets(M_coeffs_.begin(), M_coeffs_.end());
}


void Solver::update() {
  R_.setZero();
  M_coeffs_.clear();

  assemble(); 

  iterSolver.compute(M_);
  ftmp_ = iterSolver.solve(R_);

  for (std::size_t i=0; i<m.nx(); ++i){
    for (std::size_t j=0; j<m.ny(); ++j) {
      for (std::size_t k=0; k<m.nz(); ++k){
        f_(i,j,k) = ftmp_(m.flatten_index({i, j, k})) * exp(-m.dt()/eq.tau({i,j,k}));
      }
    }
  }

  istep_ += 1;
  eq.update(t()); 
  update_Lambda();
  update_vertex_f();
}

void Solver::update_vertex_f(){

  for (std::size_t i=1; i<m.nx(); ++i){
    for (std::size_t j=1; j<m.ny(); ++j){
      for (std::size_t k=1; k<m.nz(); ++k){
        vertex_f_(i,j,k) = (f_(i-1,j-1,k-1) + f_(i-1,j-1,k) + f_(i-1,j,k) + f_(i,j,k) + 
            f_(i-1,j,k-1) + f_(i,j-1,k-1) + f_(i,j-1,k) + f_(i,j,k-1)) / 8.0; 
      }
    }
  }

  eq.apply_bcs(&vertex_f_);

}
