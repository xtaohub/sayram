/*
 * File:        Mesh.h
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 *              Peng Peng <pp140594@mail.ustc.edu.cn>
 * Date:        05/12/2024 
 * 
 * Copyright (c) Xin Tao 
 *
 */

#ifndef MESH_H_
#define MESH_H_

#include <vector>
#include "common.h"
#include "Parameters.h" 
#include "Face.h"

struct Ind{
  std::size_t i;
  std::size_t j;
  std::size_t k;
}; 

class Mesh {
  public:
    Mesh(const Parameters& p); 

    const Eigen::VectorXd& x() const { return x_; }
    const Eigen::VectorXd& y() const { return y_; }
    const Eigen::VectorXd& z() const { return z_; }

    double x(int i) const { return x_(i); }
    double y(int j) const { return y_(j); }
    double z(int k) const { return z_(k); }

    double xO() const { return xO_; }
    double yO() const { return yO_; }
    double zO() const { return zO_; }

    std::size_t nx() const { return nx_; }
    std::size_t ny() const { return ny_; }
    std::size_t nz() const { return nz_; }

    double dx() const { return dx_; }
    double dy() const { return dy_; }
    double dz() const { return dz_; }
    double dt() const { return dt_; }

    double cell_volume() const { return dx() * dy() * dz(); }
    double cell_volume_dt() const { return cell_volume() / dt();}

    void indO(const Point& A, Ind* indp) const { 
      // calculate the i,j coordinate relative to the Origin
      // Note: not the cell index.
      // This function is useful to calculate fA and fB from interpolation
      indp->i = round((A(0) - xO()) / dx());
      indp->j = round((A(1) - yO()) / dy()); 
      indp->k = round((A(2) - zO()) / dz());
      
    }

    std::size_t flatten_index(const Ind& ind) const { // map 3d indices to 1
      return ind.k*nx()*ny() + ind.j*nx() + ind.i;
    }

    std::size_t nnbrs() const { return 6; } // each cell has 6 nbrs

    // define the neighbor # of six adjacent cells
    // im -- (i-1, j, k); jp -- (i, j+1, k)
    // ip -- (i+1, j, k); jm -- (i, j-1, k)
    // km -- (i, j, k-1); kp -- (i, j, k+1)
    int inbr_im() const { return 0; }
    int inbr_jp() const { return 1; }
    int inbr_ip() const { return 2; }
    int inbr_jm() const { return 3; }
    int inbr_km() const { return 4; }
    int inbr_kp() const { return 5; }

    int rinbr(int inbr) const { return rinbr_(inbr); }                                    

    void get_nbr_ind(const Ind& ind, int inbr, Ind* nbr_indp) const {
      *nbr_indp = nbr_inds(ind.i,ind.j,ind.k,inbr); 
    }

    void get_nbr_face(const Ind& ind, int inbr, Face* facep) const {
      *facep = faces(ind.i,ind.j,ind.k,inbr); 
    }

  private:
    std::size_t nx_; 
    std::size_t ny_; 
    std::size_t nz_; 
    double dx_;
    double dy_;
    double dz_;
    double dt_; 

    // coordinate origin: corresponds to i-0.5, j-0.5
    double xO_; 
    double yO_;
    double zO_;

    Eigen::VectorXd x_; 
    Eigen::VectorXd y_;
    Eigen::VectorXd z_;

    Eigen::Matrix<int, 6, 1> rinbr_; 

    xt::xtensor<Ind,4> nbr_inds;
    xt::xtensor<Face,4> faces;

    void build_connectivity();
};

#endif /* MESH_H */

