/*
 * File:        Albert_Young_IO.h
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 * Date:        11/27/2024 
 * 
 * Copyright (c) Xin Tao 
 *
 */

#ifndef ALBERT_YOUNG_IO_H_
#define ALBERT_YOUNG_IO_H_


#include "common.h"
#include "Parameters.h"

class Albert_Young_IO{
  public:
    Albert_Young_IO(const Parameters& paras);

    Xarray1d x_D; // angle(°) convert to radial
    Xarray1d y_D; // Energy(keV)
    Xarray1d z_D;

    Xtensor3d Dxx_raw;
    Xtensor3d Dxy_raw;
    Xtensor3d Dyy_raw;

    std::size_t nx_D() const { return nx_D_; }
    std::size_t ny_D() const { return ny_D_; }
    std::size_t nz_D() const { return nz_D_; }

    double xmin_D() const { return xmin_D_; }
    double xmax_D() const { return xmax_D_; }
    double ymin_D() const { return ymin_D_; }
    double ymax_D() const { return ymax_D_; }
    double zmin_D() const { return zmin_D_; }
    double zmax_D() const { return zmax_D_; }

    void update(double t) {} // do nothing

  private:

    const Parameters& paras; 

    std::size_t nx_D_, ny_D_, nz_D_;

    double xmin_D_, xmax_D_, ymin_D_, ymax_D_, zmin_D_, zmax_D_;

    void read_D(); 
}; 



#endif /* ALBERT_YOUNG_IO_H */

