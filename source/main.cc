/*
 * fin:        main.cc
 * Author:      Xin Tao <xtao@ustc.edu.cn>
 *              Peng Peng <pp140594@mail.ustc.edu.cn>
 * Date:        05/12/2024 
 * 
 * Copyright (c) Xin Tao 
 *
 */

#include <fstream>
#include <iostream>
#include "Parameters.h"
#include "Mesh.h"
#include "Albert_Young.h"
#include "Storm.h"
#include "Radial.h"
#include "Solver.h"
#include <ctime>
#include "omp.h"


int main(int argc, char** argv) {

  Parameters paras(argc,argv); 
 
  // Create mesh 
  Mesh m(paras);

  // Create diffusion coefficients object
  Albert_Young eq(paras, m); 
  // Storm eq(paras, m); 
  // Radial eq(paras, m);

  Solver solver(paras, m, &eq);

  std::string filename;
  std::ofstream out; 

  // Create output H5 file
  filename = paras.output_path() + "/" + paras.run_id() + "_data.h5";
  HighFive::File file(filename, HighFive::File::Overwrite);

  const Eigen::VectorXd alpha0 = m.x() / gPI * 180.0;
  const Eigen::VectorXd logEN = m.y().array() - log(gE0);

  xt::dump(file, "/alpha0", alpha0, xt::dump_mode::overwrite);
  xt::dump(file, "/logEN", logEN, xt::dump_mode::overwrite);
  xt::dump(file, "/L", m.z(), xt::dump_mode::overwrite);

  // The timer
  clock_t start, end;
  double cpu_time;
  start = clock();

  xt::dump(file, "/f/0", solver.f(), xt::dump_mode::overwrite);

  // Time loop for solving
  for (int tstep = 1; tstep <= paras.nsteps(); ++tstep) {

    // Solve using FVM solver
    solver.update();

    if(tstep % paras.save_every_step() == 0){
      xt::dump(file, "/f/" + std::to_string(int((tstep) / paras.save_every_step())), solver.f(), xt::dump_mode::overwrite);
    }
  }

  Xarray1d t = xt::linspace(0.0, paras.T(), paras.nplots()+1);

  xt::dump(file, "/t", t, xt::dump_mode::overwrite);

  end = clock();
  cpu_time = ((double) (end - start)) / CLOCKS_PER_SEC;
  std::cout << "CPU time used " << cpu_time << " seconds" << std::endl;

  return 0;
}

