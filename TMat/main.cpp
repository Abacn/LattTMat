//
//  main.cpp
//  XTMat
//
//  Created by Yi Hu on 1/17/20.
//

#include <iostream>
#include <chrono>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "read_input.hpp"
#include "pfunc.hpp"

using clock_type = std::chrono::time_point<std::chrono::steady_clock>;
using duration_type = std::chrono::duration<double>;

int test_main()
{
  int test_gentmat(void), test_genmatprod(void), test_geninternaltmat(void), test_genxammat(void), test_2dbitoperation(void), test_chkeigvec(void), test_chkreducedeigvec(void);
  // test_gentmat();
  // test_geninternaltmat();
  // test_genmatprod();
  // test_genxammat();
  // test_2dbitoperation();
  test_chkeigvec();
  return 0;
}

int compute_main(int argc, const char * argv[])
{
  const char *infile, *dftfname = "input.dat";
  if(argc > 1)
  {
    infile = argv[1];
  }
  else
  {
    infile = dftfname;
  }
  Read_Input input;
  int err = input.read(infile);
  if(err) return err;
#ifdef _OPENMP
  if(omp_get_max_threads() > 1) std::cout << "Num of thread: " << omp_get_max_threads() << std::endl;
#endif
  clock_type t1 = std::chrono::steady_clock::now();
  pfunc(input);
  clock_type t2 = std::chrono::steady_clock::now();
  float time_consumed = static_cast<duration_type>(t2-t1).count();
  if(time_consumed>=1.0f)
    std::cout << "Time consumed: " << time_consumed << " s" << std::endl;
  int opcount = OpCounter::get();
  if(opcount > 0)
    std::cout << "Number of iteration: " << opcount << std::endl;
  return 0;
}

int main(int argc, const char * argv[])
{
  //test_main();
  compute_main(argc, argv);
  return 0;
}
