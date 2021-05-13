//
//  main.cpp
//  DumpConfig
//  Planting equilibrium configurations by eigenvector of transfer matrix method
//
//  Created by Yi Hu on 8/5/20.
//

#include <cstdlib>
#include <string>
#include <chrono>
#include <iostream>

#include "../TMat/read_input.hpp"
#include "../TMat/commondef.hpp"
#include "utilityb.hpp"
#include "dumpconfig.hpp"
#include "dumpreduce.hpp"

using clock_type = std::chrono::time_point<std::chrono::steady_clock>;
using duration_type = std::chrono::duration<double>;

int test_main()
{
  
  return 0;
}

// Usage: DumpConfig filename, num_of_file length_of_config
int compute_main(int argc, const char * argv[])
{
  const char *infile, *dftfname = "input.dat";
  int ndump=1, length=100;
  if(argc > 1)
  {
    if(argc > 2)
    {
      if(argc > 3)
      {
        length = atoi(argv[3]);
      }
      ndump = atoi(argv[2]);
    }
    
    infile = argv[1];
  }
  else
  {
    infile = dftfname;
  }
  MyUtility::setseed();
  Read_Input input;
  int err = input.read(infile);
  if(err) return err;
  clock_type t1, t2;
  t1 = std::chrono::steady_clock::now();
  // get parameters
  const MParameter mpara{input.betaJ, input.kappa, input.kappa2, input.h, input.size};
  // dump config
  DumpMat::DumpConfig_Base *bs = nullptr;
  const int multieigs = input.neigs>1? 2: 0;
  if(2 == input.dim)
  {
    if('A' == input.model)
    {
      if(1 == input.axialdir)
      {
        bs = new DumpMat::XAT2DDumpConfig(mpara,
               std::string(input.outprefix)+"_vec.dat");
      }
      else if(2 == input.axialdir)
      {
        bs = new DumpMat::YAT2DDumpConfig(mpara,
               std::string(input.outprefix)+"_vec.dat",
               std::string(input.outprefix)+"_vecl.dat");
      }
      else if(5 == input.axialdir)
      {
        bs = new DumpMat::XAM2DDumpConfig(mpara,
               std::string(input.outprefix)+"_vec.dat", multieigs);
      }
      else if(6 == input.axialdir)
      {
        bs = new DumpMat::YAM2DDumpConfig(mpara,
               std::string(input.outprefix)+"_vec.dat",
               std::string(input.outprefix)+"_vecl.dat");
      }
      else
      {
        std::cerr << "Error: axial direction not supported yet" << std::endl;
      }
    }
    else if('B' == input.model)
    {
      if(6 == input.axialdir)
      {
        bs = new DumpMat::YBM2DDumpConfig(mpara,
                std::string(input.outprefix)+"_vec.dat",
                std::string(input.outprefix)+"_vecl.dat");
      }
      else if(7 == input.axialdir)
      {
        bs = new DumpMat::ZBM2DDumpConfig(mpara,
                std::string(input.outprefix)+"_vec.dat",
                std::string(input.outprefix)+"_vecl.dat");
      }
      else
      {
        std::cerr << "Error: axial direction not supported yet" << std::endl;
      }
    }
    else if('D' == input.model)
    {
      if(5 == input.axialdir)
      {
        bs = new DumpMat::XDM2DDumpConfig(mpara,
               std::string(input.outprefix)+"_vec.dat", multieigs);
      }
      else
      {
        std::cerr << "Error: axial direction not supported yet" << std::endl;
      }
    }
    else
    {
      std::cerr << "Error: model not supported yet" << std::endl;
    }
  }
  else
  {
    std::cerr << "Error: dimension not supported yet" << std::endl;
  }
  
  if(bs)
  {
    bs->dump(ndump, length, input.outprefix);
    delete bs;
  }
  t2 = std::chrono::steady_clock::now();
  auto time_consumed = static_cast<duration_type>(t2-t1).count();
  if(time_consumed>=1.0)
    std::cout << "Time consumed: " << time_consumed << " s" << std::endl;
  return 0;
}

int main(int argc, const char * argv[])
{
  int ret_val =  compute_main(argc, argv);
  return ret_val;
}
