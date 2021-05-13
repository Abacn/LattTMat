//
//  tests.cpp
//  XTMat
//
//  Created by Yi Hu on 1/17/20.
//

#include <iostream>
#include <fstream>
#include <Spectra/SymEigsSolver.h>
#include "tmat.hpp"
#include "tmatreduce.hpp"
#include "utility.hpp"
#include "read_input.hpp"
#include "pfunc.hpp"

int test_gentmat()
{
  // test generating whole transform matrix
  int size = 10;
  double J = 2.5, kappa = 1.0;
  auto result = TMat::genXATmat2D({J, kappa, 0., 0., size});
  std::ofstream ofs("fulltmat.dat");
  ofs << result << std::endl;
  return 0;
}

int test_geninternaltmat()
{
  // test generating whole transform matrix
  MParameter paras{1.3862943611, 0.5, 0., 4};
  //TMat::XATMatVecProd2D tmat(paras);
  TMat::XDMMatVecProd2D tmat(paras);
  std::ofstream ofs("internaltmatA.dat", std::ios::binary);
  uint64_t nstate = tmat.rows();
  auto const fullmat = tmat.raw_Tmat();
  ofs.write((char*) fullmat.data(), nstate*nstate*sizeof(typename decltype(fullmat)::Scalar));
  ofs.close();
  ofs.open("internaltmatB.dat", std::ios::binary);
  auto const fullmatt = TMat::genXDTmat2D(paras);
  nstate = fullmatt.rows();
  ofs.write((char*) fullmatt.data(), nstate*nstate*sizeof(typename decltype(fullmatt)::Scalar));
  ofs.close();
  return 0;
}

int test_calceigfromfile(const char *filename)
{
  std::ifstream in(filename, std::ifstream::ate | std::ifstream::binary);
  int nstate = (int)sqrt((double)(in.tellg()));
  in.close();
  Eigen::MatrixXd tmat(nstate, nstate);
  std::ifstream ifs(filename, std::ifstream::binary);
  ifs.read( (char *)tmat.data() , nstate*nstate*sizeof(typename Eigen::MatrixXd::Scalar) );
  Spectra::DenseSymMatProd<double> op(tmat);
  Spectra::SymEigsSolver<double, Spectra::LARGEST_MAGN, Spectra::DenseSymMatProd<double> > eigs(&op, 1, 3);
  eigs.init();
  long nconv = eigs.compute();
  std::cout << "Return status: " << nconv << "\n" << "Eigenvalue: " << eigs.eigenvalues() << std::endl;
  return 0;
}

int test_genmatprod()
{
  // test Transfer matrix vector multiplication class
  int size = 10;
  double J = 2.5, kappa = 1.0;
  double dtmp, apdv = 0.;
  auto whole_tmat = TMat::genXATmat2D({J, kappa, 0.,  0., size});
  SPIN_TYPE nstate = 1ULL << size;
  TMat::XATMatVecProd2D mvp({J, kappa, 0., 0., size});
  Eigen::MatrixXd whole_tmat2 = mvp.raw_Tmat();
  for(SPIN_TYPE rp=0; rp<nstate; ++rp)
  {
    for(SPIN_TYPE rq=0; rq<nstate; ++rq)
    {
      dtmp = whole_tmat2(rq, rp) - whole_tmat(rq, rp);
      apdv += dtmp*dtmp;
    }
  }
  //std::cout << "TMat A:\n" << whole_tmat;
  //std::cout << "\n\nTMat B:\n" << whole_tmat2;
  std::cout << "\nSq dev per row: " << apdv/nstate;
  std::cout <<std::endl;
  return 0;
}

int test_genxammat()
{
  // test generating whole reduced transform matrix
  MParameter paras{0.693147180559945, 0, 0., 2}; // ln 2
  //TMat::XATMatVecProd2D tmat(paras);
  TMat::Y2ATMatVecProd2D tmat(paras);
  std::ofstream ofs("internaltmatA.dat", std::ios::binary);
  uint64_t nstate = tmat.rows();
  auto const fullmat = tmat.raw_Tmat();
  ofs.write((char*) fullmat.data(), nstate*nstate*sizeof(typename decltype(fullmat)::Scalar));
  ofs.close();
  return 0;
}

int test_2dbitoperation()
{
  // test 2d bit opearation
  SPIN_TYPE input = 0x000000001248edb7ull;
  std::cout << std::hex;
  std::cout << "number\t" << input << "\n"
  << "rffset\t" << Spin::rffsetState_3D(input, 8, 8, 1) << "\n"
  << "dffset\t" << Spin::dffsetState_3D(input, 8, 8, 1) << "\n"
  << "reflect\t" << Spin::reflectState_3D(input, 8, 8, 0) << "\n"
  << "deflect\t" << Spin::reflectState_3D(input, 8, 8, 1) << "\n"
  << "rotate\t" << Spin::rotateState_3D(input, 8) << std::endl;
  return 0;
}

// check eigenvector
int test_chkeigvec()
{
  Read_Input input;
  int err = input.read("input_zbnnni.dat");
  if(err) return err;
  pfunc(input);
  // then read from vectors
  TMat::ZBTMatVecProd2D tmat({input.betaJ, input.kappa, input.kappa2, input.h, input.size});
  SPIN_TYPE niter = tmat.cols();
  double *vec = new double[niter], *vecl = new double[niter];
  std::ifstream ifs( "2DZBT_vec.dat", std::ios::binary );
  ifs.read((char*)(vec), sizeof(double)*niter );
  ifs.close();
  ifs.open( "2DZBT_vecl.dat" );
  ifs.read((char*)(vecl), sizeof(double)*niter );
  ifs.close();
  // output the probability of chekerboard and diagonal
  SPIN_TYPE cbspin = 0, dgspin = 0, dgspin2 = 0, dgspin3 = 0;
  for(int rp = input.size; rp>0; rp -= 4)
  {
    dgspin2 <<= 4;
    dgspin2 |= 5;
    dgspin3 <<= 4;
    dgspin3 |= 3;
    cbspin <<= 4;
    cbspin |= 1;
  }
  std::cout << "probability checkerboard, diag1, diag2: \n";
  std::cout <<  vec[cbspin] << '*' << vecl[cbspin] << '=' << vec[cbspin]*vecl[cbspin] << std::endl;
  std::cout <<  vec[dgspin] << '*' << vecl[dgspin] << '=' << vec[dgspin]*vecl[dgspin] << std::endl;
  std::cout <<  vec[dgspin2] << '*' << vecl[dgspin2] << '=' << vec[dgspin2]*vecl[dgspin2] << std::endl;
  std::cout <<  vec[dgspin3] << '*' << vecl[dgspin3] << '=' << vec[dgspin3]*vecl[dgspin3] << std::endl;
  delete [] vec;
  delete [] vecl;
   /*
  //print largest
  for(SPIN_TYPE rp=0; rp<niter; ++rp)
  {
    if(vec[rp]>0.1 || vecl[rp]>0.1)
    {
      std::cout << rp << " " << vec[rp] << " " << vecl[rp] << std::endl;
    }
  }*/
  return 0;
}

// check reduced eigenvector
int test_chkreducedeigvec()
{
  Read_Input input;
  int err = input.read("input_ybmnnni.dat");
  if(err) return err;
  pfunc(input);
  // then read from vectors
  TMat::YBMMatVecProd2D tmat({input.betaJ, input.kappa, input.kappa2, input.h, input.size});
  SPIN_TYPE niter = tmat.cols();
  double *vec = new double[niter], *vecl = new double[niter];
  std::ifstream ifs( "2DYBM_vec.dat", std::ios::binary );
  ifs.read((char*)(vec), sizeof(double)*niter );
  ifs.close();
  ifs.open( "2DYBM_vecl.dat" );
  ifs.read((char*)(vecl), sizeof(double)*niter );
  ifs.close();
  
  // output the probability of chekerboard and diagonal
  SPIN_TYPE cbspin = 0;
  for(int rp = input.size; rp>0; rp -= 4)
  {
    cbspin <<= 4;
    cbspin |= 3;
  }
  SPIN_TYPE dgspin = (cbspin << input.size) | Spin::offsetState(cbspin, input.size, 1);
  cbspin |= (cbspin << input.size);
  SPIN_TYPE dgidx = tmat.getIrrepStateIdx(dgspin), cbidx = tmat.getIrrepStateIdx(cbspin);
  std::cout << "probability Diagonal/checkerboard: \n";
  std::cout <<  vec[dgidx]*vecl[dgidx] << "/" << vec[cbidx]*vecl[cbidx] << "=" <<  (vec[dgidx]*vecl[dgidx])/(vec[cbidx]*vecl[cbidx]) << std::endl;
  delete [] vec;
  delete [] vecl;
  return 0;
}
