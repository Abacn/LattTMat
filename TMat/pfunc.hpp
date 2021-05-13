//
//  pfunc.hpp
//  TMat
//
//  Created by Yi Hu on 2/12/20.
//

#ifndef pfunc_hpp
#define pfunc_hpp

#include <fstream>
#include <iomanip>
#include <cmath>
#include <type_traits>
#include <Eigen/Core>
#include <Spectra/SymEigsSolver.h>
#include <Spectra/GenEigsSolver.h>
#include "read_input.hpp"
#include "utility.hpp"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"

/** get model parameters from input class
 */
MParameter get_MParameter(const Read_Input &input);

// if constexpr alternative
template <typename MatVecProdType, bool B>
struct moment_dummy;

// symmetric
template <typename MatVecProdType>
struct moment_dummy<MatVecProdType, true> {
  /** @param op matrix vector multiplication operation instance
   *  @param q store leading the eigenvector. Will new memory if needed
   *  @param vl store leading eigenvalue
   */
  static void getlefteigvec(MatVecProdType &op, typename MatVecProdType::Scalar* &q, typename MatVecProdType::Scalar* vl)
  {
    q = nullptr;
    if(vl != nullptr) *vl = NAN;
  }
};

// general
template <typename MatVecProdType>
struct moment_dummy<MatVecProdType, false> {
  static void getlefteigvec(MatVecProdType &op, typename MatVecProdType::Scalar* &q, typename MatVecProdType::Scalar* vl)
  {
    typename MatVecProdType::MatVecProd_T t_op = op.get_t_instance(); // transpose matrix vector multiplication instance
    Spectra::GenEigsSolver<double, Spectra::LARGEST_MAGN, typename MatVecProdType::MatVecProd_T> teigs(&t_op, 1, 7);
    teigs.init();
    long ncv = teigs.compute();
    if(ncv == 0)
    {
      std::cout << "Warning: left eigenvalue calculation failed." << std::endl;
      q = nullptr;
      if(vl != nullptr) *vl = NAN;
    }
    else
    {
      q = new typename MatVecProdType::Scalar[t_op.rows()];
      std::copy_n(teigs.eigenvectors().data(), t_op.rows(), q);
      if(vl != nullptr) *vl = teigs.eigenvalues()(0);
    }
  }
};

template <typename MatVecProdType>
void getlefteigvec(MatVecProdType &op, typename MatVecProdType::Scalar* &q, typename MatVecProdType::Scalar* vl=nullptr)
{
  moment_dummy<MatVecProdType, std::is_same<typename MatVecProdType::Scalar, double>::value>::getlefteigvec(op, q, vl);
}

/** template partition function solver
 @tparam MatVecProdType matrix-vector multiplication class
 */
template <class MatVecProdType>
void pfunc_inner(const Read_Input &input, MatVecProdType &&op)
{
  using Scalar = typename MatVecProdType::Scalar;
  constexpr bool issymtmat = std::is_same<Scalar, double>::value;
  using EigenSolver = typename std::conditional<issymtmat,
  Spectra::SymEigsSolver<double, Spectra::LARGEST_MAGN, MatVecProdType>,
  Spectra::GenEigsSolver<double, Spectra::LARGEST_MAGN, MatVecProdType> >::type;
  int maxtrial = 3, ntrial = 0;
  int itpara = op.suggested_nev(input.neigs); // parameter controlling algorithm convergence rate
  // very small matrix
  if(itpara >= op.cols()) itpara = op.cols()-1;
  EigenSolver* eigs;
  long ncv = 0; // number of converged eigenvalues
  while(1)
  {
    try {
      eigs = new EigenSolver(&op, input.neigs, itpara);
      eigs->init();
      ncv = eigs->compute();
    } catch (std::exception &e) {
      std::cerr << "Catched " << typeid( e ).name( ) << std::endl;
      std::cerr << "what()  " << e.what( ) << std::endl;
    }
    if(ncv == input.neigs) break;
    ++ntrial;
    std::cout << "eigensolver not converged (" << ntrial << "/" << maxtrial << ")" << std::endl;
    if(ntrial == maxtrial) break;
    itpara += 2;
    delete eigs;
  }
  if(ncv==0)
  {
    std::cerr << "Number of converged eigenvalue is zero" << std::endl;
    exit(-1);
  }
  // save eigen values
  Scalar eigvl = eigs->eigenvalues()(0);
  std::ofstream ofs(input.outprefix + "_val.dat");
  ofs << std::setprecision(16) << eigs->eigenvalues()/op.getnormfactor();
  ofs.close();
  // save eigen vectors
  auto eigvecs = eigs->eigenvectors();
  delete eigs; // release memory of eigensolver
  eigs = nullptr;
  if(input.savevec[0])
  {
    //std::cout << eigvs.data().rows() << "\t" << eigvs.data().cols() << std::endl;
    std::ofstream ofv(input.outprefix + "_vec.dat", std::ios::binary);
    if(1==ncv && !issymtmat)
    {
      // leading eigenvector is real
      const Scalar *rtmpvec = eigvecs.data();
      for(uint64_t rp=0; rp<op.rows(); ++rp)
      {
        ofv.write((char*)(rtmpvec+rp), sizeof(double));
      }
    }
    else
    {
      ofv.write((char*)(eigvecs.data()), ncv*op.rows()*sizeof(Scalar));
    }
    ofv.close();
  }
  if(input.savevec[1])
  {
    // physical properties
    // right eigenvector setup
    Scalar *qr = new Scalar[op.rows()];
    std::copy_n(eigvecs.data(), op.rows(), qr);
    eigvecs.resize(0,0);  // release memory for whole eigenvectors
    // left eigenvalue setup
    Scalar eigvlb;
    Scalar *tql = nullptr;
    getlefteigvec(op, tql, &eigvlb); // tql pass by reference
    // test if eigensolver consistent
    if(eigvlb==eigvlb) // not nan
    {
      if(input.savevec[0])
      {
        std::ofstream ofv(input.outprefix + "_vecl.dat", std::ios::binary);
        for(uint64_t rp=0; rp<op.rows(); ++rp)
        {
          ofv.write((char*)(tql+rp), sizeof(double));
        }
        ofv.close();
      }
      std::streamsize ss = std::cout.precision();
      std::cout.precision(10);
      std::cout << "Eigensolver check: " << std::abs(eigvl/eigvlb) << std::endl;
      std::cout.precision (ss);
    }
    const Scalar *ql = nullptr;
    if(tql == nullptr) ql = qr;
    else ql = tql;
    std::ofstream ofp(input.outprefix + "_phys.dat");
    // magnetization, susceptibility and higher order moments
    double mmoment[TMAT_NMOMENTS];
    op.moments(op.getLayerSize(), ql, qr, mmoment);
    // other observables
    std::vector<double> extend_obs;
    op.phys_extend(op.getLayerSize(), ql, qr, extend_obs);
    // free eigenvectors
    if(tql != nullptr)
    {
      delete [] tql;
      tql = nullptr;
    }
    delete [] qr;
    ofp << std::setprecision(10);
    for(int rp=0; rp<TMAT_NMOMENTS; ++rp)
    {
      if(rp>0) ofp << '\t';
      ofp << mmoment[rp];
    }
    for(int rp=0; rp<extend_obs.size(); ++rp)
    {
      if(rp==0) ofp << '\n';
      else if(rp>0) ofp << '\t';
      ofp << extend_obs[rp];
    }
    ofp.close();
  }
}

void pfunc(const Read_Input &input);

#pragma GCC diagnostic pop

#endif /* pfunc_hpp */
