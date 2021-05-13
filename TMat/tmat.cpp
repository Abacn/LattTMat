//
//  xtmat.cpp
//  XTMat
//
//  Created by Yi Hu on 1/20/20.
//

#include "tmat.hpp"

#include <iostream>

namespace TMat
{
  // base class constructors
  XMatVecProd2D::XMatVecProd2D(int size_, SPIN_TYPE nstate_, int boundary_cond/*=0*/): SymMatPhys(boundary_cond), MatVecProd2D(size_, nstate_){}
  
  YMatVecProd2D::YMatVecProd2D(int size_, SPIN_TYPE nstate_, int boundary_cond/*=0*/): GenMatPhys(boundary_cond), MatVecProd2D(size_, nstate_){}

  ZMatVecProd2D::ZMatVecProd2D(int size_, SPIN_TYPE nstate_, int boundary_cond/*=0*/): GenMatPhys(boundary_cond), MatVecProd2D(size_, nstate_){}

  XMatVecProd3D::XMatVecProd3D(int size_, SPIN_TYPE nstate_, int boundary_cond/*=0*/): SymMatPhys(boundary_cond), MatVecProd3D(size_ >> 16, size_ & 0xffff, nstate_)
  {
    if(0==nstate) nstate = 1ull << size;
  }

  YMatVecProd3D::YMatVecProd3D(int size_, SPIN_TYPE nstate_, int boundary_cond/*=0*/): GenMatPhys(boundary_cond), MatVecProd3D(size_ >> 16, size_ & 0xffff, nstate_)
  {
    if(0==nstate) nstate = 1ull << (2*size);
  }

  SymMatPhys::SymMatPhys(int boundary_cond/*=0*/): MatPhys(boundary_cond) {}

  GenMatPhys::GenMatPhys(int boundary_cond/*=0*/): MatPhys(boundary_cond)
  {
    if(boundary_cond != 0) std::cerr << "Warning: open boundary for y-models have not been implemented yet" << std::endl;
  }
  /////////////////////////////////
  // 2D X-ANNNI matrix-vector multiplication class
  XATMatVecProd2D::XATMatVecProd2D(const MParameter &mpara, int boundary_cond/*=0*/): XMatVecProd2D(mpara.size, ((SPIN_TYPE)1) << mpara.size, boundary_cond),
  mpara(mpara), boundary_cond(boundary_cond), Tvec(genTvec({mpara.betaJ*0.5, mpara.kappa, 0., mpara.h, mpara.size}, boundary_cond)) // Tvec is Tv^(1/2)
  {
    bzwj1[0] = exp(mpara.betaJ);
    bzwj1[1] = 1. / bzwj1[0];
    bzwj2[0] = exp(-mpara.betaJ*mpara.kappa); // exp(-J*kappa)
    bzwj2[1] = 1. / bzwj2[0];
    // T(internal) = T*normfactor
    // normfactor is 1/ground state boltzmann weight
    if(mpara.kappa<0.5)
    {
      normfactor = exp(size*mpara.betaJ*(mpara.kappa-2.0));
    }
    else
    {
      normfactor = exp(-size*mpara.betaJ*(mpara.kappa+1.0));
    }
  }

  // One operation. y_out = (P*T1) * x_in
  void XATMatVecProd2D::permute_mult_T1(const double *x_in, double *y_out) const
  {
    // x_in and y_out: nstate vector
    for(SPIN_TYPE rp=0; rp<nstate; ++rp)
    {
      SPIN_TYPE rq = Spin::offsetState(rp, size, 1); // permuted idx
      SPIN_TYPE rpf = rp ^ 1; // flip last bit
      y_out[rq] = bzwj1[0]*x_in[rp] + bzwj1[1]*x_in[rpf];
    }
  }

  // y_out = M * x_in
  void XATMatVecProd2D::perform_op(const double *x_in, double *y_out) const
  {
    // Tv^(1/2) * x
#pragma omp parallel for
    for(SPIN_TYPE rp=0; rp<nstate; ++rp)
    {
      y_out[rp] = Tvec[rp]*x_in[rp];
    }
    // Th * Tv^(1/2) x
    for(int it=0; it<size; ++it)
    {
      SPIN_TYPE maskshift = 1<<it;
#pragma omp parallel for
      for(SPIN_TYPE rp=0; rp<nstate; ++rp)
      {
        if(rp & maskshift)
        {
          double dtmp = y_out[rp]; // 1 spin
          y_out[rp] = bzwj1[0]*dtmp + bzwj1[1]*y_out[rp ^ maskshift];
          y_out[rp ^ maskshift] = bzwj1[1]*dtmp + bzwj1[0]*y_out[rp ^ maskshift];
        }
      }
    }
    // after the loop y_p1 stores the result
    // Tv^(1/2) * Th Tv^(1/2) x
#pragma omp parallel for
    for(SPIN_TYPE rp=0; rp<nstate; ++rp)
    {
      y_out[rp] = normfactor*Tvec[rp]*y_out[rp];
    }
    OpCounter::inc();
  }

  
  /////////////////////////////////
  // 2D Y-ANNNI transpose matrix-vector multiplication class
  YATMatVecProd2D_T::YATMatVecProd2D_T(const MParameter &mpara): YMatVecProd2D(mpara.size, ((SPIN_TYPE)1) << (2*mpara.size)),
  mpara(mpara), Tvec(genTvec({mpara.betaJ, 0., 0., mpara.h, mpara.size})) // Tvec is Tv
  {
    bzwj1[0] = exp(mpara.betaJ);
    bzwj1[1] = 1. / bzwj1[0];
    bzwj2[0] = exp(-mpara.betaJ*mpara.kappa); // exp(-J*kappa)
    bzwj2[1] = 1. / bzwj2[0];
    // T(internal) = T*normfactor
    // normfactor is 1/ground state boltzmann weight
    if(mpara.kappa<0.5)
    {
      normfactor = exp(size*mpara.betaJ*(mpara.kappa-2.0));
    }
    else
    {
      normfactor = exp(-size*mpara.betaJ*(mpara.kappa+1.0));
    }
  }
  
  // One operation. y_out = T1^T * P^T * x_in
  void YATMatVecProd2D_T::permute_mult_T1(const double *x_in, double *y_out) const
  {
    // x_in and y_out: nstate vector
    SPIN_TYPE layermask = (1ull << size)-1ull;
    for(SPIN_TYPE rq=0; rq<nstate; ++rq)
    {
      SPIN_TYPE s1 = rq & layermask,
      s2 = rq >> size;
      SPIN_TYPE rp = (s2 << size) | Spin::lffsetState(s1, size, 1); // permuted idx
      SPIN_TYPE rpf = (s2 << size) | Spin::lffsetState(s1 ^ 1ull, size, 1); // flip last bit
      y_out[rq] = bzwj2[0]*x_in[rp] + bzwj2[1]*x_in[rpf];
    }
  }
  
  // y_out = M * x_in
  void YATMatVecProd2D_T::perform_op(const double *x_in, double *y_out) const
  {
    SPIN_TYPE layermask = (1ull << size)-1ull;
    double *y_tmp = new double[nstate];
    double *y_p1 = y_tmp, *y_p2 = y_out, *y_swap;
    // Tv^T x
    for(SPIN_TYPE rp=0; rp<nstate; ++rp)
    {
      //SPIN_TYPE s1 = rp >> size,
      //s2 = rp & layermask;
      //SPIN_TYPE rq = (s2 << size) | s1 ; // restore idx
      y_p1[rp] = Tvec[rp]*x_in[rp];
    }
    // Th^T * Tv^T x
    for(int rp=0; rp<size; ++rp)
    {
      permute_mult_T1(y_p1, y_p2);
      y_swap = y_p1;
      y_p1 = y_p2;
      y_p2 = y_swap;
    }
    // after the loop y_p1 stores the result
    // Tv^T * Th^T * x
    for(SPIN_TYPE rp=0; rp<nstate; ++rp)
    {
      SPIN_TYPE s1 = rp >> size,
      s2 = rp & layermask;
      SPIN_TYPE rq = (s2 << size) | s1 ; // restore idx
      y_out[rq] = normfactor*y_p1[rp];
    }
    OpCounter::inc();
    delete [] y_tmp;
  }
  
  /////////////////////////////////
  // 2D Y-ANNNI matrix-vector multiplication class
  YATMatVecProd2D::YATMatVecProd2D(const MParameter &mpara): YMatVecProd2D(mpara.size, ((SPIN_TYPE)1) << (2*mpara.size)),
  mpara(mpara), Tvec(genTvec({mpara.betaJ, 0., 0., mpara.h, mpara.size})) // Tvec is Tv
  {
    bzwj1[0] = exp(mpara.betaJ);
    bzwj1[1] = 1. / bzwj1[0];
    bzwj2[0] = exp(-mpara.betaJ*mpara.kappa); // exp(-J*kappa)
    bzwj2[1] = 1. / bzwj2[0];
    // T(internal) = T*normfactor
    // normfactor is 1/ground state boltzmann weight
    if(mpara.kappa<0.5)
    {
      normfactor = exp(size*mpara.betaJ*(mpara.kappa-2.0));
    }
    else
    {
      normfactor = exp(-size*mpara.betaJ*(mpara.kappa+1.0));
    }
  }
  
  // One operation. y_out = (P*T1) * x_in
  void YATMatVecProd2D::permute_mult_T1(const double *x_in, double *y_out) const
  {
    // x_in and y_out: nstate vector
    SPIN_TYPE layermask = (1ull << size)-1ull;
    for(SPIN_TYPE rp=0; rp<nstate; ++rp)
    {
      SPIN_TYPE s1 = rp >> size,
                s3 = Spin::offsetState(s1, size, 1),
                s2 = rp & layermask;
      SPIN_TYPE rq = (s3 << size) | s2 ; // permuted idx
      SPIN_TYPE rpf = rp ^ (1ull << size); // flip last bit
      y_out[rq] = bzwj2[0]*x_in[rp] + bzwj2[1]*x_in[rpf];
    }
  }
  
  // y_out = M * x_in
  void YATMatVecProd2D::perform_op(const double *x_in, double *y_out) const
  {
    SPIN_TYPE layermask = (1ull << size)-1ull;
    double *y_tmp = new double[nstate];
    double *y_p1 = y_tmp, *y_p2 = y_out, *y_swap;
    // x
    for(SPIN_TYPE rp=0; rp<nstate; ++rp)
    {
      y_p1[rp] = x_in[rp];
    }
    // Tv * x
    for(int rp=0; rp<size; ++rp)
    {
      permute_mult_T1(y_p1, y_p2);
      y_swap = y_p1;
      y_p1 = y_p2;
      y_p2 = y_swap;
    }
    // after the loop y_p1 stores the result
    // Th * Tv * x
    for(SPIN_TYPE rp=0; rp<nstate; ++rp)
    {
      SPIN_TYPE s1 = rp >> size,
      s2 = rp & layermask;
      SPIN_TYPE rq = (s2 << size) | s1 ; // restore idx
      y_out[rp] = normfactor*Tvec[rp]*y_p1[rq];
    }
    OpCounter::inc();
    delete [] y_tmp;
  }
  
  
  /////////////////////////////////
  // 2D X2-ANNNI matrix-vector multiplication class
  Y2ATMatVecProd2D::Y2ATMatVecProd2D(const MParameter &mpara): YMatVecProd2D(mpara.size, 1ull << 2*mpara.size),
  mpara(mpara), Tvec(genTvec({mpara.betaJ*0.5, mpara.kappa, 0., mpara.h, mpara.size})) // Tvec is Tv^(1/2)
  {
    bzwj1[0] = exp(mpara.betaJ);
    bzwj1[1] = 1. / bzwj1[0];
    bzwj2[0] = exp(mpara.betaJ*(-2.0*mpara.kappa)); // exp(-2*J*kappa)
    bzwj2[1] = 1. / bzwj2[0];
    // T(internal) = T*normfactor
    // normfactor is 1/ground state boltzmann weight
    if(mpara.kappa<0.5)
    {
      normfactor = exp(2*size*mpara.betaJ*(mpara.kappa-2.0));
    }
    else
    {
      normfactor = exp(-2*size*mpara.betaJ*(mpara.kappa+1.0));
    }
  }
  
  // One operation. y_out = (P*T1) * x_in
  void Y2ATMatVecProd2D::permute_mult_T1(const double *x_in, double *y_out) const
  {
    // x_in and y_out: nstate vector
    SPIN_TYPE highbit=1<<mpara.size, bothbit=highbit|1, mask=highbit-1;
    SPIN_TYPE rpleft, rpright, rq, jbit;
    for(SPIN_TYPE rp=0; rp<nstate; ++rp)
    {
      rpleft = rp&mask;
      rpright = rp>>mpara.size;
      rq = Spin::offsetState(rpright, size, 1)<<mpara.size | Spin::offsetState(rpleft, size, 1); // permuted idx
      jbit = (rpleft ^ rpright)&1;
      y_out[rq] = bzwj1[jbit]*(bzwj2[0]*x_in[rp] + x_in[rp ^ highbit]) +
      bzwj1[jbit ^ 1]*(bzwj2[1]*x_in[rp ^ bothbit] + x_in[rp ^ 1]);
    }
  }
  
  // y_out = M * x_in
  void Y2ATMatVecProd2D::perform_op(const double *x_in, double *y_out) const
  {
    double *y_tmp = new double[nstate];
    double *y_p1 = y_tmp, *y_p2 = y_out, *y_swap;
    // Tv^(1/2) * x
    for(SPIN_TYPE rp=0; rp<nstate; ++rp)
    {
      y_p1[rp] = Tvec[rp]*x_in[rp];
    }
    // Th * Tv^(1/2) x
    for(int rp=0; rp<size; ++rp)
    {
      permute_mult_T1(y_p1, y_p2);
      y_swap = y_p1;
      y_p1 = y_p2;
      y_p2 = y_swap;
    }
    // after the loop y_p1 stores the result
    // Tv^(1/2) * Th Tv^(1/2) x
    for(SPIN_TYPE rp=0; rp<nstate; ++rp)
    {
      y_out[rp] = normfactor*Tvec[rp]*y_p1[rp];
    }
    OpCounter::inc();
    delete [] y_tmp;
  }
  
  int Y2ATMatVecProd2D::moments(int size, const Scalar *ql, const Scalar *qr, double results[TMAT_NMOMENTS]) const
  {
    // TODO
    for(int rp=0; rp<TMAT_NMOMENTS; ++rp) results[rp] = 0.0;
    return -1;
  }
  
  
  /////////////////////////
  // 2D DNNI matrix-vector multiplication class
  XDTMatVecProd2D::XDTMatVecProd2D(const MParameter &mpara, int boundary_cond/*=0*/):
  XMatVecProd2D(mpara.size, ((SPIN_TYPE)1) << mpara.size),
  mpara(mpara), boundary_cond(boundary_cond),
  Tvec(genTvec({mpara.betaJ*0.5, mpara.kappa, 0., mpara.h, mpara.size}, boundary_cond)) // Tvec is Tv^(1/2)
  {
    bzwj1[0] = exp(mpara.betaJ);
    bzwj1[1] = 1. / bzwj1[0];
    bzwj2[0] = exp(-mpara.betaJ*mpara.kappa); // exp(-J*kappa)
    bzwj2[1] = 1. / bzwj2[0];
    // T(internal) = T*normfactor
    // normfactor is 1/ground state boltzmann weight
    if(mpara.kappa<0.5)
    {
      normfactor = exp(size*mpara.betaJ*2.0*(mpara.kappa-1.0));
    }
    else
    {
      normfactor = exp(-size*mpara.betaJ*2.0*mpara.kappa);
    }
  }
  
  // One operation. y_out = (P*T1) * x_in
  void XDTMatVecProd2D::permute_mult_T1(const double *x_in, double *y_out) const
  {
    SPIN_TYPE nstateb = nstate<<1;
    SPIN_TYPE mask = nstate-1, maskb = nstate - 2;
    // x_in and y_out: nstate vector
#pragma omp parallel for
    for(SPIN_TYPE rpb=0; rpb<nstateb; ++rpb) // row
    {
      SPIN_TYPE auxp = rpb>>size;
      SPIN_TYPE rp = rpb & mask;
      SPIN_TYPE rpp = auxp<<size | Spin::offsetState(rp, size, 1); // permuted idx
      SPIN_TYPE rq = (rp&maskb) | auxp;
      SPIN_TYPE sbit = rp&1;
      y_out[rpp] = bzwj1[sbit ^ auxp]*bzwj2[(rq>>1&1) ^ sbit]*(bzwj2[sbit]*x_in[rq] + bzwj2[sbit ^ 1]*x_in[rq | nstate]); // J term
    }
  }
  
  void XDTMatVecProd2D::permute_mult_T1_last(const double *x_in, double *y_out, SPIN_TYPE lastbit) const
  {
    SPIN_TYPE nstateb = nstate<<1;
    SPIN_TYPE mask = nstate-1, maskb = nstate - 2;
    // x_in and y_out: nstate vector
#pragma omp parallel for
    for(SPIN_TYPE rpb=0; rpb<nstateb; ++rpb) // row
    {
      SPIN_TYPE auxp = rpb>>size;
      SPIN_TYPE rp = rpb & mask;
      SPIN_TYPE rpp = auxp<<size | Spin::offsetState(rp, size, 1); // permuted idx
      SPIN_TYPE rq = (rp&maskb) | auxp;
      SPIN_TYPE sbit = rp&1;
      y_out[rpp] = bzwj1[sbit ^ auxp]*bzwj2[lastbit ^ sbit]*(bzwj2[sbit]*x_in[rq] + bzwj2[sbit ^ 1]*x_in[rq | nstate]); // J term
    }
  }
  
  // y_out = M * x_in
  void XDTMatVecProd2D::perform_op(const double *x_in, double *y_out) const
  {
    // introduce one auxiliary spin
    SPIN_TYPE nstateb = nstate<<1, nstatec = nstate>>1;
    double *y_tmp1 = new double[nstateb], *y_tmp2 = new double[nstateb];
    double *y_p1 = y_tmp1, *y_p2 = y_tmp2, *y_swap;
    if(0 == boundary_cond)
    {
      // Tv^(1/2) * x and with two auxiliary bits
      for(SPIN_TYPE lastbit: {0, 1})
      {
        y_p1 = y_tmp1;
        y_p2 = y_tmp2;
#pragma omp parallel for
        for(SPIN_TYPE rp=0; rp<nstate; ++rp)
        {
          SPIN_TYPE rpb = (rp&nstatec)<<1 | rp;
          if((rp&1)==lastbit)
          {
            y_p1[rpb] = Tvec[rp]*x_in[rp]; // Tvec[rp]*
          }
          else
          {
            y_p1[rpb] = 0.; // Tvec[rp]*
          }
          y_p1[rpb ^ nstate] = 0.;
        }
        for(int rp=0; rp<size-1; ++rp)
        {
          permute_mult_T1(y_p1, y_p2);
          y_swap = y_p1;
          y_p1 = y_p2;
          y_p2 = y_swap;
        }
        // last trun
        permute_mult_T1_last(y_p1, y_p2, lastbit);
        // after the loop y_p1 stores the result
        // Tv^(1/2) * Th Tv^(1/2) x
        if(0==lastbit)
        {
#pragma omp parallel for
          for(SPIN_TYPE rp = 0; rp < nstate; ++rp)
          {
            y_out[rp] = normfactor*Tvec[rp]*(y_p2[rp]+y_p2[rp ^ nstate]);
          }
        }
        else
        {
#pragma omp parallel for
          for(SPIN_TYPE rp = 0; rp < nstate; ++rp)
          {
            y_out[rp] += normfactor*Tvec[rp]*(y_p2[rp]+y_p2[rp ^ nstate]);
          }
        }
      }
    }
    else if(5==boundary_cond || 9==boundary_cond)
    {
      // constraint boundary conditions
      y_p1 = y_tmp1;
      y_p2 = y_tmp2;
      if(5 == boundary_cond)
      {
        // two up
#pragma omp parallel for
        for(SPIN_TYPE rp=0; rp<nstate; ++rp)
        {
          y_p1[rp] = Tvec[rp]*x_in[rp]*bzwj2[rp&1]*bzwj2[rp>>(size-1)];
          y_p1[rp | nstate] = 0.;
        }
      }
      else
      {
#pragma omp parallel for
        for(SPIN_TYPE rp=0; rp<nstate; ++rp)
        {
          y_p1[rp | nstate] = Tvec[rp]*x_in[rp]*bzwj2[(rp&1)^1]*bzwj2[rp>>(size-1)];
          y_p1[rp] = 0.;
        }
      }
      for(int rp=0; rp<size-1; ++rp)
      {
        permute_mult_T1(y_p1, y_p2);
        y_swap = y_p1;
        y_p1 = y_p2;
        y_p2 = y_swap;
      }
      permute_mult_T1_last(y_p1, y_p2, 0);
#pragma omp parallel for
      for(SPIN_TYPE rp = 0; rp < nstate; ++rp)
      {
        y_out[rp] = normfactor*Tvec[rp]*(y_p2[rp]+y_p2[rp | nstate]);
      }
    }
    OpCounter::inc();
    delete [] y_tmp1;
    delete [] y_tmp2;
  }
  
  /////////////////////////
  // 2D BNNNI (and 3NN) matrix-vector multiplication class, propagated in the diagonal direction
  ZBTMatVecProd2D::ZBTMatVecProd2D(const MParameter &mpara): ZMatVecProd2D(mpara.size, 1ull << mpara.size),
  mpara(mpara), Tvec(genTvec({mpara.betaJ*0.5, mpara.kappa, mpara.kappa2, mpara.h, mpara.size}))
  {
    bzwj1[0] = exp(mpara.betaJ);
    bzwj1[1] = 1. / bzwj1[0];
    bzwj2[0] = exp(-mpara.betaJ*mpara.kappa); // exp(-J*kappa)
    bzwj2[1] = 1. / bzwj2[0];
    bzwj3[0] = exp(-mpara.betaJ*mpara.kappa2); // exp(-J*kappa2)
    bzwj3[1] = 1. / bzwj3[0];
    // T(internal) = T*normfactor
    // normfactor is 1/ground state boltzmann weight
    double normfactorexp;
    if(mpara.kappa + mpara.kappa2 < 0.5)
    {
      // ferromagnatic
      normfactorexp = size*2.0*(mpara.kappa+mpara.kappa2-1.0);
    }
    else if(mpara.kappa < 0.5 && mpara.kappa2-mpara.kappa<0.5)
    {
      // <2>
      normfactorexp = -size;
    }
    else if(mpara.kappa2 < 2*mpara.kappa)
    {
      // checkerboard
      normfactorexp = -size*2.0*mpara.kappa;
    }
    else
    {
      // <1>
      normfactorexp = size*2.0*(mpara.kappa - mpara.kappa2);
    }
    double halfnormfactor = exp(0.5*mpara.betaJ*normfactorexp);
    normfactor = halfnormfactor*halfnormfactor;
    for(int rp=0; rp<nstate; ++rp)
    {
      Tvec(rp) *= halfnormfactor;
    }
  }

  // One operation. y_out = (P*T1) * x_in
  void ZBTMatVecProd2D::permute_mult_T1(const double *x_in, double *y_out, int optdirection) const
  {
    SPIN_TYPE nstateb = nstate<<2;
    SPIN_TYPE mask = nstate-1, maskb = nstate - 4; // 11..1, 11..00
    // x_in and y_out: nstate vector
    for(SPIN_TYPE rpb=0; rpb<nstateb; ++rpb) // row
    {
      SPIN_TYPE auxp = rpb>>size;
      SPIN_TYPE rp = rpb & mask;
      SPIN_TYPE rpp = auxp<<size | Spin::offsetState(rp, size, 2); // permuted idx
      SPIN_TYPE rq = (rp&maskb) | auxp;
      SPIN_TYPE slbit = rp>>1 & 1, srbit = rp & 1, // current left and right bits
                plbit = auxp>>1, prbit = auxp & 1, // previous left and right bits
                flbit = rp>>3 & 1, frbit = rp>>2 & 1; // two forward bits
      // get interaction energy not involved back bits
      if(!optdirection)
        y_out[rpp] = bzwj1[srbit^plbit] * bzwj2[srbit^frbit] * bzwj2[slbit^flbit] * bzwj3[srbit^prbit] * bzwj3[slbit^plbit];
      else
        y_out[rpp] = bzwj1[slbit^frbit] * bzwj1[slbit^prbit] * bzwj2[srbit^frbit] * bzwj2[slbit^flbit] * bzwj3[srbit^prbit] * bzwj3[slbit^plbit];
      // then combine with energy involved back bits
      double dtmp = 0.;
      for(SPIN_TYPE auxq=0; auxq<4; ++auxq)
      {
        SPIN_TYPE albit = auxq>>1, arbit = auxq&1;
        SPIN_TYPE rqq = auxq<<size | rq; // autq at leading bits
        if(!optdirection)
          dtmp += x_in[rqq]*bzwj1[srbit^albit]*bzwj2[srbit^arbit]*bzwj2[slbit^albit];
        else
          dtmp += x_in[rqq]*bzwj2[srbit^arbit]*bzwj2[slbit^albit];
      }
      y_out[rpp] *= dtmp; // J term
    }
  }
  
  // one operation at last run, use lastbit to replace flbit
  void ZBTMatVecProd2D::permute_mult_T1_last(const double *x_in, double *y_out, SPIN_TYPE lastbit, int optdirection) const
  {
    SPIN_TYPE nstateb = nstate<<2;
    SPIN_TYPE mask = nstate-1, maskb = nstate - 4; // 11..1, 11..00
    SPIN_TYPE flbit = lastbit>>1, frbit = lastbit & 1;
    // x_in and y_out: nstate vector
    for(SPIN_TYPE rpb=0; rpb<nstateb; ++rpb) // row
    {
      SPIN_TYPE auxp = rpb>>size;
      SPIN_TYPE rp = rpb & mask;
      SPIN_TYPE rpp = auxp<<size | Spin::offsetState(rp, size, 2); // permuted idx
      SPIN_TYPE rq = (rp&maskb) | auxp;
      SPIN_TYPE slbit = rp>>1 & 1, srbit = rp & 1, // current left and right bits
                plbit = auxp>>1, prbit = auxp & 1; // previous left and right bits
      // get interaction energy not involved back bits
      if(!optdirection)
        y_out[rpp] = bzwj1[srbit^plbit] * bzwj2[srbit^frbit] * bzwj2[slbit^flbit] * bzwj3[srbit^prbit] * bzwj3[slbit^plbit];
      else
        y_out[rpp] = bzwj1[slbit^frbit] * bzwj1[slbit^prbit] * bzwj2[srbit^frbit] * bzwj2[slbit^flbit] * bzwj3[srbit^prbit] * bzwj3[slbit^plbit];
      // then combine with energy involved back bits
      double dtmp = 0.;
      for(SPIN_TYPE auxq=0; auxq<4; ++auxq)
      {
        SPIN_TYPE albit = auxq>>1, arbit = auxq&1;
        SPIN_TYPE rqq = auxq<<size | rq; // autq at leading bits
        if(!optdirection)
          dtmp += x_in[rqq]*bzwj1[srbit^albit]*bzwj2[srbit^arbit]*bzwj2[slbit^albit];
        else
          dtmp += x_in[rqq]*bzwj2[srbit^arbit]*bzwj2[slbit^albit];
      }
      y_out[rpp] *= dtmp; // J term
    }
  }

  void ZBTMatVecProd2D::perform_op(const double *x_in, double *y_out) const
  {
    perform_op_internal(x_in, y_out, 0);
  }
  
  // matrix-vector multiplication
  void ZBTMatVecProd2D::perform_op_internal(const double *x_in, double *y_out, int optdirection) const
  {
    // introduce one auxiliary spin
    SPIN_TYPE nstateb = nstate<<2, nstatec = nstate>>1 | nstate>>2;
    double *y_tmp1 = new double[nstateb], *y_tmp2 = new double[nstateb];
    double *y_p1 = y_tmp1, *y_p2 = y_tmp2, *y_swap;
    // Tv^(1/2) * x and with two auxiliary bits
    for(SPIN_TYPE lastbit: {0, 1, 2, 3})
    {
      y_p1 = y_tmp1;
      y_p2 = y_tmp2;
      for(SPIN_TYPE rp=0; rp<nstate; ++rp)
      {
        SPIN_TYPE rpb = (rp&nstatec)<<2 | rp;
        if((rp&3)==lastbit)
        {
          y_p1[rpb] = Tvec[rp]*x_in[rp]; // Tvec[rp]*
        }
        else
        {
          y_p1[rpb] = 0.; // Tvec[rp]*
        }
        y_p1[rpb ^ nstate] = y_p1[rpb ^ (nstate<<1) ] = y_p1[rpb ^ (nstatec << 2)] = 0.;
      }
      for(int rp=0; rp<size/2-1; ++rp)
      {
        permute_mult_T1(y_p1, y_p2, optdirection);
        y_swap = y_p1;
        y_p1 = y_p2;
        y_p2 = y_swap;
      }
      // last trun
      permute_mult_T1_last(y_p1, y_p2, lastbit, optdirection);
      // after the loop y_p1 stores the result
      // Tv^(1/2) * Th Tv^(1/2) x
      if(0==lastbit)
      {
        for(SPIN_TYPE rp = 0; rp < nstate; ++rp)
        {
          y_out[rp] = Tvec[rp]*(y_p2[rp]+y_p2[rp ^ nstate]+y_p2[rp ^ (nstate<<1)]+y_p2[rp ^ (nstatec << 2)]);
        }
      }
      else
      {
        for(SPIN_TYPE rp = 0; rp < nstate; ++rp)
        {
          y_out[rp] += Tvec[rp]*(y_p2[rp]+y_p2[rp ^ nstate]+y_p2[rp ^ (nstate<<1)]+y_p2[rp ^ (nstatec << 2)]);
        }
      }
    }
    OpCounter::inc();
    delete [] y_tmp1;
    delete [] y_tmp2;
  }


  ZBTMatVecProd2D_T::ZBTMatVecProd2D_T(const ZBTMatVecProd2D &tmathandle):
  ZMatVecProd2D(tmathandle.size, 1ull << tmathandle.size),
  tmathandle(tmathandle)
  {
  }

  // matrix-vector multiplication
  void ZBTMatVecProd2D_T::perform_op(const double *x_in, double *y_out) const
  {
    tmathandle.perform_op_internal(x_in, y_out, 1);
  }


  //  Generate raw (unsymmetrized) TMat independent of matrix-vector multiplication. For internal and debug usage
  /** Generate unsymmetrized 2D x-ANNNI transfer matrix
   */
  Eigen::MatrixXd genXATmat2D(const MParameter &mpara)
  {
    SPIN_TYPE nstate = 1ULL << mpara.size;
    double J1 = -mpara.betaJ;
    Eigen::VectorXd V = XATMatVecProd2D::genTvec(mpara);
    Eigen::MatrixXd M(nstate, nstate);
    for(SPIN_TYPE nowstate = 0ULL; nowstate < nstate; ++nowstate)
    {
      for(SPIN_TYPE nextstate = 0ULL; nextstate < nstate; ++nextstate)
      {
        int n3anti = Spin::nAntiSpin(nowstate, nextstate);
        double bzn = Spin::spinBzExpn(mpara.size, n3anti, J1);
        M(nowstate, nextstate) = V(nowstate)*exp(bzn);
      }
    }
    return M;
  }
  
  /** Generate unsymmetrized 2D y-ANNNI transfer matrix
   */
  Eigen::MatrixXd genYATmat2D(const MParameter &mpara)
  {
    SPIN_TYPE nstate = 1ULL << 2*mpara.size, nlstate = 1ULL << mpara.size;
    double J2 = mpara.betaJ*mpara.kappa;
    Eigen::VectorXd V = YATMatVecProd2D::genTvec(mpara);
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(nstate, nstate);
    for(SPIN_TYPE nowstate = 0ULL; nowstate < nlstate; ++nowstate)
    {
      for(SPIN_TYPE centerstate = 0ULL; centerstate < nlstate; ++centerstate)
      {
        SPIN_TYPE rowidx = centerstate<<mpara.size | nowstate;
        for(SPIN_TYPE rightstate = 0ULL; rightstate < nlstate; ++rightstate)
        {
          SPIN_TYPE colidx = rightstate<<mpara.size | centerstate;
          int n3anti = Spin::nAntiSpin(nowstate, rightstate);
          double bzn = Spin::spinBzExpn(mpara.size, n3anti, J2);
          M(rowidx, colidx) = V(rowidx)*exp(bzn);
        }
      }
    }
    return M;
  }
  
  /** Generate unsymmetrized 2D x-DNNI transfer matrix
   */
  Eigen::MatrixXd genXDTmat2D(const MParameter &mpara)
  {
    SPIN_TYPE nstate = 1ULL << mpara.size;
    double J1 = -mpara.betaJ, J2 = mpara.betaJ*mpara.kappa;
    Eigen::VectorXd V = XDTMatVecProd2D::genTvec(mpara);
    Eigen::MatrixXd M(nstate, nstate);
    for(SPIN_TYPE nowstate = 0ULL; nowstate < nstate; ++nowstate)
    {
      for(SPIN_TYPE nextstate = 0ULL; nextstate < nstate; ++nextstate)
      {
        int n1anti = Spin::nAntiSpin(nowstate, nextstate),
            n2anti = Spin::nAntiSpin(nowstate, Spin::offsetState(nextstate, mpara.size, 1)),
            n3anti = Spin::nAntiSpin(nowstate, Spin::lffsetState(nextstate, mpara.size, 1));
        double bzn1 = Spin::spinBzExpn(mpara.size, n1anti, J1),
               bzn2 = Spin::spinBzExpn(mpara.size, n2anti, J2),
               bzn3 = Spin::spinBzExpn(mpara.size, n3anti, J2);
        // for raw matrix, V(nowstate); for symmetric matrix: sqrt(V(nowstate)*V(nextstate))
        M(nowstate, nextstate) = V(nowstate)*exp(bzn1+bzn2+bzn3);
      }
    }
    return M;
  }
}
