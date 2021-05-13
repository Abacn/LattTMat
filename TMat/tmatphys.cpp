//
//  tmatphys.cpp
//  TMat
//
//  Created by Yi Hu on 4/17/20.
//

#include <iostream>
#include "tmat.hpp"

namespace TMat{
  // moments
  int SymMatPhys::moments(int size, const Scalar *q, double results[TMAT_NMOMENTS]) const
  {
    SPIN_TYPE nstate = 1ULL << size, rp;
    double dtmp, dPara, dmmt;
    int rq;
    //double rpara=0., ranti=0.;
    for(rp=0; rp<TMAT_NMOMENTS; ++rp)
    {
      results[rp] = 0.;
    }
    for(SPIN_TYPE rp=0ULL; rp<nstate; ++rp)
    {
      dmmt = dPara = Spin::spinBzexp(rp, size);
      dtmp = std::norm(q[rp]);
      for(rq=0; rq<TMAT_NMOMENTS; ++rq)
      {
        results[rq] += dmmt*dtmp;
        dmmt *= dPara;
      }
    }
    for(rq=0; rq<TMAT_NMOMENTS; ++rq)
    {
      results[rq] /= size;
    }
    return 0;
  }
  
  int SymMatPhys::phys_extend(int size, const Scalar *q, std::vector<double> &results) const
  {
    SPIN_TYPE nstate = 1ULL << size, rp;
    double dtmp, dmmt;
    // results: <node>, <node^2>
    constexpr int RESULTSIZE = 2;
    results.resize(RESULTSIZE);
    for(rp=0; rp<RESULTSIZE; ++rp)
    {
      results[rp] = 0.;
    }
    for(SPIN_TYPE rp=0ULL; rp<nstate; ++rp)
    {
      dmmt = Spin::nnode(rp, size, boundary_cond);
      dtmp = std::norm(q[rp]);
      results[0] += dmmt*dtmp;
      results[1] += dmmt*dmmt*dtmp;
    }
    return 0;
  }

  int SymMatPhys::moments(int size, const Scalar *ql, const Scalar *qr, double results[TMAT_NMOMENTS]) const
  {
    return moments(size, ql, results);
  }
  
  int SymMatPhys::phys_extend(int size, const Scalar *ql, const Scalar *qr, std::vector<double> &results) const
  {
    return phys_extend(size, ql, results);
  }

  int GenMatPhys::moments(int size, const Scalar *ql, const Scalar *qr, double results[TMAT_NMOMENTS]) const
  {
    SPIN_TYPE nstate = 1ULL << (2*size), layermask = (1ull << size)-1ull, rp;
    double dtmp, dPara, dmmt;
    double dnorm = 0.;
    int rq;
    //double rpara=0., ranti=0.;
    for(rp=0; rp<TMAT_NMOMENTS; ++rp)
    {
      results[rp] = 0.;
    }
    for(SPIN_TYPE rp=0ULL; rp<nstate; ++rp)
    {
      dmmt = dPara = Spin::spinBzexp(rp & layermask, size);
      dtmp = (ql[rp]*qr[rp]).real();
      dnorm += dtmp;
      for(rq=0; rq<TMAT_NMOMENTS; ++rq)
      {
        results[rq] += dmmt*dtmp;
        dmmt *= dPara;
      }
    }
    for(rq=0; rq<TMAT_NMOMENTS; ++rq)
    {
      results[rq] /= size*dnorm;
    }
    return 0;
  }

  int GenMatPhys::phys_extend(int size, const Scalar *ql, const Scalar *qr, std::vector<double> &results) const
  {
    // not implemented yet
    constexpr int RESULTSIZE = 0;
    results.resize(RESULTSIZE);
    return 0;
  }

  int ZMatVecProd2D::moments(int size, const Scalar *ql, const Scalar *qr, double results[TMAT_NMOMENTS]) const
  {
    SPIN_TYPE nstate = 1ULL << size, rp;
    double dtmp, dPara, dmmt;
    double dnorm = 0.;
    int rq;
    //double rpara=0., ranti=0.;
    for(rp=0; rp<TMAT_NMOMENTS; ++rp)
    {
      results[rp] = 0.;
    }
    for(SPIN_TYPE rp=0ULL; rp<nstate; ++rp)
    {
      dmmt = dPara = Spin::spinBzexp(rp, size);
      dtmp = (ql[rp]*qr[rp]).real();
      dnorm += dtmp;
      for(rq=0; rq<TMAT_NMOMENTS; ++rq)
      {
        results[rq] += dmmt*dtmp;
        dmmt *= dPara;
      }
    }
    for(rq=0; rq<TMAT_NMOMENTS; ++rq)
    {
      results[rq] /= size*dnorm;
    }
    return 0;
  }
  
};
