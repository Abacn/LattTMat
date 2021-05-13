//
//  tmat3d.cpp
//  TMat
//
//  Created by Yi Hu on 1/4/21.
//

#include <iostream>

#include "tmat3d.hpp"

namespace TMat{

/** \brief Create equivalent state list for 3D transfer matrix
 *
 */
void genX3DDuplist(int size1, int size2, int boundary_cond, int censym_flag, int ising_flag, int rotate_flag, std::vector<SPIN_TYPE> &idxlist, std::vector<unsigned short> &duplist)
{
  int size = size1*size2;
  SPIN_TYPE nstate = 1ull << size, mask = nstate-1;
  SPIN_TYPE nowstate;
  equivbits buf(size);
  idxlist.clear();
  duplist.clear();
  int reducefactor=4; // reflect layers
  //open boundary condition not implemented yet
  //if(0==boundary_cond)
  //{
  reducefactor *= size;
  //}
  if(censym_flag) reducefactor *= 2;
  if(rotate_flag) reducefactor *= 2;
  idxlist.reserve(nstate/reducefactor);
  duplist.reserve(nstate/reducefactor);
  for(nowstate = 0LL; nowstate < nstate; ++nowstate)
  {
    if(buf.getbit(nowstate)) continue;
    idxlist.push_back(nowstate);
    duplist.push_back(0);
    // resolve rotate
    SPIN_TYPE nextrstate = nowstate;
    int rot_status = 0;
    do{
      // resolve flip all spins
      SPIN_TYPE nextcstate = nextrstate;
      int cen_status = 0;
      do{
        // resolve reflectx
        SPIN_TYPE nextm1state = nextcstate;
        int m1_status = 0;
        do{
          // resolve reflecty
          SPIN_TYPE nextm2state = nextm1state;
          int m2_status = 0;
          do{
            // resolve right offsetState
            for(SPIN_TYPE nexto1state=nextm2state;
                !buf.getbit(nexto1state);
                nexto1state=Spin::rffsetState_3D(nexto1state, size1, size2, 1))
            {
              for(SPIN_TYPE nexto2state=nexto1state;
                  !buf.getbit(nexto2state);
                  nexto2state=Spin::dffsetState_3D(nexto2state, size1, size2, 1))
              {
                buf.setbit(nexto2state);
                ++duplist.back();
              }
            }
            if(0==m2_status)
            {
              m2_status = 1;
              nextm2state = Spin::reflectState_3D(nextm2state, size1, size2, 1);
            }
            else break;
          }while(!buf.getbit(nextm2state));
          if(0==m1_status)
          {
            m1_status = 1;
            nextm1state = Spin::reflectState_3D(nextm1state, size1, size2, 0);
          }
          else break;
        }while(!buf.getbit(nextm1state));
        if(0==cen_status && censym_flag)
        {
          cen_status = 1;
          nextcstate = ~nextcstate & mask;
        }
        else break;
      }while(!buf.getbit(nextcstate));
      if(0==rot_status && rotate_flag)
      {
        rot_status = 1;
        nextrstate = Spin::rotateState_3D(nextrstate, size1);
      }
      else break;
    }while(!buf.getbit(nextrstate));
  }
}

/////////////////////////////////
// 3D X-ANNNI matrix-vector multiplication class
XAMMatVecProd3D::XAMMatVecProd3D(const MParameter &mpara, const int extend_cond/*=0*/): XMatVecProd3D(mpara.size, 0, extend_cond & 1),
mpara(mpara),
nleadbit(size<10 ? 1 : 2+(int)log2(size)), // Tvec is Tv^(1/2)
censym_flag( mpara.h==0 ), ising_flag(mpara.kappa==0), rotate_flag(ising_flag && size1==size2)
{
  assert(nleadbit>=1);
  if(boundary_cond)
  {
    std::cerr << "Warning: open boundary condition not implemented yet" << std::endl;
  }
  bzwj1[0] = exp(mpara.betaJ);
  bzwj1[1] = 1. / bzwj1[0];
  bzwj2[0] = exp(-mpara.betaJ*mpara.kappa); // exp(-J*kappa)
  bzwj2[1] = 1. / bzwj2[0];
  // T(internal) = T*normfactor
  double normfactorexp;
  // normfactor is 1/ground state boltzmann weight
  if(mpara.kappa<0.5)
  {
    normfactorexp = size*(mpara.kappa-3.0);
  }
  else
  {
    normfactorexp = -size*(mpara.kappa+2.0);
  }
  normfactor = exp(mpara.betaJ*normfactorexp);
  // generate duplist and idxlist. genX2DDuplist is non-member function
  genX3DDuplist(size1, size2, boundary_cond, censym_flag, ising_flag, rotate_flag, idxlist, duplist);
  // set matrix size
  ncolsize = idxlist.size();
  genTvec(normfactorexp);
  std::cout << "col size: " << ncolsize << std::endl;
}


int XAMMatVecProd3D::moments(int size, const Scalar *q, double results[TMAT_NMOMENTS]) const
{
  SPIN_TYPE nowstate;
  double dtmp, dPara, dmmt;
  SPIN_TYPE rp, rq;
  //double rpara=0., ranti=0.;
  for(rp=0; rp<TMAT_NMOMENTS; ++rp)
  {
    results[rp] = 0.;
  }
  for(rp=0; rp<cols(); ++rp)
  {
    nowstate = idxlist[rp];
    dmmt = dPara = Spin::spinBzexp(nowstate, size);
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

int XAMMatVecProd3D::phys_extend(int size, const Scalar *q, std::vector<double> &results) const
{
  SPIN_TYPE rp;
  // results: <node>, <node^2>
  constexpr int RESULTSIZE = 2;
  results.resize(RESULTSIZE);
  for(rp=0; rp<RESULTSIZE; ++rp)
  {
    // TODO
    results[rp] = NAN;
  }
  return -1;
}

void XAMMatVecProd3D::genTvec(double normfactorexp)
{
  // H = J1 <s1 s2> + J2 [s1 s2] - hJ1 s1 = -J <s1 s2> + kappa*J [s1 s2] + h*J*s1
  double halfbetaJ = 0.5*mpara.betaJ;
  double J1 = -halfbetaJ, J2 = halfbetaJ*mpara.kappa, hJ = -halfbetaJ*mpara.h, Jnormf = halfbetaJ*normfactorexp;
  Tvec = Eigen::VectorXd(ncolsize);
#pragma omp parallel for
  for(uint64_t rp=0; rp<idxlist.size(); ++rp)
  {
    SPIN_TYPE shift1astate, shift1bstate, shift2state;
    int n1aanti, n1banti, n2anti, n0anti;
    SPIN_TYPE nowstate = idxlist[rp];
    shift1astate = Spin::rffsetState_3D(nowstate, size1, size2, 1);
    shift1bstate = Spin::dffsetState_3D(nowstate, size1, size2, 1);
    shift2state = Spin::rffsetState_3D(nowstate, size1, size2, 2);
    n1aanti = Spin::nAntiSpin(nowstate, shift1astate);
    n1banti = Spin::nAntiSpin(nowstate, shift1bstate);
    n2anti = Spin::nAntiSpin(nowstate, shift2state);
    n0anti = Spin::nAntiSpin(nowstate);
    double bzexp = 2*J1*(n1aanti+n1banti-size) + J2*(2*n2anti-size) + hJ*(2*n0anti-size) + Jnormf;
    Tvec(rp) = exp(bzexp)/sqrt((double)(duplist[rp]));
  }
}

void XAMMatVecProd3D::perform_op(const double *x_in, double *y_out) const
{
  SPIN_TYPE nloop = 1ull<<nleadbit, vecsize=nstate>>nleadbit;
  SPIN_TYPE mask=(1ull<<size)-1, maskb=vecsize-1, maskc=nloop-1, maskbinv = ~maskb & (nstate-1);
  double *y_tmp = new double[vecsize]; // memory bottleneck here: nstate*4 byte memory allocated
  const int nsubleadbit = size - nleadbit;
  // transform M->T
  if(censym_flag)
  {
    // In zero field, matrix is centrosymmetric
    // only needs to consider half of the matrix
    // nleadbit always >=1, guarantees nloop>=1
    nloop >>= 1;
  }
  for(SPIN_TYPE rdup=0; rdup<nloop; ++rdup)
  {
    SPIN_TYPE rduphigh = rdup<<nsubleadbit;
    // Tv^(1/2) * x
#pragma omp parallel for
    for(SPIN_TYPE rp=0; rp<idxlist.size(); ++rp)
    {
      double dtmp = Tvec[rp]*x_in[rp];
      SPIN_TYPE nowstate = idxlist[rp];
      // find all equalbits
      SPIN_TYPE nextrstate = nowstate;
      int rot_status = 0;
      do{
        // resolve flip all spins
        SPIN_TYPE nextcstate = nextrstate;
        int cen_status = 0;
        do{
          // resolve reflectx
          SPIN_TYPE nextm1state = nextcstate;
          int m1_status = 0;
          do{
            // resolve reflecty
            SPIN_TYPE nextm2state = nextm1state;
            int m2_status = 0;
            do{
              // resolve right offsetState
              SPIN_TYPE nexto1state=nextm2state;
              do
              {
                SPIN_TYPE nexto2state=nexto1state;
                do{
                  if((nexto2state & maskbinv) == rduphigh)
                    y_tmp[nexto2state & maskb] = dtmp;
                  nexto2state=Spin::dffsetState_3D(nexto2state, size1, size2, 1);
                }while(nexto2state!=nexto1state);
                nexto1state=Spin::rffsetState_3D(nexto1state, size1, size2, 1);
              }while(nexto1state!=nextm2state);
              if(0==m2_status)
              {
                m2_status = 1;
                nextm2state = Spin::reflectState_3D(nextm2state, size1, size2, 1);
              }
              else break;
            }while(1);
            if(0==m1_status)
            {
              m1_status = 1;
              nextm1state = Spin::reflectState_3D(nextm1state, size1, size2, 0);
            }
            else break;
          }while(1);
          if(0==cen_status && censym_flag)
          {
            cen_status = 1;
            nextcstate = ~nextcstate & mask;
          }
          else break;
        }while(1);
        if(0==rot_status && rotate_flag)
        {
          rot_status = 1;
          nextrstate = Spin::rotateState_3D(nextrstate, size1);
        }
        else break;
      }while(1);
    }
    // Th * Tv^(1/2) x
    permute_mult_TL(y_tmp);
    // transform T->M
    if(censym_flag)
    {
#pragma omp parallel for
      for(SPIN_TYPE rp=0; rp<idxlist.size(); ++rp)
      {
        SPIN_TYPE nowstate = idxlist[rp];
        double dtmp = Spin::powdi(bzwj1, Spin::spinBzexp(nowstate>>nsubleadbit, rdup, nleadbit));
        // make use of the centrosymmetricity. Reduce computation by half.
        if(rdup)
          y_out[rp] += dtmp*y_tmp[nowstate & maskb] + y_tmp[~nowstate & maskb]/dtmp;
        else
          y_out[rp] = dtmp*y_tmp[nowstate & maskb] + y_tmp[~nowstate & maskb]/dtmp;
      }
    }
    else
    {
  #pragma omp parallel for
      for(SPIN_TYPE rp=0; rp<idxlist.size(); ++rp)
      {
        SPIN_TYPE nowstate = idxlist[rp];
        double dtmp = Spin::powdi(bzwj1, Spin::spinBzexp(nowstate>>nsubleadbit, rdup, nleadbit));
        if(rdup)
          y_out[rp] += dtmp*y_tmp[nowstate & maskb];
        else
          y_out[rp] = dtmp*y_tmp[nowstate & maskb];
      }
    }
  }
#pragma omp parallel for
  for(SPIN_TYPE rp=0; rp<idxlist.size(); ++rp)
  {
    y_out[rp] *= duplist[rp]*Tvec[rp];
  }
  OpCounter::inc();
  delete [] y_tmp;
}

void XAMMatVecProd3D::permute_mult_TL(double *y_tmp) const
{
  SPIN_TYPE partialnstate = nstate >> nleadbit;
  for(int it=0; it<size-nleadbit; ++it)
  {
    SPIN_TYPE maskshift = 1<<it;
#pragma omp parallel for
    for(SPIN_TYPE rp=0; rp<partialnstate; ++rp)
    {
      if(rp & maskshift)
      {
        double dtmp = y_tmp[rp]; // 1 spin
        y_tmp[rp] = bzwj1[0]*dtmp + bzwj1[1]*y_tmp[rp ^ maskshift];
        y_tmp[rp ^ maskshift] = bzwj1[1]*dtmp + bzwj1[0]*y_tmp[rp ^ maskshift];
      }
    }
  }
}

}
