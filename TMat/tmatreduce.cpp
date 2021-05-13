//
//  tmatreduce.cpp
//  TMat
//
//  Created by Yi Hu on 8/14/20.
//

#include <iostream>
#include <algorithm>
#include <unordered_set>
#include <stdexcept>
#include <omp.h>

#include "tmatreduce.hpp"

namespace TMat{

  /// generate irrep list for X-2D models
  void genX2DDuplist(int size, int boundary_cond, int censym_flag, std::vector<SPIN_TYPE> &idxlist, std::vector<unsigned char> &duplist)
  {
    SPIN_TYPE nstate = 1ull << size, mask = nstate-1;
    SPIN_TYPE nowstate, nextstate, nextmstate;
    equivbits buf(size);
    idxlist.clear();
    duplist.clear();
    if(size>=10)
    {
      int reducefactor;
      if(0==boundary_cond)
      {
        if(censym_flag) reducefactor = size*4;
        else reducefactor = size*2;
      }
      else
      {
        if(censym_flag) reducefactor = 4;
        else reducefactor = 2;
      }
      idxlist.reserve(nstate/reducefactor);
      duplist.reserve(nstate/reducefactor);
    }
    for(nowstate = 0LL; nowstate < nstate; ++nowstate)
    {
      if(!buf.getbit(nowstate))
      {
        idxlist.push_back(nowstate);
        duplist.push_back(0);
        if(0==boundary_cond)
        {
          for(nextstate = nowstate; !buf.getbit(nextstate); nextstate = nextstate>>1 | (nextstate&1)<<(size-1))// shift one bit
          {
            ++duplist.back();
            buf.setbit(nextstate);
            if(censym_flag)
            {
              nextmstate = ~nextstate & mask;
              if(!buf.getbit(nextmstate))
              {
                ++duplist.back();
                buf.setbit(nextmstate); // flip all bits
              }
            }
          }
        
          // count backwards
          for(nextstate = Spin::bitreversal(nowstate, size);
              !buf.getbit(nextstate);
              nextstate = nextstate>>1 | (nextstate&1)<<(size-1))// shift one bit
          {
            ++duplist.back();
            buf.setbit(nextstate);
            if(censym_flag)
            {
              nextmstate = ~nextstate & mask;
              if(!buf.getbit(nextmstate))
              {
                ++duplist.back();
                buf.setbit(nextmstate); // flip all bits
              }
            }
          }
        }
        else
        {
          buf.setbit(nowstate);
          ++duplist.back();
          nextstate = Spin::bitreversal(nowstate, size);
          if(nextstate != nowstate)
          {
            ++duplist.back();
            buf.setbit(nextstate);
          }
          if(censym_flag)
          {
            nextmstate = ~nowstate & mask;
            if(!buf.getbit(nextmstate))
            {
              ++duplist.back();
              buf.setbit(nextmstate);
            }
            nextmstate = ~nextstate & mask;
            if(!buf.getbit(nextmstate))
            {
              ++duplist.back();
              buf.setbit(nextmstate);
            }
          }
        }
      }
    }
  }

  /// generate irrep list for Y-2D models
  void genY2DDuplist(int size, int boundary_cond, int censym_flag, std::vector<SPIN_TYPE> &idxlist, std::vector<SPIN_TYPE> &swapidxlist, std::vector<std::array<unsigned char, 4> > &leadprdlist, std::vector<unsigned char> &duplist)
  {
    const SPIN_TYPE nstate = 1ull << (2*size);
    const SPIN_TYPE layerstate = 1ull << size, layermask = layerstate-1, mask = nstate-1;
    const unsigned char ch_size = size;
    SPIN_TYPE nowstate, nextstate, nowistate, smalleststate, tmpstate;
    equivbits buf(size);
    idxlist.clear();
    duplist.clear();
    // approximate data sizes
    if(size>=10)
    {
      if(censym_flag)
      {
        // Cn + m + i
        idxlist.reserve(nstate/(size*4));
        duplist.reserve(nstate/(size*4));
        swapidxlist.reserve(nstate/(size*4));
        leadprdlist.reserve(layerstate/(size*4));
      }
      else
      {
        // Cn + i
        idxlist.reserve(nstate/(size*2));
        duplist.reserve(nstate/(size*2));
        swapidxlist.reserve(nstate/(size*2));
        leadprdlist.reserve(layerstate/(size*2));
      }
    }
    equivbits eqb(size*2);
    for(nowstate = 0ull; nowstate<nstate; ++nowstate)
    {
      if(!eqb.getbit(nowstate))
      {
        SPIN_TYPE leftbit = nowstate >> size, rightbit = nowstate & layermask;
        // leadprdlist update
        if(idxlist.size()==0 || (idxlist.back() >> size) != leftbit)
        {
          leadprdlist.push_back({ch_size, ch_size, ch_size, ch_size});
          int count = Spin::findcyclic(leftbit>>1 | (leftbit&1)<<(size-1), leftbit, size);
          if(count<size) leadprdlist.back()[0] = ++count;
          // count backwards
          leadprdlist.back()[1] = Spin::findcyclic(Spin::bitreversal(leftbit, size), leftbit, size);
          if(censym_flag && __builtin_popcountll(leftbit)*2 == size)
          {
            // flip all spin
            leadprdlist.back()[2] = Spin::findcyclic(~leftbit & layermask, leftbit, size);
            // flip all spin and count backwards
            leadprdlist.back()[3] = Spin::findcyclic(Spin::bitreversal(~leftbit & layermask, size), leftbit, size);
          }
        }
        swapidxlist.push_back(idxlist.size());
        idxlist.push_back(nowstate);
        duplist.push_back(0);
        smalleststate = rightbit << size | leftbit;
        for(nextstate = nowstate; !eqb.getbit(nextstate); )// shift one bit
        {
          tmpstate = rightbit << size | leftbit;
          if(tmpstate < smalleststate) smalleststate = tmpstate;
          ++duplist.back();
          eqb.setbit(nextstate);
          leftbit = leftbit>>1 | (leftbit&1)<<(size-1);
          rightbit = rightbit>>1 | (rightbit&1)<<(size-1);
          nextstate = leftbit << size | rightbit;
        }
        leftbit = Spin::bitreversal(nowstate >> size, size);
        rightbit = Spin::bitreversal(nowstate & layermask, size);
        nowistate = leftbit << size | rightbit;
        for(nextstate = nowistate; !eqb.getbit(nextstate); )// shift one bit
        {
          tmpstate = rightbit << size | leftbit;
          if(tmpstate < smalleststate) smalleststate = tmpstate;
          ++duplist.back();
          eqb.setbit(nextstate);
          leftbit = leftbit>>1 | (leftbit&1)<<(size-1);
          rightbit = rightbit>>1 | (rightbit&1)<<(size-1);
          nextstate = leftbit << size | rightbit;
        }
        // flip all spin
        if(censym_flag)
        {
          nowistate = ~nowstate & mask;
          leftbit = nowistate >> size;
          rightbit = nowistate & layermask;
          for(nextstate = nowistate; !eqb.getbit(nextstate); )// shift one bit
          {
            tmpstate = rightbit << size | leftbit;
            if(tmpstate < smalleststate) smalleststate = tmpstate;
            ++duplist.back();
            eqb.setbit(nextstate);
            leftbit = leftbit>>1 | (leftbit&1)<<(size-1);
            rightbit = rightbit>>1 | (rightbit&1)<<(size-1);
            nextstate = leftbit << size | rightbit;
          }
          leftbit = Spin::bitreversal(nowistate >> size, size);
          rightbit = Spin::bitreversal(nowistate & layermask, size);
          nowistate = leftbit << size | rightbit;
          for(nextstate = nowistate; !eqb.getbit(nextstate); )// shift one bit
          {
            tmpstate = rightbit << size | leftbit;
            if(tmpstate < smalleststate) smalleststate = tmpstate;
            ++duplist.back();
            eqb.setbit(nextstate);
            leftbit = leftbit>>1 | (leftbit&1)<<(size-1);
            rightbit = rightbit>>1 | (rightbit&1)<<(size-1);
            nextstate = leftbit << size | rightbit;
          }
        }
        if(smalleststate < nowstate)
        {
          auto pos = std::lower_bound(idxlist.begin(), idxlist.end(), smalleststate);
          if(*pos == smalleststate)
          {
            long idx = pos - idxlist.begin();
            swapidxlist.back() = idx;
            swapidxlist[idx] = idxlist.size()-1;
          }
          else
          {
            // should not happen
            throw std::logic_error("Swap layer idx not found. \nThis should not happen unless memory is unexpectedly modified");
          }
        }
      }
    }
  }

  /// generate irrep list for Z-2D models
  void genZ2DDuplist(int size, int boundary_cond, int censym_flag, std::vector<SPIN_TYPE> &idxlist, std::vector<unsigned char> &duplist)
  {
    SPIN_TYPE nstate = 1ull << size, mask = nstate-1;
    SPIN_TYPE nowstate, nextstate, nextmstate;
    equivbits buf(size);
    idxlist.clear();
    duplist.clear();
    if(size>=10)
    {
      int reducefactor;
      if(0==boundary_cond)
      {
        if(censym_flag) reducefactor = size*2;
        else reducefactor = size;
      }
      else
      {
        if(censym_flag) reducefactor = 4;
        else reducefactor = 2;
      }
      idxlist.reserve(nstate/reducefactor);
      duplist.reserve(nstate/reducefactor);
    }
    for(nowstate = 0LL; nowstate < nstate; ++nowstate)
    {
      if(!buf.getbit(nowstate))
      {
        idxlist.push_back(nowstate);
        duplist.push_back(0);
        if(0==boundary_cond)
        {
          for(nextstate = nowstate; !buf.getbit(nextstate); nextstate = nextstate>>2 | (nextstate&3)<<(size-2))// shift two bits
          {
            ++duplist.back();
            buf.setbit(nextstate);
            if(censym_flag)
            {
              nextmstate = ~nextstate & mask;
              if(!buf.getbit(nextmstate))
              {
                ++duplist.back();
                buf.setbit(nextmstate); // flip all bits
              }
            }
          }
        
          // count backwards
          for(nextstate = Spin::bitreversal( nowstate>>1 | (nowstate&1)<<(size-1) , size);
              !buf.getbit(nextstate);
              nextstate = nextstate>>2 | (nextstate&3)<<(size-2))// shift two bit
          {
            ++duplist.back();
            buf.setbit(nextstate);
            if(censym_flag)
            {
              nextmstate = ~nextstate & mask;
              if(!buf.getbit(nextmstate))
              {
                ++duplist.back();
                buf.setbit(nextmstate); // flip all bits
              }
            }
          }
        }
        else
        {
          buf.setbit(nowstate);
          ++duplist.back();
          nextstate = Spin::bitreversal( nowstate>>1 | (nowstate&1)<<(size-1) , size);
          if(nextstate != nowstate)
          {
            ++duplist.back();
            buf.setbit(nextstate);
          }
          if(censym_flag)
          {
            nextmstate = ~nowstate & mask;
            if(!buf.getbit(nextmstate))
            {
              ++duplist.back();
              buf.setbit(nextmstate);
            }
            nextmstate = ~nextstate & mask;
            if(!buf.getbit(nextmstate))
            {
              ++duplist.back();
              buf.setbit(nextmstate);
            }
          }
        }
      }
    }
  }

  /////////////////////////////////
  // 2D X-ANNNI matrix-vector multiplication class
  XAMMatVecProd2D::XAMMatVecProd2D(const MParameter &mpara, const int extend_cond/*=0*/): XMatVecProd2D(mpara.size, ((SPIN_TYPE)1) << mpara.size, extend_cond & 1),
  mpara(mpara),
  nleadbit(mpara.size<10 || boundary_cond ? 1 : 2+(int)log2(mpara.size)), // Tvec is Tv^(1/2)
  censym_flag( mpara.h==0  && (extend_cond & 2) == 0 )
  {
    assert(nleadbit>=1);
    bzwj1[0] = exp(mpara.betaJ);
    bzwj1[1] = 1. / bzwj1[0];
    bzwj2[0] = exp(-mpara.betaJ*mpara.kappa); // exp(-J*kappa)
    bzwj2[1] = 1. / bzwj2[0];
    // T(internal) = T*normfactor
    double normfactorexp;
    // normfactor is 1/ground state boltzmann weight
    if(mpara.kappa<0.5)
    {
      normfactorexp = size*(mpara.kappa-2.0);
    }
    else
    {
      normfactorexp = -size*(mpara.kappa+1.0);
    }
    normfactor = exp(mpara.betaJ*normfactorexp);
    // generate duplist and idxlist. genX2DDuplist is non-member function
    genX2DDuplist(size, boundary_cond, censym_flag, idxlist, duplist);
    genTvec(normfactorexp);
    // set matrix size
    ncolsize = idxlist.size();
  }
  
void XAMMatVecProd2D::genTvec(double normfactorexp)
  {
    SPIN_TYPE mask = nstate-1;
#ifdef DEBUG
    std::cout << "TMat size: " << idxlist.size() << std::endl;
#endif
    // idxlist.shrink_to_fit();
    // duplist.shrink_to_fit();
    // H = J1 <s1 s2> + J2 [s1 s2] - hJ1 s1 = -J <s1 s2> + kappa*J [s1 s2] + h*J*s1
    double halfbetaJ = 0.5*mpara.betaJ;
    double J1 = -halfbetaJ, J2 = halfbetaJ*mpara.kappa, hJ = -halfbetaJ*mpara.h, Jnormf = halfbetaJ*normfactorexp;
    Tvec = Eigen::VectorXd(idxlist.size());
#pragma omp parallel for
    for(uint64_t rp=0; rp<idxlist.size(); ++rp)
    {
      SPIN_TYPE shift1state, shift2state;
      int n1anti, n2anti, n0anti;
      SPIN_TYPE nowstate = idxlist[rp];
      if(!boundary_cond)
      {
        shift1state = Spin::offsetState(nowstate, size, 1);
        shift2state = Spin::offsetState(nowstate, size, 2);
        n1anti = Spin::nAntiSpin(nowstate, shift1state);
        n2anti = Spin::nAntiSpin(nowstate, shift2state);
      }
      else
      {
        // open boundary
        shift1state = nowstate >> 1;
        shift2state = nowstate >> 2;
        n1anti = Spin::nAntiSpin(nowstate & mask>>1, shift1state);
        n2anti = Spin::nAntiSpin(nowstate & mask>>2, shift2state);
      }
      n0anti = Spin::nAntiSpin(nowstate);
      double bzexp = J1*(2*n1anti-size) + J2*(2*n2anti-size) + hJ*(2*n0anti-size) + Jnormf;
      Tvec(rp) = exp(bzexp)/sqrt((double)(duplist[rp]));
    }
  }
  
  void XAMMatVecProd2D::perform_op(const double *x_in, double *y_out) const
  {
    SPIN_TYPE nloop = 1<<nleadbit, vecsize=nstate>>nleadbit;
    SPIN_TYPE maskb=vecsize-1, maskc=nloop-1, maskbinv = ~maskb & (nstate-1);
    double *y_tmp = new double[vecsize]; // memory usage O(2^L / 2^(log L) ) = O(2^L/L)
    const int nsubleadbit = size - nleadbit;
    // transform M->T
    if(censym_flag)
    {
      // In zero field, matrix is centrosymmetric
      // only needs to consider half of the matrix
      // nleadbit always >=1, guarantees nloop>=1
      nloop >>= 1;
      for(SPIN_TYPE rdup=0; rdup<nloop; ++rdup)
      {
        SPIN_TYPE rduphigh = rdup<<nsubleadbit;
        SPIN_TYPE rduprev = Spin::bitreversal(rdup, nleadbit);
        // Tv^(1/2) * x
#pragma omp parallel for
        for(SPIN_TYPE rp=0; rp<idxlist.size(); ++rp)
        {
          double dtmp = Tvec[rp]*x_in[rp];
          SPIN_TYPE nowstate = idxlist[rp];
          // shift one bit
          SPIN_TYPE nextstate = nowstate;
          if(boundary_cond)
          {
            if((nextstate & maskbinv) == rduphigh)
              y_tmp[nextstate & maskb] = dtmp;
            else if((~nextstate & maskbinv) == rduphigh)
              y_tmp[~nextstate & maskb] = dtmp;  // flip all bits
            if((nextstate & maskc) == rduprev)
              y_tmp[Spin::bitreversal(nextstate, size) & maskb] = dtmp;
            else if((~nextstate & maskc) == rduprev)
              y_tmp[Spin::bitreversal(~nextstate, size) & maskb] = dtmp;
          }
          else
          {
            do{
              if((nextstate & maskbinv) == rduphigh)
                y_tmp[nextstate & maskb] = dtmp;
              else if((~nextstate & maskbinv) == rduphigh)
                y_tmp[~nextstate & maskb] = dtmp;  // flip all bits
              if((nextstate & maskc) == rduprev)
                y_tmp[Spin::bitreversal(nextstate, size) & maskb] = dtmp;
              else if((~nextstate & maskc) == rduprev)
                y_tmp[Spin::bitreversal(~nextstate, size) & maskb] = dtmp;
              nextstate = nextstate>>1 | (nextstate&1)<<(size-1);
            }while(nextstate != nowstate);
          }
        }
        // Th * Tv^(1/2) x
        permute_mult_TL(y_tmp);
        // transform T->M
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
    }
    else
    {
      // nonzero field
      for(SPIN_TYPE rdup=0; rdup<nloop; ++rdup)
      {
        SPIN_TYPE rduphigh = rdup<<nsubleadbit;
        SPIN_TYPE rduprev = Spin::bitreversal(rdup, nleadbit);
        // Tv^(1/2) * x
#pragma omp parallel for
        for(SPIN_TYPE rp=0; rp<idxlist.size(); ++rp)
        {
          double dtmp = Tvec[rp]*x_in[rp];
          SPIN_TYPE nowstate = idxlist[rp];
          // shift one bit
          SPIN_TYPE nextstate = nowstate;
          if(boundary_cond)
          {
            if((nextstate&maskbinv) == rduphigh)
              y_tmp[nextstate & maskb] = dtmp;
            if((nextstate & maskc) == rduprev)
              y_tmp[Spin::bitreversal(nextstate, size) & maskb] = dtmp;
          }
          else
          {
            do{
              if((nextstate&maskbinv) == rduphigh)
                y_tmp[nextstate & maskb] = dtmp;
              if((nextstate & maskc) == rduprev)
                y_tmp[Spin::bitreversal(nextstate, size) & maskb] = dtmp;
              nextstate = nextstate>>1 | (nextstate&1)<<(size-1);
            }while(nextstate != nowstate);
          }
        }
        // Th * Tv^(1/2) x
        permute_mult_TL(y_tmp);
        // after the loop y_p1 stores the result
        // transform T->M
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
  
  // same as XATMatVecProd2D
  void XAMMatVecProd2D::permute_mult_T1(const double *x_in, double *y_out) const
  {
    // x_in and y_out: nstate vector
#pragma omp parallel for
    for(SPIN_TYPE rp=0; rp<nstate; ++rp)
    {
      SPIN_TYPE rq = Spin::offsetState(rp, size, 1); // permuted idx
      SPIN_TYPE rpf = rp ^ 1; // flip last bit
      y_out[rq] = bzwj1[0]*x_in[rp] + bzwj1[1]*x_in[rpf];
    }
  }
  
  // conduct half part of T1*x for (L-1) times
  // To save memory of a factor of 2
  void XAMMatVecProd2D::permute_mult_TL(double *y_tmp) const
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
  
  int XAMMatVecProd2D::moments(int size, const Scalar *q, double results[TMAT_NMOMENTS]) const
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
  
  int XAMMatVecProd2D::phys_extend(int size, const Scalar *q, std::vector<double> &results) const
  {
    SPIN_TYPE nowstate, rp;
    double dtmp, dmmt;
    // results: <node>, <node^2>
    constexpr int RESULTSIZE = 2;
    results.resize(RESULTSIZE);
    for(rp=0; rp<RESULTSIZE; ++rp)
    {
      results[rp] = 0.;
    }
    for(SPIN_TYPE rp=0ULL; rp<cols(); ++rp)
    {
      nowstate = idxlist[rp];
      dmmt = Spin::nnode(nowstate, size, boundary_cond);
      dtmp = std::norm(q[rp]);
      results[0] += dmmt*dtmp;
      results[1] += dmmt*dmmt*dtmp;
    }
    return 0;
  }

  std::vector<SPIN_TYPE> XAMMatVecProd2D::getEqvStates(SPIN_TYPE nowstate) const
  {
    std::unordered_set<SPIN_TYPE> results;
    SPIN_TYPE nextstate, nextmstate;
    SPIN_TYPE mask = nstate-1;
    if(0==boundary_cond)
    {
      for(nextstate = nowstate; results.find(nextstate)==results.end(); nextstate = nextstate>>1 | (nextstate&1)<<(size-1))// shift one bit
      {
        results.insert(nextstate);
        if(censym_flag)
        {
          nextmstate = ~nextstate & mask;
          results.insert(nextmstate);
        }
      }
      // count backwards
      for(nextstate = Spin::bitreversal(nowstate, mpara.size);
          results.find(nextstate)==results.end();
          nextstate = nextstate>>1 | (nextstate&1)<<(size-1))// shift one bit
      {
        results.insert(nextstate);
        if(censym_flag)
        {
          nextmstate = ~nextstate & mask;
          results.insert(nextmstate);
        }
      }
    }
    else
    {
      results.insert(nextstate);
      nextstate = Spin::bitreversal(nowstate, size);
      results.insert(nextstate);
      if(censym_flag)
      {
        nextmstate = ~nowstate & mask;
        results.insert(nextmstate);
      }
    }
    // construct vector from set
    return std::vector<SPIN_TYPE>(results.begin(), results.end());
  }

  /////////////////////////////////
  // 2D y-ANNNI matrix-vector multiplication class
  YAMMatVecProd2D::YAMMatVecProd2D(const MParameter &mpara,int extend_cond/*=0*/):
  YMatVecProd2D(mpara.size, ((SPIN_TYPE)1) << (2*mpara.size)),
  mpara(mpara), // Tvec is Tv
  censym_flag(mpara.h==0.0 && (extend_cond & 2) == 0),
  layertype(extend_cond&16? 1 : 0)
  {
    bzwj1[0] = exp(mpara.betaJ);
    bzwj1[1] = 1. / bzwj1[0];
    bzwj2[0] = exp(-mpara.betaJ*mpara.kappa); // exp(-J*kappa)
    bzwj2[1] = 1. / bzwj2[0];
    // for all frustrated case
    if(1 == layertype)
    {
      bzwj3[0] = exp(-mpara.betaJ*mpara.kappa2); // exp(-J*kappa)
      bzwj3[1] = 1. / bzwj3[0];
    }
    genY2DDuplist(size, boundary_cond, censym_flag, idxlist, swapidxlist, leadprdlist, duplist);
    
    // T(internal) = T*normfactor
    // normfactor is 1/ground state boltzmann weight
    double normfactorexp;
    // normfactor is 1/ground state boltzmann weight

    if(0 == layertype)
    {
      if(mpara.kappa<0.5)
      {
        normfactorexp = size*(mpara.kappa-2.0);
      }
      else
      {
        normfactorexp = -size*(mpara.kappa+1.0);
      }
      genTvec(normfactorexp);
    }
    else
    {
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
      genTvec1(normfactorexp);
    }
    normfactor = exp(mpara.betaJ*normfactorexp);
    // set matrix size
    ncolsize = idxlist.size();
  }
  
  // generate Tvec, idxlist and duplist
  void YAMMatVecProd2D::genTvec(double normfactorexp)
  {
    const SPIN_TYPE layerstate = 1ull << size, layermask = layerstate-1;
    // set Tvec
#ifdef DEBUG
    std::cout << "TMat size: " << idxlist.size() << std::endl;
    std::cout << "leadprd size: " << leadprdlist.size() << std::endl;
#endif
    // idxlist.shrink_to_fit();
    // duplist.shrink_to_fit();
    // H = J1 <s1 s2> + J2 [s1 s2] - hJ1 s1 = -J <s1 s2> + kappa*J [s1 s2] + h*J*s1
    double halfbetaJ = 0.5*mpara.betaJ;
    double J1 = -halfbetaJ, hJ = -halfbetaJ*mpara.h, Jnormf = halfbetaJ*normfactorexp;
    Tvec = Eigen::VectorXd(idxlist.size());
//#pragma omp parallel for
    for(uint64_t rp=0; rp<idxlist.size(); ++rp)
    {
      SPIN_TYPE nowstate = idxlist[rp];
      SPIN_TYPE leftbit = nowstate >> size, rightbit = nowstate & layermask;
      SPIN_TYPE shift1state = Spin::offsetState(leftbit, size, 1);
      int n1anti = Spin::nAntiSpin(leftbit, shift1state);
      int n2anti = Spin::nAntiSpin(leftbit, rightbit);
      int n0anti = Spin::nAntiSpin(leftbit);
      double bzexp = J1*(2*(n1anti+n2anti-size)) + hJ*(2*n0anti-size) + Jnormf;
      Tvec(rp) = exp(bzexp)/sqrt((double)(duplist[rp]));
    }
  }
  
  // generate Tvec, idxlist and duplist, all frustrated
  void YAMMatVecProd2D::genTvec1(double normfactorexp)
  {
    const SPIN_TYPE layerstate = 1ull << size, layermask = layerstate-1;
    // set Tvec
  #ifdef DEBUG
    std::cout << "TMat size: " << idxlist.size() << std::endl;
    std::cout << "leadprd size: " << leadprdlist.size() << std::endl;
  #endif
    // idxlist.shrink_to_fit();
    // duplist.shrink_to_fit();
    // H = J1 <s1 s2> + J2 [s1 s2] - hJ1 s1 = -J <s1 s2> + kappa*J [s1 s2] + h*J*s1
    double halfbetaJ = 0.5*mpara.betaJ;
    double J1 = -halfbetaJ, hJ = -halfbetaJ*mpara.h, Jnormf = halfbetaJ*normfactorexp;
    double kappa1J = halfbetaJ*mpara.kappa, kappa2J = halfbetaJ*mpara.kappa2;
    Tvec = Eigen::VectorXd(idxlist.size());
//#pragma omp parallel for
    for(uint64_t rp=0; rp<idxlist.size(); ++rp)
    {
      SPIN_TYPE nowstate = idxlist[rp];
      SPIN_TYPE leftbit = nowstate >> size, rightbit = nowstate & layermask;
      SPIN_TYPE shift1state = Spin::offsetState(leftbit, size, 1);
      SPIN_TYPE shift2state = Spin::offsetState(leftbit, size, 2);
      SPIN_TYPE cross1state = Spin::offsetState(rightbit, size, 1);
      SPIN_TYPE cross2state = Spin::lffsetState(rightbit, size, 1);
      int n1anti = Spin::nAntiSpin(leftbit, shift1state);
      int n2anti = Spin::nAntiSpin(leftbit, rightbit);
      int n0anti = Spin::nAntiSpin(leftbit);
      int b1anti = Spin::nAntiSpin(leftbit, shift2state);
      int d1anti = Spin::nAntiSpin(leftbit, cross1state);
      int d2anti = Spin::nAntiSpin(leftbit, cross2state);
      double bzexp = J1*(2*(n1anti+n2anti-size)) + kappa1J*(2*b1anti - size) + kappa2J*(2*(d1anti+d2anti-size)) + hJ*(2*n0anti-size) + Jnormf;
      Tvec(rp) = exp(bzexp)/sqrt((double)(duplist[rp]));
    }
  }

  // One operation. y_out = T1^T * P^T * x_in
  void YAMMatVecProd2D::permute_mult_T1(const double *x_in, double *y_out) const
  {
    // x_in and y_out: nstate vector
    SPIN_TYPE layermask = (1ull << size)-1ull;
//#pragma omp parallel for
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
  
  // x_tmp: 2^L * 2^L
  void YAMMatVecProd2D::permute_mult_TL(double *y_tmp) const
  {
    const SPIN_TYPE layersize = 1ull << size;
    for(int it=0; it<size; ++it)
    {
      SPIN_TYPE maskshift = 1<<it;
//#pragma omp parallel for
      for(SPIN_TYPE rp=0; rp<layersize; ++rp)
      {
        if(rp & maskshift)
        {
          double dtmp = y_tmp[rp]; // 1 spin
          y_tmp[rp] = bzwj2[0]*dtmp + bzwj2[1]*y_tmp[rp ^ maskshift];
          y_tmp[rp ^ maskshift] = bzwj2[1]*dtmp + bzwj2[0]*y_tmp[rp ^ maskshift];
        }
      }
    }
  }
  
  // y_out = M * x_in
  void YAMMatVecProd2D::perform_op(const double *x_in, double *y_out) const
  {
    const SPIN_TYPE layersize = 1ull << size, masklayer = layersize - 1;
    const unsigned char ch_size = size;
    SPIN_TYPE leadrp = 0, rp, rq;
    double *y_tmp = new double[layersize];
    // centre bits are in higher bits, right bits are in lower bits
    for(rp=0; rp<idxlist.size(); rp=rq, ++leadrp)
    {
      SPIN_TYPE leftbit = idxlist[rp] >> size;
      for(rq = rp; rq<idxlist.size() && idxlist[rq]>>size == leftbit; ++rq)
      {
        double dtmp = Tvec(rq)*x_in[rq];
        // expand rq
        SPIN_TYPE rightbit = idxlist[rq] & masklayer;
        for(unsigned char ch_rp=0; ch_rp<4; ++ch_rp)
        {
          unsigned char ch_shift = leadprdlist[leadrp][ch_rp];
          SPIN_TYPE rightnewbit = rightbit;
          if(!ch_rp || ch_shift != ch_size)
          {
            switch(ch_rp)
            {
              case 0: break;
              case 1: rightnewbit = Spin::bitreversal(rightbit, size); break;
              case 2: rightnewbit = ~rightbit & masklayer; break;
              case 3: rightnewbit = Spin::bitreversal(~rightbit & masklayer, size); break;
            }
            if(ch_rp)
            {
              rightnewbit = Spin::offsetState(rightnewbit, size, ch_shift);
            }
            for(unsigned char sizeA=0; sizeA<ch_size; sizeA+=leadprdlist[leadrp][0])
            {
              y_tmp[rightnewbit] = dtmp;
              rightnewbit = Spin::offsetState(rightnewbit, size, leadprdlist[leadrp][0]);
            }
          }
        }
      }
      permute_mult_TL(y_tmp);
      // send y_tmp to y_out
//#pragma omp parallel for //-- OMP does not gain speed in Y-TM
      for(SPIN_TYPE rr=rp; rr<rq; ++rr)
      {
        SPIN_TYPE swapidx = swapidxlist[rr];
        SPIN_TYPE rightbit = idxlist[rr] & masklayer;
        y_out[swapidx] = Tvec(swapidx)*duplist[swapidx]*y_tmp[rightbit];
      }
    }
    OpCounter::inc();
    delete [] y_tmp;
  }
  
  int YAMMatVecProd2D::moments(int size, const Scalar *ql, const Scalar *qr, double results[TMAT_NMOMENTS]) const
  {
    SPIN_TYPE rp;
    double dtmp, dPara, dmmt;
    double dnorm = 0.;
    int rq;
    //double rpara=0., ranti=0.;
    for(rq=0; rq<TMAT_NMOMENTS; ++rq)
    {
      results[rq] = 0.;
    }
    for(rp=0ULL; rp<idxlist.size(); ++rp)
    {
      dmmt = dPara = Spin::spinBzexp(idxlist[rp] >> size, size);
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

#define SCANCYCLE(nowstate) for(nextstate = nowstate; results.find(nextstate) == results.end(); )/*shift one bit*/\
{\
  results.insert(nextstate);\
  leftbit = leftbit>>1 | (leftbit&1)<<(size-1);\
  rightbit = rightbit>>1 | (rightbit&1)<<(size-1);\
  nextstate = leftbit << size | rightbit;\
}
  std::vector<SPIN_TYPE> YAMMatVecProd2D::getEqvStates(SPIN_TYPE nowstate) const
  {
    std::unordered_set<SPIN_TYPE> results;
    const SPIN_TYPE layerstate = 1ull << size, layermask = layerstate-1, mask = nstate-1;
    SPIN_TYPE nextstate, nowistate;
    SPIN_TYPE leftbit = nowstate >> size, rightbit = nowstate & layermask;
    SCANCYCLE(nowstate)
    leftbit = Spin::bitreversal(nowstate >> size, size);
    rightbit = Spin::bitreversal(nowstate & layermask, size);
    nowistate = leftbit << size | rightbit;
    SCANCYCLE(nowistate)
    // flip all spin
    if(censym_flag)
    {
      nowistate = ~nowstate & mask;
      leftbit = nowistate >> size;
      rightbit = nowistate & layermask;
      SCANCYCLE(nowistate)
      leftbit = Spin::bitreversal(nowistate >> size, size);
      rightbit = Spin::bitreversal(nowistate & layermask, size);
      nowistate = leftbit << size | rightbit;
      SCANCYCLE(nowistate)
    }
    return std::vector<SPIN_TYPE>(results.begin(), results.end());
  }
#undef SCANCYCLE

  uint64_t YAMMatVecProd2D::getIrrepStateIdx(uint64_t nowstate) const
  {
    auto allstate = getEqvStates(nowstate);
    std::sort(allstate.begin(), allstate.end() );
    auto it = std::lower_bound(idxlist.begin(), idxlist.end(), allstate[0] );
    uint64_t index = std::distance(idxlist.begin(), it);
    return index;
  }

// 2D y-ANNNI matrix-vector transport multiplication class
YAMMatVecProd2D_T::YAMMatVecProd2D_T(YAMMatVecProd2D const &tmathandle):
YMatVecProd2D(tmathandle.size, ((SPIN_TYPE)1) << (2*tmathandle.size)), tmathandle(tmathandle)
{
  ncolsize = tmathandle.cols();
}

/// transpose matrix-vector multiplication
void YAMMatVecProd2D_T::perform_op(const double *x_in, double *y_out) const
{
  const SPIN_TYPE layersize = 1ull << size, masklayer = layersize - 1;
  const unsigned char ch_size = size;
  SPIN_TYPE leadrp = 0, rp, rq;
  double *y_tmp = new double[layersize];
  // centre bits are in higher bits, right bits are in lower bits
  for(rp=0; rp<tmathandle.idxlist.size(); rp=rq, ++leadrp)
  {
    SPIN_TYPE leftbit = tmathandle.idxlist[rp] >> size;
    
    for(rq = rp; rq<tmathandle.idxlist.size() && tmathandle.idxlist[rq]>>size == leftbit; ++rq)
    {
      SPIN_TYPE rightbit = tmathandle.idxlist[rq] & masklayer;
      // swapidx equivalent to rightbit . leftbit
      SPIN_TYPE swapidx = tmathandle.swapidxlist[rq];
      double dtmp = tmathandle.Tvec(swapidx)*x_in[swapidx];
//#pragma omp parallel for
      for(unsigned char ch_rp=0; ch_rp<4; ++ch_rp)
      {
        unsigned char ch_shift = tmathandle.leadprdlist[leadrp][ch_rp];
        SPIN_TYPE rightnewbit = rightbit;
        if(!ch_rp || ch_shift != ch_size)
        {
          switch(ch_rp)
          {
            case 0: break;
            case 1: rightnewbit = Spin::bitreversal(rightbit, size); break;
            case 2: rightnewbit = ~rightbit & masklayer; break;
            case 3: rightnewbit = Spin::bitreversal(~rightbit & masklayer, size); break;
          }
          if(ch_rp)
          {
            rightnewbit = Spin::offsetState(rightnewbit, size, ch_shift);
          }
          for(unsigned char sizeA=0; sizeA<ch_size; sizeA+=tmathandle.leadprdlist[leadrp][0])
          {
            y_tmp[rightnewbit] = dtmp;
            rightnewbit = Spin::offsetState(rightnewbit, size, tmathandle.leadprdlist[leadrp][0]);
          }
        }
      }
    }
    tmathandle.permute_mult_TL(y_tmp);
    // send y_tmp to y_out
//#pragma omp parallel for
    for(SPIN_TYPE rr=rp; rr<rq; ++rr)
    {
      SPIN_TYPE rightbit = tmathandle.idxlist[rr] & masklayer;
      y_out[rr] = tmathandle.Tvec(rr)*tmathandle.duplist[rr]*y_tmp[rightbit];
    }
  }
  OpCounter::inc();
  delete [] y_tmp;
}

YBMMatVecProd2D::YBMMatVecProd2D(const MParameter &mpara, int extend_cond/*=0*/):
YAMMatVecProd2D(mpara, extend_cond | 16)
{
}

// 2D DNNI reduced matrix-vector multiplication class
XDMMatVecProd2D::XDMMatVecProd2D(const MParameter &mpara, int extend_cond):
XMatVecProd2D(mpara.size, ((SPIN_TYPE)1) << mpara.size, 0),
mpara(mpara),
nleadbit(mpara.size<10 ? 1 : 2+(int)log2(mpara.size)), // Tvec is Tv^(1/2)
censym_flag( mpara.h==0 && (extend_cond & 2) == 0)
{
  assert(nleadbit>=1);
  bzwj1[0] = exp(mpara.betaJ);
  bzwj1[1] = 1. / bzwj1[0];
  bzwj2[0] = exp(-mpara.betaJ*mpara.kappa); // exp(-J*kappa)
  bzwj2[1] = 1. / bzwj2[0];
  // T(internal) = T*normfactor
  double normfactorexp;
  // normfactor is 1/ground state boltzmann weight
  if(mpara.kappa<0.5)
  {
    normfactorexp = size*2.0*(mpara.kappa-1.0);
  }
  else
  {
    normfactorexp = -size*2.0*mpara.kappa;
  }
  normfactor = exp(mpara.betaJ*normfactorexp);
  // generate duplist and idxlist. genX2DDuplist is non-member function
  // boundary condition is set to 0 (periodic boundary). Open boundary not implemented yet
  genX2DDuplist(size, boundary_cond, censym_flag, idxlist, duplist);
  // set matrix size
  ncolsize = idxlist.size();
  genTvec(normfactorexp);
}

void XDMMatVecProd2D::perform_op(const double *x_in, double *y_out) const
{
  // nleadbit always >= 1
  SPIN_TYPE nloop = 1<<nleadbit;
  const SPIN_TYPE vecsize=nstate>>nleadbit;
  const SPIN_TYPE mask=nstate-1, maskb=vecsize-1, maskc=nloop-1, maskd = maskc>>1;
  const SPIN_TYPE maskbinv = ~maskb & (nstate-1);
  const int nsubleadbit = size - nleadbit;
  SPIN_TYPE nsubleadstate = 1ull << nsubleadbit;
  // memory usage O(2^L / 2^(log L) ) = O(2^L/L)
  // 2 stands for auxiliary bit
  double *x_tmp = new double[vecsize], *y_tmp = new double[2*vecsize];
  // transform M->T
#pragma omp parallel for
  for(SPIN_TYPE rp=0; rp<ncolsize; ++rp)
  {
    y_out[rp] = 0.0;
  }
  if(censym_flag)
  {
    // In zero field, matrix is centrosymmetric
    // only needs to consider half of the matrix
    // nleadbit always >=1, guarantees nloop>=1
    nloop >>= 1;
    for(SPIN_TYPE rdup=0; rdup<nloop; ++rdup)
    {
      SPIN_TYPE rduphigh = rdup<<nsubleadbit;
      SPIN_TYPE rduprev = Spin::bitreversal(rdup, nleadbit);
      // Tv^(1/2) * x
#pragma omp parallel for
      for(SPIN_TYPE rp=0; rp<ncolsize; ++rp)
      {
        double dtmp = Tvec[rp]*x_in[rp];
        SPIN_TYPE nowstate = idxlist[rp];
        // shift one bit
        SPIN_TYPE nextstate = nowstate;
        if(boundary_cond)
        {
          if((nextstate & maskbinv) == rduphigh)
            x_tmp[nextstate & maskb] = dtmp;
          else if((~nextstate & maskbinv) == rduphigh)
            x_tmp[~nextstate & maskb] = dtmp;  // flip all bits
          if((nextstate & maskc) == rduprev)
            x_tmp[Spin::bitreversal(nextstate, size) & maskb] = dtmp;
          else if((~nextstate & maskc) == rduprev)
            x_tmp[Spin::bitreversal(~nextstate, size) & maskb] = dtmp;
        }
        else
        {
          do{
            if((nextstate & maskbinv) == rduphigh)
              x_tmp[nextstate & maskb] = dtmp;
            else if((~nextstate & maskbinv) == rduphigh)
              x_tmp[~nextstate & maskb] = dtmp;  // flip all bits
            if((nextstate & maskc) == rduprev)
              x_tmp[Spin::bitreversal(nextstate, size) & maskb] = dtmp;
            else if((~nextstate & maskc) == rduprev)
              x_tmp[Spin::bitreversal(~nextstate, size) & maskb] = dtmp;
            nextstate = nextstate>>1 | (nextstate&1)<<(size-1);
          }while(nextstate != nowstate);
        }
      }
      // fbit: highest bit of new layer
      for(int fbit: {0, 1})
      {
#pragma omp parallel for
        for(SPIN_TYPE rp=0; rp<vecsize; ++rp)
        {
          // in centro-symmetric case, highest bit of old layer is always 0
          y_tmp[rp] = x_tmp[rp]*bzwj2[(rp&1) ^ fbit];
          y_tmp[rp | nsubleadstate ] = 0.0;
        }
        // Th * Tv^(1/2) x
        permute_mult_TL(y_tmp, (rdup & 1) );
        // transform T->M
#pragma omp parallel for
        for(SPIN_TYPE rp=0; rp<idxlist.size(); ++rp)
        {
          SPIN_TYPE nowstate = idxlist[rp];
          if(fbit)
          {
            // flip all spins
            nowstate = ~nowstate & mask;
          }
          // calculate leading bits interaction
          SPIN_TYPE rdupb = nowstate >> nsubleadbit;
          // Ising interaction
          double dtmp = Spin::powdi(bzwj1, Spin::spinBzexp(rdup, rdupb, nleadbit));
          // diagonal frustration
          dtmp *= Spin::powdi(bzwj2, Spin::spinBzexp((rdup>>1 ^ rdupb) & maskd, nleadbit-1) + Spin::spinBzexp((rdup ^ rdupb>>1) & maskd, nleadbit-1));
          // leading/subl boundary diagonal interaction /
          y_out[rp] += dtmp*(y_tmp[nowstate & maskb]*bzwj2[rdupb&1] + y_tmp[(nowstate & maskb) | vecsize]*bzwj2[(rdupb&1) ^ 1]);
        }
      }
    }
  }
  else
  {
    // not imposing centrosymmetric
    for(SPIN_TYPE rdup=0; rdup<nloop; ++rdup)
    {
      SPIN_TYPE rduphigh = rdup<<nsubleadbit;
      SPIN_TYPE rduprev = Spin::bitreversal(rdup, nleadbit);
      // Tv^(1/2) * x
#pragma omp parallel for
      for(SPIN_TYPE rp=0; rp<idxlist.size(); ++rp)
      {
        double dtmp = Tvec[rp]*x_in[rp];
        SPIN_TYPE nowstate = idxlist[rp];
        // shift one bit
        SPIN_TYPE nextstate = nowstate;
        if(boundary_cond)
        {
          if((nextstate&maskbinv) == rduphigh)
            x_tmp[nextstate & maskb] = dtmp;
          if((nextstate & maskc) == rduprev)
            x_tmp[Spin::bitreversal(nextstate, size) & maskb] = dtmp;
        }
        else
        {
          do{
            if((nextstate&maskbinv) == rduphigh)
              x_tmp[nextstate & maskb] = dtmp;
            if((nextstate & maskc) == rduprev)
              x_tmp[Spin::bitreversal(nextstate, size) & maskb] = dtmp;
            nextstate = nextstate>>1 | (nextstate&1)<<(size-1);
          }while(nextstate != nowstate);
        }
      }
      // fbit: highest bit of new layer
      for(int fbit: {0, 1})
      {
        // gbit: highest bit of old layer
        SPIN_TYPE gbit = rdup >> (nleadbit-1) << nsubleadbit;
#pragma omp parallel for
        for(SPIN_TYPE rp=0; rp<vecsize; ++rp)
        {
          // in centro-symmetric case, highest bit of old layer is always 0
          y_tmp[rp | gbit] = x_tmp[rp]*bzwj2[(rp&1) ^ fbit];
          y_tmp[rp | (gbit ^ nsubleadstate) ] = 0.0;
        }
        // Th * Tv^(1/2) x
        permute_mult_TL(y_tmp, (rdup & 1) );
        // transform T->M
#pragma omp parallel for
        for(SPIN_TYPE rp=0; rp<idxlist.size(); ++rp)
        {
          SPIN_TYPE nowstate = idxlist[rp];
          // highest bit match
          if(!(nowstate >> (size-1) ^ fbit))
          {
            // calculate leading bits interaction
            SPIN_TYPE rdupb = nowstate >> nsubleadbit;
            // Ising interaction
            double dtmp = Spin::powdi(bzwj1, Spin::spinBzexp(rdup, rdupb, nleadbit));
            // diagonal frustration
            dtmp *= Spin::powdi(bzwj2, Spin::spinBzexp((rdup>>1 ^ rdupb) & maskd, nleadbit-1) + Spin::spinBzexp((rdup ^ rdupb>>1) & maskd, nleadbit-1));
            // leading/subl boundary diagonal interaction /
            y_out[rp] += dtmp*(y_tmp[nowstate & maskb]*bzwj2[rdupb&1] + y_tmp[(nowstate & maskb) | vecsize]*bzwj2[(rdupb&1) ^ 1]);
          }

        }
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
  delete [] x_tmp;
}

void XDMMatVecProd2D::permute_mult_TL(double *y_tmp, int outbit) const
{
  // ytmp has an auxiliary bit
  SPIN_TYPE partialnstate = nstate >> nleadbit;
  const int nloop = size-nleadbit;
  for(int it=0; it<nloop; ++it)
  {
    SPIN_TYPE maskshift = 1<<it, bshift = ~maskshift & (nstate-1);
#pragma omp parallel for
    for(SPIN_TYPE rp=0; rp<partialnstate; ++rp)
    {
      if(rp & maskshift)
      {
        int nextbit = outbit;
        if(it < nloop-1)
          nextbit = rp>>(it+1) & 1;
        SPIN_TYPE i00 = rp & bshift, i01 = rp, i10 = (rp & bshift) | partialnstate, i11 = rp | partialnstate;
        double dtmp_00 = y_tmp[i00], dtmp_01 = y_tmp[i01], dtmp_10 = y_tmp[i10], dtmp_11 = y_tmp[i11];
        // old 0, new 0
        y_tmp[i00] = bzwj1[0]*bzwj2[nextbit]*(bzwj2[0]*dtmp_00 + bzwj2[1]*dtmp_10);
        // old 0, new 1
        y_tmp[i01] = bzwj1[1]*bzwj2[nextbit ^ 1]*(bzwj2[1]*dtmp_00 + bzwj2[0]*dtmp_10);
        // old 1, new 0
        y_tmp[i10] = bzwj1[1]*bzwj2[nextbit]*(bzwj2[0]*dtmp_01 + bzwj2[1]*dtmp_11);
        // old 1, new 1
        y_tmp[i11] = bzwj1[0]*bzwj2[nextbit ^ 1]*(bzwj2[1]*dtmp_01 + bzwj2[0]*dtmp_11);
      }
    }
  }
}

int XDMMatVecProd2D::moments(int size, const Scalar *q, double results[TMAT_NMOMENTS]) const
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

int XDMMatVecProd2D::phys_extend(int size, const Scalar *q, std::vector<double> &results) const
{
  const SPIN_TYPE nstate = 1ULL << size, mask = nstate-1, stripflip = 0xaaaaaaaaaaaaaaaaULL;
  SPIN_TYPE rp;
  double dtmp, dmmt, dmmtb, dbgd;
  // results: <node>, <node^2>
  constexpr int RESULTSIZE = 4;
  results.resize(RESULTSIZE);
  for(rp=0; rp<RESULTSIZE; ++rp)
  {
    results[rp] = 0.;
  }
  dbgd = 0;
  for(rp=0ULL; rp<cols(); ++rp)
  {
    // m
    SPIN_TYPE nowstate = idxlist[rp];
    dmmt = Spin::spinBzexp((nowstate ^ stripflip) & mask, size);
    dmmtb = Spin::spinBzexp(nowstate, size);
    dtmp = std::norm(q[rp]);
    dbgd += dtmp;
    results[0] += abs(dmmt)*dtmp;
    results[1] += abs(dmmtb)*dtmp;
    results[2] += dmmt*dmmt*dtmp;
    results[3] += dmmtb*dmmtb*dtmp;
  }
  results[0] /= size;
  results[1] /= size;
  results[2] /= size*size;
  results[3] /= size*size;
  return 0;
}

std::vector<SPIN_TYPE> XDMMatVecProd2D::getEqvStates(SPIN_TYPE nowstate) const
{
  // Same as 2DXAM
  std::unordered_set<SPIN_TYPE> results;
  SPIN_TYPE nextstate, nextmstate;
  SPIN_TYPE mask = nstate-1;
  if(0==boundary_cond)
  {
    for(nextstate = nowstate; results.find(nextstate)==results.end(); nextstate = nextstate>>1 | (nextstate&1)<<(size-1))// shift one bit
    {
      results.insert(nextstate);
      if(censym_flag)
      {
        nextmstate = ~nextstate & mask;
        results.insert(nextmstate);
      }
    }
    // count backwards
    for(nextstate = Spin::bitreversal(nowstate, mpara.size);
        results.find(nextstate)==results.end();
        nextstate = nextstate>>1 | (nextstate&1)<<(size-1))// shift one bit
    {
      results.insert(nextstate);
      if(censym_flag)
      {
        nextmstate = ~nextstate & mask;
        results.insert(nextmstate);
      }
    }
  }
  else
  {
    results.insert(nextstate);
    nextstate = Spin::bitreversal(nowstate, size);
    results.insert(nextstate);
    if(censym_flag)
    {
      nextmstate = ~nowstate & mask;
      results.insert(nextmstate);
    }
  }
  // construct vector from set
  return std::vector<SPIN_TYPE>(results.begin(), results.end());
}

void XDMMatVecProd2D::genTvec(double normfactorexp)
{
  SPIN_TYPE mask = nstate-1;
  SPIN_TYPE nowstate;
#ifdef DEBUG
  std::cout << "TMat size: " << idxlist.size() << std::endl;
#endif
  // idxlist.shrink_to_fit();
  // duplist.shrink_to_fit();
  // H = J1 <s1 s2> + J2 [s1 s2] - hJ1 s1 = -J <s1 s2> + kappa*J [s1 s2] + h*J*s1
  double halfbetaJ = 0.5*mpara.betaJ;
  double J1 = -halfbetaJ, hJ = -halfbetaJ*mpara.h, Jnormf = halfbetaJ*normfactorexp;
  SPIN_TYPE shift1state;
  int n1anti, n0anti;
  Tvec = Eigen::VectorXd(ncolsize);
  for(uint64_t rp=0; rp<ncolsize; ++rp)
  {
    nowstate = idxlist[rp];
    if(!boundary_cond)
    {
      shift1state = Spin::offsetState(nowstate, size, 1);
      n1anti = Spin::nAntiSpin(nowstate, shift1state);
    }
    else
    {
      // open boundary
      shift1state = nowstate >> 1;
      n1anti = Spin::nAntiSpin(nowstate & mask>>1, shift1state);
    }
    n0anti = Spin::nAntiSpin(nowstate);
    double bzexp = J1*(2*n1anti-size) + hJ*(2*n0anti-size) + Jnormf;
    Tvec(rp) = exp(bzexp)/sqrt((double)(duplist[rp]));
  }
}

// BNNNI and 3NN model propagated in the diagonal direction
ZBMMatVecProd2D::ZBMMatVecProd2D(const MParameter &mpara, int extend_cond):
ZMatVecProd2D(mpara.size, 1ull<<mpara.size), mpara(mpara),
nleadbit(mpara.size<10 ? 2 : 2+((int)log2(mpara.size) & ~1)), // Tvec is Tv^(1/2)
censym_flag( mpara.h==0 && (extend_cond & 2) == 0)
{
  assert(nleadbit>=2);
  bzwj1[0] = exp(mpara.betaJ);
  bzwj1[1] = 1. / bzwj1[0];
  bzwj2[0] = exp(-mpara.betaJ*mpara.kappa); // exp(-J*kappa1)
  bzwj2[1] = 1. / bzwj2[0];
  bzwj3[0] = exp(-mpara.betaJ*mpara.kappa2); // exp(-J*kappa2)
  bzwj3[1] = 1. / bzwj3[0];
  bzwj2b[0] = exp(-2*mpara.betaJ*mpara.kappa); // exp(-2 J kappa1);
  bzwj2b[3] = 1. / bzwj2b[0];
  bzwj2b[1] = bzwj2b[2] = 1.;
  // oneside kappa and diagonal kappa2
  for(SPIN_TYPE rp=0; rp<16; ++rp)
  {
    SPIN_TYPE ra = rp >> 2, rb = rp & 3;
    bzwj2c[rp] = exp( mpara.betaJ*(mpara.kappa*Spin::spinBzexp(ra, 2) + mpara.kappa2*Spin::spinBzexp(rb, 2)) );
  }
  // T(internal) = T*normfactor
  double normfactorexp;
  // normfactor is 1/ground state boltzmann weight
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
  normfactor = exp(mpara.betaJ*normfactorexp);
  // generate duplist and idxlist. genX2DDuplist is non-member function
  // boundary condition is set to 0 (periodic boundary). Open boundary not implemented yet
  genZ2DDuplist(size, boundary_cond, censym_flag, idxlist, duplist);
  // set matrix size
  ncolsize = idxlist.size();
  genTvec(normfactorexp);
}

void ZBMMatVecProd2D::perform_op(const double *x_in, double *y_out) const
{
  perform_op_internal(x_in, y_out, 0);
}

void ZBMMatVecProd2D::perform_op_internal(const double *x_in, double *y_out, int optdirection) const
{
  // nleadbit always >= 1
  SPIN_TYPE nloop = 1ull<<nleadbit;
  const SPIN_TYPE vecsize=nstate >> nleadbit;
  const SPIN_TYPE mask=nstate-1, maskb=vecsize-1, maskc=nloop-1;
  const SPIN_TYPE maskbinv = ~maskb & (nstate-1);
  const int nsubleadbit = size - nleadbit;
  // memory usage O(2^L / 2^(log L) ) = O(2^L/L)
  // 4 stands for 2 auxiliary bits
  double *x_tmp = new double[vecsize], *y_p1 = new double[4*vecsize], *y_p2 = new double[4*vecsize], *y_swap;
  // transform M->T
#pragma omp parallel for
  for(SPIN_TYPE rp=0; rp<ncolsize; ++rp)
  {
    y_out[rp] = 0.0;
  }
  if(censym_flag)
  {
    // In zero field, matrix is centrosymmetric
    // only needs to consider half of the matrix
    // nleadbit always >=2, guarantees nloop>=1
    nloop >>= 1;
    for(SPIN_TYPE rdup=0; rdup<nloop; ++rdup)
    {
      SPIN_TYPE rduphigh = rdup<<nsubleadbit;
      SPIN_TYPE rduphinv = (~rdup&maskc)<<nsubleadbit;
      // Tv^(1/2) * x
#pragma omp parallel for
      for(SPIN_TYPE rp=0; rp<ncolsize; ++rp)
      {
        double dtmp = Tvec[rp]*x_in[rp];
        SPIN_TYPE nowstate = idxlist[rp];
        // shift one bit
        SPIN_TYPE nextstate = nowstate;
        do{
          if((nextstate & maskbinv) == rduphigh)
            x_tmp[nextstate & maskb] = dtmp;
          else if((nextstate & maskbinv) == rduphinv)
            x_tmp[~nextstate & maskb] = dtmp;  // flip all bits
          SPIN_TYPE nextistate = Spin::bitreversal( nextstate>>1 | (nextstate&1)<<(size-1) , size);
          if((nextistate & maskbinv) == rduphigh)
            x_tmp[nextistate & maskb] = dtmp;
          else if((nextistate & maskbinv) == rduphinv)
            x_tmp[~nextistate & maskb] = dtmp;
          nextstate = nextstate>>2 | (nextstate&3)<<(size-2);
        }while(nextstate != nowstate);
      }
      // first two bits
      // for centrosym. fbit is either 00 or 01
      SPIN_TYPE fbit = rdup >> (nleadbit-2);
      for(SPIN_TYPE lastbit: {0, 1, 2, 3})
      {
#pragma omp parallel for
        for(SPIN_TYPE rp=0; rp<vecsize; ++rp)
        {
          // set x with auxiliary bits
          SPIN_TYPE rrp = rp<<2 | fbit;
          if((rp&3)==lastbit)
            y_p1[rrp] = x_tmp[rp];
          else
            y_p1[rrp] = 0.0;
          y_p1[rrp ^ 1] = y_p1[rrp ^ 2] = y_p1[rrp ^ 3] = 0.0;
        }
        // Th * Tv^(1/2) x
        for(int rp=0; rp<nsubleadbit/2-1; ++rp)
        {
          permute_mult_T1(y_p1, y_p2, optdirection);
          y_swap = y_p1;
          y_p1 = y_p2;
          y_p2 = y_swap;
        }
        // last trun
        permute_mult_T1_last(y_p1, y_p2, rdup&3, optdirection);
        // transform T->M
#pragma omp parallel for
        for(SPIN_TYPE rp=0; rp<idxlist.size(); ++rp)
        {
          SPIN_TYPE nowstate = idxlist[rp];
          double dtmp = leadingbz(y_p2, nowstate, rdup, lastbit, optdirection);
          SPIN_TYPE invstate = ~nowstate & mask;
          dtmp += leadingbz(y_p2, invstate, rdup, lastbit, optdirection);
          y_out[rp] += dtmp;
        }
      }
    }
  }
  else
  {
    // flip all spin symmetry is not considered
    for(SPIN_TYPE rdup=0; rdup<nloop; ++rdup)
    {
      SPIN_TYPE rduphigh = rdup<<nsubleadbit;
      // Tv^(1/2) * x
#pragma omp parallel for
      for(SPIN_TYPE rp=0; rp<ncolsize; ++rp)
      {
        double dtmp = Tvec[rp]*x_in[rp];
        SPIN_TYPE nowstate = idxlist[rp];
        // shift one bit
        SPIN_TYPE nextstate = nowstate;
        do{
          if((nextstate & maskbinv) == rduphigh)
            x_tmp[nextstate & maskb] = dtmp;
          SPIN_TYPE nextistate = Spin::bitreversal( nextstate>>1 | (nextstate&1)<<(size-1) , size);
          if((nextistate & maskbinv) == rduphigh)
            x_tmp[nextistate & maskb] = dtmp;
          nextstate = nextstate>>2 | (nextstate&3)<<(size-2);
        }while(nextstate != nowstate);
      }
      // first two bits
      // for centrosym. fbit is either 00 or 01
      SPIN_TYPE fbit = rdup >> (nleadbit-2);
      for(SPIN_TYPE lastbit: {0, 1, 2, 3})
      {
#pragma omp parallel for
        for(SPIN_TYPE rp=0; rp<vecsize; ++rp)
        {
          // set x with auxiliary bits
          SPIN_TYPE rrp = rp<<2 | fbit;
          if((rp&3)==lastbit)
            y_p1[rrp] = x_tmp[rp];
          else
            y_p1[rrp] = 0.0;
          y_p1[rrp ^ 1] = y_p1[rrp ^ 2] = y_p1[rrp ^ 3] = 0.0;
        }
        // Th * Tv^(1/2) x
        for(int rp=0; rp<nsubleadbit/2-1; ++rp)
        {
          permute_mult_T1(y_p1, y_p2, optdirection);
          y_swap = y_p1;
          y_p1 = y_p2;
          y_p2 = y_swap;
        }
        // last trun
        permute_mult_T1_last(y_p1, y_p2, rdup&3, optdirection);
        // transform T->M
#pragma omp parallel for
        for(SPIN_TYPE rp=0; rp<idxlist.size(); ++rp)
        {
          SPIN_TYPE nowstate = idxlist[rp];
          double dtmp = leadingbz(y_p2, nowstate, rdup, lastbit, optdirection);
          y_out[rp] += dtmp;
        }
      }
    }
  }
#pragma omp parallel for
  for(SPIN_TYPE rp=0; rp<idxlist.size(); ++rp)
  {
    y_out[rp] *= duplist[rp]*Tvec[rp];
  }
  OpCounter::inc();
  delete [] y_p1;
  delete [] y_p2;
  delete [] x_tmp;
}


void ZBMMatVecProd2D::permute_mult_T1(double *x_in, double *y_out, int optdirection) const
{
  const int nsubleadbit = size - nleadbit;
  const SPIN_TYPE vecsize=nstate >> nleadbit;
  SPIN_TYPE nstateb = vecsize<<2;
  SPIN_TYPE maskb = vecsize - 4; // 11..1, 11..00
  // x_in and y_out: nstate vector
#pragma omp parallel for
  for(SPIN_TYPE rpb=0; rpb<nstateb; ++rpb) // row
  {
    SPIN_TYPE auxp = rpb & 3;
    SPIN_TYPE rp = rpb >> 2;
    SPIN_TYPE rpp = Spin::offsetState(rp, nsubleadbit, 2)<<2 | auxp; // permuted idx
    SPIN_TYPE rq = ((rp&maskb) | auxp) << 2;
    unsigned char sbit = rp & 3, j2bidx = ((rp & 0xc) | auxp) ^ (sbit<<2 | sbit); // four interaction quick index for bzwj2b[]
    double dtmp;
    if(!optdirection)
    {
      // current right bit, previous left bit
      unsigned char srbit = rp & 1, plbit = auxp>>1;
      // get interaction energy not involved back bits
      y_out[rpp] = bzwj1[srbit^plbit] * bzwj2c[j2bidx];
      // then combine with energy involved back bits
      dtmp = x_in[rq]*bzwj1[srbit]*bzwj2b[sbit]; // 00
      dtmp += x_in[rq+1]*bzwj1[srbit]*bzwj2b[sbit^1]; // 01
      dtmp += x_in[rq+2]*bzwj1[srbit^1]*bzwj2b[sbit^2]; // 10
      dtmp += x_in[rq+3]*bzwj1[srbit^1]*bzwj2b[sbit^3]; // 11
    }
    else
    {
      // current left bit, forward right bit, previous right bit
      unsigned char slbit = rp>>1 & 1, frbit = rp>>2 & 1, prbit = auxp & 1;
      // get interaction energy not involved back bits
      y_out[rpp] = bzwj1[slbit^frbit] * bzwj1[slbit^prbit] * bzwj2c[j2bidx];
      // then combine with energy involved back bits
      dtmp = x_in[rq]*bzwj2b[sbit];  // 00
      dtmp += x_in[rq+1]*bzwj2b[sbit^1];  // 01
      dtmp += x_in[rq+2]*bzwj2b[sbit^2];  // 10
      dtmp += x_in[rq+3]*bzwj2b[sbit^3];  // 11
    }
    y_out[rpp] *= dtmp; // J term
  }
}

void ZBMMatVecProd2D::permute_mult_T1_last(double *x_in, double *y_out, SPIN_TYPE lastbit, int optdirection) const
{
  const int nsubleadbit = size - nleadbit;
  const SPIN_TYPE vecsize=nstate >> nleadbit;
  SPIN_TYPE nstateb = vecsize<<2;
  SPIN_TYPE maskb = vecsize - 4; // 11..1, 11..00
  // x_in and y_out: nstate vector
#pragma omp parallel for
  for(SPIN_TYPE rpb=0; rpb<nstateb; ++rpb) // row
  {
    SPIN_TYPE auxp = rpb & 3;
    SPIN_TYPE rp = rpb >> 2;
    SPIN_TYPE rpp = Spin::offsetState(rp, nsubleadbit, 2)<<2 | auxp; // permuted idx
    SPIN_TYPE rq = ((rp&maskb) | auxp) << 2;
    unsigned char sbit = rp & 3, j2bidx = ((lastbit << 2) | auxp) ^ (sbit<<2 | sbit); // four interaction quick index for bzwj2b[]
    double dtmp;
    if(!optdirection)
    {
      // current right bit, previous left bit
      unsigned char srbit = rp & 1, plbit = auxp>>1;
      // get interaction energy not involved back bits
      y_out[rpp] = bzwj1[srbit^plbit] * bzwj2c[j2bidx];
      // then combine with energy involved back bits
      dtmp = x_in[rq]*bzwj1[srbit]*bzwj2b[sbit]; // 00
      dtmp += x_in[rq+1]*bzwj1[srbit]*bzwj2b[sbit^1]; // 01
      dtmp += x_in[rq+2]*bzwj1[srbit^1]*bzwj2b[sbit^2]; // 10
      dtmp += x_in[rq+3]*bzwj1[srbit^1]*bzwj2b[sbit^3]; // 11
    }
    else
    {
      // current left bit, forward right bit, previous right bit
      unsigned char slbit = rp>>1 & 1, frbit = lastbit & 1, prbit = auxp & 1;
      // get interaction energy not involved back bits
      y_out[rpp] = bzwj1[slbit^frbit] * bzwj1[slbit^prbit] * bzwj2c[j2bidx];
      // then combine with energy involved back bits
      dtmp = x_in[rq]*bzwj2b[sbit];  // 00
      dtmp += x_in[rq+1]*bzwj2b[sbit^1];  // 01
      dtmp += x_in[rq+2]*bzwj2b[sbit^2];  // 10
      dtmp += x_in[rq+3]*bzwj2b[sbit^3];  // 11
    }
    y_out[rpp] *= dtmp; // J term
  }
}

double ZBMMatVecProd2D::leadingbz(double *y_tmp, SPIN_TYPE nowstate, SPIN_TYPE leadbits, SPIN_TYPE lastbit, int optdirection) const
{
  static const SPIN_TYPE MASKEVEN = 0x5555555555555555ULL; // even bits (start from 0)
  const int nsubleadbit = size-nleadbit;
  const SPIN_TYPE vecsize=nstate >> nleadbit;
  const SPIN_TYPE maskb=vecsize-1, maskc=(1ull << nleadbit)-1, maskd=maskc>>2;
  SPIN_TYPE nowleadbits = nowstate >> nsubleadbit, nowshead = (nowstate & maskb)<<2;
  // interaction of nowstate's leadbits with (leadbits(2))(leadbits(nleadbit))(lastbit2(2bits))
  // first calculate interaction that newlastbit not involved
  double df1 = Spin::powdi(bzwj3, Spin::spinBzexp(nowleadbits, leadbits, nleadbit)); // for kappa2
  df1 *= Spin::powdi(bzwj2, Spin::spinBzexp(nowleadbits, leadbits>>2 | lastbit<<(nleadbit-2), nleadbit)
                     + Spin::spinBzexp(nowleadbits>>2, leadbits & maskd, nleadbit-2) ); // for kappa diagonal one side
  SPIN_TYPE isingnow, isingnxt; // ising interaction bits
  int isingint;
  double df2;
  if(!optdirection)
  {
    isingnow = nowleadbits & MASKEVEN;
    isingnow |= isingnow >> 1;
    isingnxt = (leadbits >> 1) & MASKEVEN;
    isingnxt |= (isingnxt << 1) & (maskc>>1);
    isingint = Spin::spinBzexp(isingnow, isingnxt, nleadbit-1);
    unsigned char lastldbit = nowleadbits&1;
    df2 = bzwj1[lastldbit] * (bzwj2b[(nowleadbits&3)] * y_tmp[nowshead] + bzwj2b[(nowleadbits&3) ^ 1] * y_tmp[nowshead+1])
         + bzwj1[lastldbit^1] * (bzwj2b[(nowleadbits&3) ^ 2] * y_tmp[nowshead+2] + bzwj2b[(nowleadbits&3) ^ 3] * y_tmp[nowshead+3]);
  }
  else
  {
    isingnow = (nowleadbits >> 1) & MASKEVEN;
    isingnow |= isingnow << 1;
    isingnxt = (leadbits | (lastbit << nleadbit) ) & MASKEVEN;
    isingnxt = (isingnxt & maskc ) | (isingnxt >> 1);
    isingint = Spin::spinBzexp(isingnow, isingnxt, nleadbit);
    df2 = bzwj2b[(nowleadbits&3)] * y_tmp[nowshead] + bzwj2b[(nowleadbits&3) ^ 1] * y_tmp[nowshead+1]
         + bzwj2b[(nowleadbits&3) ^ 2] * y_tmp[nowshead+2] + bzwj2b[(nowleadbits&3) ^ 3] * y_tmp[nowshead+3];
  }
  df1 *= Spin::powdi(bzwj1, isingint);
  return df1 * df2;
}

int ZBMMatVecProd2D::moments(int size, const Scalar *ql, const Scalar *qr, double results[TMAT_NMOMENTS]) const
{
  SPIN_TYPE rp;
  double dtmp, dPara, dmmt;
  double dnorm = 0.;
  int rq;
  //double rpara=0., ranti=0.;
  for(rq=0; rq<TMAT_NMOMENTS; ++rq)
  {
    results[rq] = 0.;
  }
  for(rp=0ULL; rp<idxlist.size(); ++rp)
  {
    dmmt = dPara = Spin::spinBzexp(idxlist[rp], size);
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
  return -1; // changed to 0 after YAMMatVecProd2D_T finished
}

int ZBMMatVecProd2D::phys_extend(int size, const Scalar *ql, const Scalar *qr, std::vector<double> &results) const
{
  // TODO
  return -1;
}

std::vector<SPIN_TYPE> ZBMMatVecProd2D::getEqvStates(SPIN_TYPE nowstate) const
{
  std::unordered_set<SPIN_TYPE> results;
  SPIN_TYPE nextstate, nextmstate;
  SPIN_TYPE mask = nstate-1;
  for(nextstate = nowstate; results.find(nextstate)==results.end(); nextstate = nextstate>>2 | (nextstate&3)<<(size-2))// shift two bits
  {
    results.insert(nextstate);
    if(censym_flag)
    {
      nextmstate = ~nextstate & mask;
      results.insert(nextmstate);
    }
  }
  // count backwards
  for(nextstate = Spin::bitreversal( nowstate>>1 | (nowstate&1)<<(size-1) , size);
      results.find(nextstate)==results.end();
      nextstate = nextstate>>2 | (nextstate&3)<<(size-2))// shift one bit
  {
    results.insert(nextstate);
    if(censym_flag)
    {
      nextmstate = ~nextstate & mask;
      results.insert(nextmstate);
    }
  }
  // construct vector from set
  return std::vector<SPIN_TYPE>(results.begin(), results.end());
}

void ZBMMatVecProd2D::genTvec(double normfactorexp)
{
  SPIN_TYPE mask = nstate-1;
#ifdef DEBUG
  std::cout << "TMat size: " << idxlist.size() << std::endl;
#endif
  // H = J1 <s1 s2> + J2 [s1 s2] - hJ1 s1 = -J <s1 s2> + kappa*J [s1 s2] + h*J*s1
  double halfbetaJ = 0.5*mpara.betaJ;
  double J1 = -halfbetaJ, J3 = halfbetaJ*mpara.kappa2, hJ = -halfbetaJ*mpara.h, Jnormf = halfbetaJ*normfactorexp;
  Tvec = Eigen::VectorXd(ncolsize);
#pragma omp parallel for
  for(uint64_t rp=0; rp<ncolsize; ++rp)
  {
    int n1anti, n2anti, n0anti;
    double bz1, bz2, bz0;
    SPIN_TYPE shift1state, shift2state;
    SPIN_TYPE nowstate = idxlist[rp];
    if(!boundary_cond)
    {
      shift1state = Spin::offsetState(nowstate, mpara.size, 1);
      shift2state = Spin::offsetState(nowstate, mpara.size, 2);
      n1anti = Spin::nAntiSpin(nowstate, shift1state);
      n2anti = Spin::nAntiSpin(nowstate, shift2state);
      n0anti = Spin::nAntiSpin(nowstate);
      bz1 = Spin::spinBzExpn(mpara.size, n1anti, J1);
      bz2 = Spin::spinBzExpn(mpara.size, n2anti, J3);
      bz0 = Spin::spinBzExpn(mpara.size, n0anti, hJ);
    }
    else
    {
      // open boundary
      shift1state = nowstate >> 1;
      shift2state = nowstate >> 2;
      n1anti = Spin::nAntiSpin(nowstate & mask>>1, shift1state);
      n2anti = Spin::nAntiSpin(nowstate & mask>>2, shift2state);
      n0anti = Spin::nAntiSpin(nowstate);
      bz1 = Spin::spinBzExpn(mpara.size-1, n1anti, J1);
      bz2 = Spin::spinBzExpn(mpara.size-2, n2anti, J3);
      bz0 = Spin::spinBzExpn(mpara.size, n0anti, hJ);
    }
    Tvec(rp) = exp(bz1+bz2+bz0+Jnormf)/sqrt((double)(duplist[rp]));
  }
}

ZBMMatVecProd2D_T::ZBMMatVecProd2D_T(const ZBMMatVecProd2D &tmathandle):
ZMatVecProd2D(tmathandle.size, 1ull << tmathandle.size ),
tmathandle(tmathandle)
{
  ncolsize = tmathandle.cols();
}

// y_out = M * x_in
void ZBMMatVecProd2D_T::perform_op(const double *x_in, double *y_out) const
{
  tmathandle.perform_op_internal(x_in, y_out, 1);
}

}
