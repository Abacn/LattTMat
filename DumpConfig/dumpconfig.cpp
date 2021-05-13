//
//  dumpconfig.cpp
//  DumpConfig
//
//  Created by Yi Hu on 8/5/20.
//

#include <algorithm>
#include <fstream>
#include <sstream>
#include "utilityb.hpp"
#include "dumpconfig.hpp"
#include "../TMat/utility.hpp"

namespace DumpMat{

SymDumpConfig::SymDumpConfig(const int N, const uint64_t nstate, const std::string &vecfname/*=""*/):
DumpConfig<double>(N, nstate)
{
  // initialize qmargin in constructor (for full transfer matrices) or later (for reduced transfer matrices)
  if(vecfname.size()>0) initq(vecfname);
}

void SymDumpConfig::initq(const std::string &vecfname, uint64_t fsize)
{
  if(fsize==-1) fsize = nstate;
  vec = readfile(vecfname, fsize);
  qmargin = std::vector<double>(fsize);
  double dtmp = 0.;
  for(uint64_t rp=0; rp<fsize; ++rp)
  {
    dtmp += vec[rp]*vec[rp];
    qmargin[rp] = dtmp;
  }
}

GenDumpConfig::GenDumpConfig(const int N, const uint64_t nstate, const std::string &vecrfname/*=""*/, const std::string &veclfname/*=""*/):
DumpConfig<std::complex<double> >(N, nstate)
{
  // initialize qmargin in constructor (for full transfer matrices) or later (for reduced transfer matrices)
  if(vecrfname.size()>0 && veclfname.size()>0) initq(vecrfname, veclfname);
}

void GenDumpConfig::initq(const std::string &vecrfname, const std::string &veclfname, uint64_t fsize/*=-1*/)
{
  if(fsize==-1) fsize = nstate;
  vec = readfile(vecrfname, fsize);
  std::vector<double> vecl = readfile(veclfname, fsize);
  double dtmp = 0.;
  qmargin = std::vector<double>(fsize);
  for(uint64_t rp=0; rp<fsize; ++rp)
  {
    dtmp += vec[rp]*vecl[rp];
    qmargin[rp] = dtmp;
  }
}

XAT2DDumpConfig::XAT2DDumpConfig(const MParameter &mpara, const std::string &vecfname):
SymDumpConfig(mpara.size, 1ull<<mpara.size, vecfname), mpara(mpara)
{
  interacton[0] = mpara.betaJ;
  interacton[1] = -mpara.kappa*mpara.betaJ;
  interacton[2] = mpara.h*mpara.betaJ;
}


double XAT2DDumpConfig::kernel(SPIN_TYPE a, SPIN_TYPE b) const
{
  return exp(selfintr(a)+crossintr(a, b));
}


double XAT2DDumpConfig::selfintr(SPIN_TYPE ra) const
{
  double ttint = interacton[2]*(size-2.0*Spin::nAntiSpin(ra)); // external field
  // NN same layer
  ttint += interacton[0]*(size-2.0*Spin::nAntiSpin(ra, Spin::offsetState(ra, size, 1)));
  // NN cross layer
  ttint += interacton[1]*(size-2.0*Spin::nAntiSpin(ra, Spin::offsetState(ra, size, 2)));
  return ttint;
}

double XAT2DDumpConfig::crossintr(SPIN_TYPE ra, SPIN_TYPE rb) const
{
  return interacton[0]*(size-2.0*Spin::nAntiSpin(ra, rb));
}

void XAT2DDumpConfig::dump(const int N, const int length, const std::string &fprefix)
{
  using arr_type = std::vector<double>;
  double dtmp, intr1, intr2;
  SPIN_TYPE ra, rc;
  SPIN_TYPE mask=(1ull<<size)-1;
  Cache<SPIN_TYPE, arr_type > qumucache(20);
  arr_type tmpvec(1ull<<size);
  const arr_type *qptr;
  int ntotal=0, ncachmiss=0;
  for(int rp=0; rp<N; ++rp)
  {
    std::stringstream sst;
    sst << fprefix << "_config_" << rp << ".dat";
    std::ofstream ofs(sst.str());
    // first generate endpoint
    dtmp = MyUtility::rand()*qmargin.back();
    auto pos = std::upper_bound(qmargin.begin(), qmargin.end(), dtmp);
    ra = pos-qmargin.begin();
    for(int rq=1; rq<length; ++rq)
    {
      DumpConfig_Base::output(ra, size, ofs);
      intr1 = selfintr(ra);
      qptr = qumucache.get(ra);
      ++ntotal;
      if(!qptr)
      {
        // construct Q(a'|a)
        dtmp = 0.;
        for(rc=0; rc<(1ull<<size); ++rc)
        {
          intr2 = crossintr(ra, rc);
          dtmp += exp(intr1+intr2)*vec[rc];
          tmpvec[rc] = dtmp;
        }
        qumucache.put(ra, tmpvec);
        qptr = &tmpvec;
        ++ncachmiss;
      }
      // then find the next state
      dtmp = MyUtility::rand()*qptr->at(mask);
      auto posb = std::upper_bound(qptr->begin(), qptr->end(), dtmp);
      ra = posb - qptr->begin();
    }
    DumpConfig_Base::output(ra, size, ofs);
    ofs.close();
  }
  std::cout << "Cache hit rate: " << (double)(ntotal-ncachmiss)/ntotal << std::endl;
}


YAT2DDumpConfig::YAT2DDumpConfig(const MParameter &mpara, const std::string &vecrfname, const std::string &veclfname):
GenDumpConfig(mpara.size, 1ull<<(2*mpara.size), vecrfname, veclfname), mpara(mpara)
{
  interacton[0] = mpara.betaJ;
  interacton[1] = -mpara.kappa*mpara.betaJ;
  interacton[2] = mpara.h*mpara.betaJ;
}

double YAT2DDumpConfig::kernel(SPIN_TYPE a, SPIN_TYPE b) const
{
  SPIN_TYPE mask = (1<<size)-1;
  SPIN_TYPE ra=a&mask, rb1=a>>size, rb2=b&mask, rc=b>>size;
  if(rb1!=rb2) return 0.0;
  else return exp(selfintr(ra, rb1)+crossintr(ra, rc));
}

double YAT2DDumpConfig::selfintr(SPIN_TYPE ra, SPIN_TYPE rb) const
{
  double ttint = interacton[2]*(size-2.0*Spin::nAntiSpin(ra)); // external field
  // NN same layer
  ttint += interacton[0]*(size-2.0*Spin::nAntiSpin(ra, Spin::offsetState(ra, size, 1)));
  // NN cross layer
  ttint += interacton[0]*(size-2.0*Spin::nAntiSpin(ra, rb));
  return ttint;
}

double YAT2DDumpConfig::crossintr(SPIN_TYPE ra, SPIN_TYPE rc) const
{
  return interacton[1]*(size-2.0*Spin::nAntiSpin(ra, rc));
}

// dump equilibrium configurations
void YAT2DDumpConfig::dump(const int N, const int length, const std::string &fprefix)
{
  using arr_type = std::vector<double>;
  double dtmp, intr1, intr2;
  SPIN_TYPE thislayer, ra, rb, rc;
  SPIN_TYPE mask=(1ull<<size)-1;
  Cache<SPIN_TYPE, arr_type > qumucache(20);
  arr_type tmpvec(1ull<<size);
  const arr_type *qptr;
  int ntotal=0, ncachmiss=0;
  for(int rp=0; rp<N; ++rp)
  {
    std::stringstream sst;
    sst << fprefix << "_config_" << rp << ".dat";
    std::ofstream ofs(sst.str());
    // first generate endpoint
    dtmp = MyUtility::rand()*qmargin.back();
    auto pos = std::upper_bound(qmargin.begin(), qmargin.end(), dtmp);
    thislayer = pos-qmargin.begin();
    DumpConfig_Base::output(thislayer, size, ofs);
    for(int rq=1; rq<length; ++rq)
    {
      ra = thislayer & mask;
      rb = thislayer >> size;
      DumpConfig_Base::output(rb, size, ofs);
      intr1 = selfintr(ra, rb);
      qptr = qumucache.get(thislayer);
      ++ntotal;
      if(!qptr)
      {
        // construct Q(a'|a)
        dtmp = 0.;
        for(rc=0; rc<(1ull<<size); ++rc)
        {
          intr2 = crossintr(ra, rc);
          dtmp += exp(intr1+intr2)*vec[rc<<size|rb];
          tmpvec[rc] = dtmp;
        }
        qumucache.put(thislayer, tmpvec);
        qptr = &tmpvec;
        ++ncachmiss;
      }
      // then find the next state
      dtmp = MyUtility::rand()*qptr->at(mask);
      auto posb = std::upper_bound(qptr->begin(), qptr->end(), dtmp);
      rc = posb - qptr->begin();
      thislayer = (rc<<size|rb);
    }
    DumpConfig_Base::output(thislayer>>size, size, ofs);
    ofs.close();
  }
  std::cout << "Cache hit rate: " << (double)(ntotal-ncachmiss)/ntotal << std::endl;
}

}
