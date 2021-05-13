//
//  dumpreduce.cpp
//  DumpConfig
//
//  Created by Yi Hu on 10/11/20.
//
#include <algorithm>
#include <fstream>
#include <sstream>
#include "utilityb.hpp"
#include "dumpreduce.hpp"

#include "../TMat/utility.hpp"

namespace DumpMat{

XAM2DDumpConfig::XAM2DDumpConfig(const MParameter &mpara, const std::string &vecfname, const int extend_cond/*=0*/):
SymDumpConfig(mpara.size, 1ull<<mpara.size, ""), mpara(mpara), mathandle(mpara, extend_cond)
{
  interacton[0] = mpara.betaJ;
  interacton[1] = -mpara.kappa*mpara.betaJ;
  interacton[2] = mpara.h*mpara.betaJ;
  initq(vecfname, mathandle.cols());
}

double XAM2DDumpConfig::kernel(SPIN_TYPE a, SPIN_TYPE b) const
{
  return exp(0.5*(selfintr(a)+selfintr(b))+crossintr(a, b));
}


double XAM2DDumpConfig::selfintr(SPIN_TYPE ra) const
{
  double ttint = interacton[2]*(size-2.0*Spin::nAntiSpin(ra)); // external field
  // NN same layer
  ttint += interacton[0]*(size-2.0*Spin::nAntiSpin(ra, Spin::offsetState(ra, size, 1)));
  // NN cross layer
  ttint += interacton[1]*(size-2.0*Spin::nAntiSpin(ra, Spin::offsetState(ra, size, 2)));
  return ttint;
}

double XAM2DDumpConfig::crossintr(SPIN_TYPE ra, SPIN_TYPE rb) const
{
  return interacton[0]*(size-2.0*Spin::nAntiSpin(ra, rb));
}

void XAM2DDumpConfig::dump(const int N, const int length, const std::string &fprefix)
{
  using arr_type = std::vector<double>;
  double dtmp;
  SPIN_TYPE ra, qa, qc;
  Cache<SPIN_TYPE, arr_type > qumucache(20);
  arr_type tmpvec(qmargin.size());
  const arr_type *qptr;
  int ntotal=0, ncachmiss=0;
  if(qmargin.size()==0) return; // not initialized yet
  assert(qmargin.size() == mathandle.rows());
  double *x_in = new double[qmargin.size()], *y_out = new double[qmargin.size()];
  for(int rp=0; rp<N; ++rp)
  {
    std::stringstream sst;
    sst << fprefix << "_config_" << rp << ".dat";
    std::ofstream ofs(sst.str());
    // first generate endpoint
    dtmp = MyUtility::rand()*qmargin.back();
    auto pos = std::upper_bound(qmargin.begin(), qmargin.end(), dtmp);
    qa = pos-qmargin.begin();
    auto qapair = mathandle.getIrrepState(qa);
    ra = qapair.first;
    for(int rq=1; rq<length; ++rq)
    {
      DumpConfig_Base::output(ra, size, ofs);
      qptr = qumucache.get(qa);
      ++ntotal;
      if(!qptr)
      {
        // construct Q(a'|a)
        // initialize in vec
        for(qc=0; qc<qmargin.size(); ++qc)
        {
          x_in[qc] = 0.;
        }
        x_in[qa] = 1.;
        mathandle.perform_op(x_in, y_out);
        dtmp = 0.;
        for(qc=0; qc<qmargin.size(); ++qc)
        {
          dtmp += y_out[qc]*vec[qc];
          tmpvec[qc] = dtmp;
        }
        qumucache.put(qa, tmpvec);
        qptr = &tmpvec;
        ++ncachmiss;
      }
      // then find the next state
      dtmp = MyUtility::rand()*qptr->back();
      auto posb = std::upper_bound(qptr->begin(), qptr->end(), dtmp);
      qa = posb - qptr->begin();
      // find next states from equivalent states aka qa
      auto qpair = mathandle.getIrrepState(qa);
      SPIN_TYPE rc1 = qpair.first;
      auto eqstates = mathandle.getEqvStates(rc1);
      assert(eqstates.size() == qpair.second);
      dtmp = 0.;
      std::vector<double> qcond(eqstates.size());
      for(qc=0; qc<eqstates.size(); ++qc)
      {
        dtmp += exp(crossintr(ra, eqstates[qc]));
        qcond[qc] = dtmp;
      }
      dtmp = MyUtility::rand()*qcond.back();
      auto posc = std::upper_bound(qcond.begin(), qcond.end(), dtmp);
      qc = posc - qcond.begin();
      ra = eqstates[qc];
    }
    DumpConfig_Base::output(ra, size, ofs);
    ofs.close();
  }
  delete [] x_in;
  delete [] y_out;
  std::cout << "Cache hit rate: " << (double)(ntotal-ncachmiss)/ntotal << std::endl;
}

YAM2DDumpConfig::YAM2DDumpConfig(const MParameter &mpara, const std::string &vecrfname, const std::string &veclfname, const int extend_cond/*=0*/):
GenDumpConfig(mpara.size, 1ull<<mpara.size, "", ""), mpara(mpara), mathandle(mpara, extend_cond)
{
  interacton[0] = mpara.betaJ;
  interacton[1] = -mpara.kappa*mpara.betaJ;
  interacton[2] = -mpara.h*mpara.betaJ;
  initq(vecrfname, veclfname, mathandle.cols());
}

double YAM2DDumpConfig::kernel(SPIN_TYPE a, SPIN_TYPE b) const
{
  SPIN_TYPE mask = (1<<size)-1;
  SPIN_TYPE ra=a&mask, rb1=a>>size, rb2=b&mask, rc=b>>size;
  if(rb1!=rb2) return 0.0;
  else return exp(0.5*(selfintr(rb1, ra)+selfintr(rc, rb2))+crossintr(ra, rc));
}

double YAM2DDumpConfig::selfintr(SPIN_TYPE ra, SPIN_TYPE rb) const
{
  double ttint = interacton[2]*(size-2.0*Spin::nAntiSpin(ra)); // external field
  // NN same layer
  ttint += interacton[0]*(size-2.0*Spin::nAntiSpin(ra, Spin::offsetState(ra, size, 1)));
  // NN cross layer
  ttint += interacton[0]*(size-2.0*Spin::nAntiSpin(ra, rb));
  return ttint;
}

double YAM2DDumpConfig::crossintr(SPIN_TYPE ra, SPIN_TYPE rc) const
{
  return interacton[1]*(size-2.0*Spin::nAntiSpin(ra, rc));
}

void YAM2DDumpConfig::dump(const int N, const int length, const std::string &fprefix)
{
  using arr_type = std::vector<double>;
  double dtmp;
  SPIN_TYPE thislayer, ra, rb, rc, qa, qc;
  SPIN_TYPE mask=(1ull<<size)-1;
  Cache<SPIN_TYPE, arr_type > qumucache(20);
  arr_type tmpvec(qmargin.size());
  const arr_type *qptr;
  int ntotal=0, ncachmiss=0;
  if(qmargin.size()==0) return; // not initialized yet
  assert(qmargin.size() == mathandle.rows());
  double *x_in = new double[qmargin.size()], *y_out = new double[qmargin.size()];
  for(int rp=0; rp<N; ++rp)
  {
    std::stringstream sst;
    sst << fprefix << "_config_" << rp << ".dat";
    std::ofstream ofs(sst.str());
    // first generate endpoint
    dtmp = MyUtility::rand()*qmargin.back();
    auto pos = std::upper_bound(qmargin.begin(), qmargin.end(), dtmp);
    qa = pos-qmargin.begin();
    auto qapair = mathandle.getIrrepState(qa);
    thislayer = qapair.first;
    for(int rq=1; rq<length; ++rq)
    {
      ra = thislayer & mask;
      rb = thislayer >> size;
      DumpConfig_Base::output(ra, size, ofs);
      qptr = qumucache.get(qa);
      ++ntotal;
      if(!qptr)
      {
        // construct Q(a'|a)
        // initialize in vec
        for(qc=0; qc<qmargin.size(); ++qc)
        {
          x_in[qc] = 0.;
        }
        x_in[qa] = 1.;
        mathandle.perform_op(x_in, y_out);
        dtmp = 0.;
        for(qc=0; qc<qmargin.size(); ++qc)
        {
          dtmp += y_out[qc]*vec[qc];
          tmpvec[qc] = dtmp;
        }
        qumucache.put(qa, tmpvec);
        qptr = &tmpvec;
        ++ncachmiss;
      }
      // then find the next state
      dtmp = MyUtility::rand()*qptr->back();
      auto posb = std::upper_bound(qptr->begin(), qptr->end(), dtmp);
      qa = posb - qptr->begin();
      // find next states from equivalent states aka qa
      auto qpair = mathandle.getIrrepState(qa);
      SPIN_TYPE rc1 = qpair.first;
      auto eqstates = mathandle.getEqvStates(rc1);
      assert(eqstates.size() == qpair.second);
      dtmp = 0.;
      std::vector<double> qcond(eqstates.size());
      for(qc=0; qc<eqstates.size(); ++qc)
      {
        SPIN_TYPE nextlayer = eqstates[qc];
        if((nextlayer & mask) == rb)
        {
          rc = nextlayer >> size;
          dtmp += exp(crossintr(ra, rc));
        }
        qcond[qc] = dtmp;
      }
      dtmp = MyUtility::rand()*qcond.back();
      auto posc = std::upper_bound(qcond.begin(), qcond.end(), dtmp);
      qc = posc - qcond.begin();
      thislayer = eqstates[qc];
    }
    DumpConfig_Base::output(thislayer & mask, size, ofs);
    ofs.close();
  }
  delete [] x_in;
  delete [] y_out;
  std::cout << "Cache hit rate: " << (double)(ntotal-ncachmiss)/ntotal << std::endl;
}

YBM2DDumpConfig::YBM2DDumpConfig(const MParameter &mpara, const std::string &vecrfname, const std::string &veclfname, const int extend_cond/*=0*/):
GenDumpConfig(mpara.size, 1ull<<mpara.size, "", ""), mpara(mpara), mathandle(mpara, extend_cond)
{
  interacton[0] = mpara.betaJ;
  interacton[1] = -mpara.kappa*mpara.betaJ;
  interacton[2] = -mpara.kappa2*mpara.betaJ;
  interacton[3] = mpara.h*mpara.betaJ;
  initq(vecrfname, veclfname, mathandle.cols());
}

double YBM2DDumpConfig::kernel(SPIN_TYPE a, SPIN_TYPE b) const
{
  SPIN_TYPE mask = (1<<size)-1;
  SPIN_TYPE ra=a&mask, rb1=a>>size, rb2=b&mask, rc=b>>size;
  if(rb1!=rb2) return 0.0;
  else return exp(0.5*(selfintr(rb1, ra)+selfintr(rc, rb2))+crossintr(ra, rc));
}

double YBM2DDumpConfig::selfintr(SPIN_TYPE ra, SPIN_TYPE rb) const
{
  double ttint = interacton[3]*(size-2.0*Spin::nAntiSpin(ra)); // external field
  // NN same layer
  ttint += interacton[0]*(size-2.0*Spin::nAntiSpin(ra, Spin::offsetState(ra, size, 1)));
  // NN cross layer
  ttint += interacton[0]*(size-2.0*Spin::nAntiSpin(ra, rb));
  // NNN same layer
  ttint += interacton[1]*(size-2.0*Spin::nAntiSpin(ra, Spin::offsetState(ra, size, 2)));
  // diagonal cross layer
  ttint += interacton[2]*(size-2.0*Spin::nAntiSpin(ra, Spin::offsetState(rb, size, 1)));
  ttint += interacton[2]*(size-2.0*Spin::nAntiSpin(ra, Spin::lffsetState(rb, size, 1)));
  return ttint;
}

double YBM2DDumpConfig::crossintr(SPIN_TYPE ra, SPIN_TYPE rc) const
{
  return interacton[1]*(size-2.0*Spin::nAntiSpin(ra, rc));
}

void YBM2DDumpConfig::dump(const int N, const int length, const std::string &fprefix)
{
  using arr_type = std::vector<double>;
  double dtmp;
  SPIN_TYPE thislayer, ra, rb, rc, qa, qc;
  SPIN_TYPE mask=(1ull<<size)-1;
  Cache<SPIN_TYPE, arr_type > qumucache(20);
  arr_type tmpvec(qmargin.size());
  const arr_type *qptr;
  int ntotal=0, ncachmiss=0;
  if(qmargin.size()==0) return; // not initialized yet
  assert(qmargin.size() == mathandle.rows());
  double *x_in = new double[qmargin.size()], *y_out = new double[qmargin.size()];
  for(int rp=0; rp<N; ++rp)
  {
    std::stringstream sst;
    sst << fprefix << "_config_" << rp << ".dat";
    std::ofstream ofs(sst.str());
    // first generate endpoint
    dtmp = MyUtility::rand()*qmargin.back();
    auto pos = std::upper_bound(qmargin.begin(), qmargin.end(), dtmp);
    qa = pos-qmargin.begin();
    auto qapair = mathandle.getIrrepState(qa);
    thislayer = qapair.first;
    for(int rq=1; rq<length; ++rq)
    {
      ra = thislayer & mask;
      rb = thislayer >> size;
      DumpConfig_Base::output(ra, size, ofs);
      qptr = qumucache.get(qa);
      ++ntotal;
      if(!qptr)
      {
        // construct Q(a'|a)
        // initialize in vec
        for(qc=0; qc<qmargin.size(); ++qc)
        {
          x_in[qc] = 0.;
        }
        x_in[qa] = 1.;
        mathandle.perform_op(x_in, y_out);
        dtmp = 0.;
        for(qc=0; qc<qmargin.size(); ++qc)
        {
          dtmp += y_out[qc]*vec[qc];
          tmpvec[qc] = dtmp;
        }
        qumucache.put(qa, tmpvec);
        qptr = &tmpvec;
        ++ncachmiss;
      }
      // then find the next state
      dtmp = MyUtility::rand()*qptr->back();
      auto posb = std::upper_bound(qptr->begin(), qptr->end(), dtmp);
      qa = posb - qptr->begin();
      // find next states from equivalent states aka qa
      auto qpair = mathandle.getIrrepState(qa);
      SPIN_TYPE rc1 = qpair.first;
      auto eqstates = mathandle.getEqvStates(rc1);
      assert(eqstates.size() == qpair.second);
      dtmp = 0.;
      std::vector<double> qcond(eqstates.size());
      for(qc=0; qc<eqstates.size(); ++qc)
      {
        SPIN_TYPE nextlayer = eqstates[qc];
        if((nextlayer & mask) == rb)
        {
          rc = nextlayer >> size;
          dtmp += exp(crossintr(ra, rc));
        }
        qcond[qc] = dtmp;
      }
      dtmp = MyUtility::rand()*qcond.back();
      auto posc = std::upper_bound(qcond.begin(), qcond.end(), dtmp);
      qc = posc - qcond.begin();
      thislayer = eqstates[qc];
    }
    DumpConfig_Base::output(thislayer & mask, size, ofs);
    ofs.close();
  }
  delete [] x_in;
  delete [] y_out;
  std::cout << "Cache hit rate: " << (double)(ntotal-ncachmiss)/ntotal << std::endl;
}

ZBM2DDumpConfig::ZBM2DDumpConfig(const MParameter &mpara, const std::string &vecrfname, const std::string &veclfname, const int extend_cond/*=0*/):
GenDumpConfig(mpara.size, 1ull<<mpara.size, "", ""), mpara(mpara), mathandle(mpara, extend_cond)
{
  interacton[0] = mpara.betaJ;
  interacton[1] = -mpara.kappa*mpara.betaJ;
  interacton[2] = -mpara.kappa2*mpara.betaJ;
  interacton[3] = mpara.h*mpara.betaJ;
  initq(vecrfname, veclfname, mathandle.cols());
}

double ZBM2DDumpConfig::kernel(SPIN_TYPE ra, SPIN_TYPE rb) const
{
  return exp( 0.5*(selfintr(ra)+selfintr(rb))+crossintr(ra, rb));
}

double ZBM2DDumpConfig::selfintr(SPIN_TYPE ra) const
{
  double ttint = interacton[3]*(size-2.0*Spin::nAntiSpin(ra)); // external field
  // NN same layer
  ttint += interacton[0]*(size-2.0*Spin::nAntiSpin(ra, Spin::offsetState(ra, size, 1)));
  // DNN same layer
  ttint += interacton[2]*(size-2.0*Spin::nAntiSpin(ra, Spin::offsetState(ra, size, 2)));
  return ttint;
}

double ZBM2DDumpConfig::crossintr(SPIN_TYPE ra, SPIN_TYPE rb) const // TODO
{
  double ttint = interacton[2]*(size-2.0*Spin::nAntiSpin(ra, rb)); // DNNI interaction
  // BNNNI interaction
  ttint += interacton[1]*(size-2.0*Spin::nAntiSpin(ra, Spin::offsetState(rb, size, 2) ));
  ttint += interacton[1]*(size-2.0*Spin::nAntiSpin(ra, Spin::lffsetState(rb, size, 2) ));
  // Ising interaction
  int halfsize = size/2, nanti = 0;
  rb = rb << 1 | rb >> (size-1);
  for(int rp=0; rp<halfsize; ++rp)
  {
    nanti += (ra ^ rb) & 1;
    rb >>= 2;
    nanti += (ra ^ rb) & 1;
    ra >>= 2;
  }
  ttint += interacton[0]*(size-2.0*nanti);
  return ttint;
}

void ZBM2DDumpConfig::dump(const int N, const int length, const std::string &fprefix) // TODO
{
  using arr_type = std::vector<double>;
  double dtmp;
  SPIN_TYPE ra, qa, qc;
  Cache<SPIN_TYPE, arr_type > qumucache(20);
  arr_type tmpvec(qmargin.size());
  const arr_type *qptr;
  int ntotal=0, ncachmiss=0;
  if(qmargin.size()==0) return; // not initialized yet
  assert(qmargin.size() == mathandle.rows());
  double *x_in = new double[qmargin.size()], *y_out = new double[qmargin.size()];
  for(int rp=0; rp<N; ++rp)
  {
    std::stringstream sst;
    sst << fprefix << "_config_" << rp << ".dat";
    std::ofstream ofs(sst.str());
    // first generate endpoint
    dtmp = MyUtility::rand()*qmargin.back();
    auto pos = std::upper_bound(qmargin.begin(), qmargin.end(), dtmp);
    qa = pos-qmargin.begin();
    auto qapair = mathandle.getIrrepState(qa);
    ra = qapair.first;
    for(int rq=1; rq<length; ++rq)
    {
      DumpConfig_Base::output(ra, size, ofs);
      qptr = qumucache.get(qa);
      ++ntotal;
      if(!qptr)
      {
        // construct Q(a'|a)
        // initialize in vec
        for(qc=0; qc<qmargin.size(); ++qc)
        {
          x_in[qc] = 0.;
        }
        x_in[qa] = 1.;
        mathandle.perform_op(x_in, y_out);
        dtmp = 0.;
        for(qc=0; qc<qmargin.size(); ++qc)
        {
          dtmp += y_out[qc]*vec[qc];
          tmpvec[qc] = dtmp;
        }
        qumucache.put(qa, tmpvec);
        qptr = &tmpvec;
        ++ncachmiss;
      }
      // then find the next state
      dtmp = MyUtility::rand()*qptr->back();
      auto posb = std::upper_bound(qptr->begin(), qptr->end(), dtmp);
      qa = posb - qptr->begin();
      // find next states from equivalent states aka qa
      auto qpair = mathandle.getIrrepState(qa);
      SPIN_TYPE rc1 = qpair.first;
      auto eqstates = mathandle.getEqvStates(rc1);
      assert(eqstates.size() == qpair.second);
      dtmp = 0.;
      std::vector<double> qcond(eqstates.size());
      for(qc=0; qc<eqstates.size(); ++qc)
      {
        dtmp += exp(crossintr(ra, eqstates[qc]));
        qcond[qc] = dtmp;
      }
      dtmp = MyUtility::rand()*qcond.back();
      auto posc = std::upper_bound(qcond.begin(), qcond.end(), dtmp);
      qc = posc - qcond.begin();
      ra = eqstates[qc];
    }
    DumpConfig_Base::output(ra, size, ofs);
    ofs.close();
  }
  delete [] x_in;
  delete [] y_out;
  std::cout << "Cache hit rate: " << (double)(ntotal-ncachmiss)/ntotal << std::endl;
}

XDM2DDumpConfig::XDM2DDumpConfig(const MParameter &mpara, const std::string &vecfname, const int extend_cond/*=0*/):
SymDumpConfig(mpara.size, 1ull<<mpara.size, ""), mpara(mpara), mathandle(mpara, extend_cond)
{
  interacton[0] = mpara.betaJ;
  interacton[1] = -mpara.kappa*mpara.betaJ;
  interacton[2] = mpara.h*mpara.betaJ;
  initq(vecfname, mathandle.cols());
}

double XDM2DDumpConfig::kernel(SPIN_TYPE a, SPIN_TYPE b) const
{
  return exp(0.5*(selfintr(a)+selfintr(b))+crossintr(a, b));
}


double XDM2DDumpConfig::selfintr(SPIN_TYPE ra) const
{
  double ttint = interacton[2]*(size-2.0*Spin::nAntiSpin(ra)); // external field
  // NN same layer
  ttint += interacton[0]*(size-2.0*Spin::nAntiSpin(ra, Spin::offsetState(ra, size, 1)));
  return ttint;
}

double XDM2DDumpConfig::crossintr(SPIN_TYPE ra, SPIN_TYPE rb) const
{
  return interacton[0]*(size-2.0*Spin::nAntiSpin(ra, rb))
         +interacton[1]*(size-2.0*Spin::nAntiSpin(ra, Spin::offsetState(rb, size, 1) ))
         +interacton[1]*(size-2.0*Spin::nAntiSpin(ra, Spin::lffsetState(rb, size, 1) ));
}

void XDM2DDumpConfig::dump(const int N, const int length, const std::string &fprefix)
{
  using arr_type = std::vector<double>;
  double dtmp;
  SPIN_TYPE ra, qa, qc;
  Cache<SPIN_TYPE, arr_type > qumucache(20);
  arr_type tmpvec(qmargin.size());
  const arr_type *qptr;
  int ntotal=0, ncachmiss=0;
  if(qmargin.size()==0) return; // not initialized yet
  assert(qmargin.size() == mathandle.rows());
  double *x_in = new double[qmargin.size()], *y_out = new double[qmargin.size()];
  for(int rp=0; rp<N; ++rp)
  {
    std::stringstream sst;
    sst << fprefix << "_config_" << rp << ".dat";
    std::ofstream ofs(sst.str());
    // first generate endpoint
    dtmp = MyUtility::rand()*qmargin.back();
    auto pos = std::upper_bound(qmargin.begin(), qmargin.end(), dtmp);
    qa = pos-qmargin.begin();
    auto qapair = mathandle.getIrrepState(qa);
    ra = qapair.first;
    for(int rq=1; rq<length; ++rq)
    {
      DumpConfig_Base::output(ra, size, ofs);
      qptr = qumucache.get(qa);
      ++ntotal;
      if(!qptr)
      {
        // construct Q(a'|a)
        // initialize in vec
        for(qc=0; qc<qmargin.size(); ++qc)
        {
          x_in[qc] = 0.;
        }
        x_in[qa] = 1.;
        mathandle.perform_op(x_in, y_out);
        dtmp = 0.;
        for(qc=0; qc<qmargin.size(); ++qc)
        {
          dtmp += y_out[qc]*vec[qc];
          tmpvec[qc] = dtmp;
        }
        qumucache.put(qa, tmpvec);
        qptr = &tmpvec;
        ++ncachmiss;
      }
      // then find the next state
      dtmp = MyUtility::rand()*qptr->back();
      auto posb = std::upper_bound(qptr->begin(), qptr->end(), dtmp);
      qa = posb - qptr->begin();
      // find next states from equivalent states aka qa
      auto qpair = mathandle.getIrrepState(qa);
      SPIN_TYPE rc1 = qpair.first;
      auto eqstates = mathandle.getEqvStates(rc1);
      assert(eqstates.size() == qpair.second);
      dtmp = 0.;
      std::vector<double> qcond(eqstates.size());
      for(qc=0; qc<eqstates.size(); ++qc)
      {
        dtmp += exp(crossintr(ra, eqstates[qc]));
        qcond[qc] = dtmp;
      }
      dtmp = MyUtility::rand()*qcond.back();
      auto posc = std::upper_bound(qcond.begin(), qcond.end(), dtmp);
      qc = posc - qcond.begin();
      ra = eqstates[qc];
    }
    DumpConfig_Base::output(ra, size, ofs);
    ofs.close();
  }
  delete [] x_in;
  delete [] y_out;
  std::cout << "Cache hit rate: " << (double)(ntotal-ncachmiss)/ntotal << std::endl;
}

}
