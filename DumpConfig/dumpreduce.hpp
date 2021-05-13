//
//  dumpreduce.hpp
//  DumpConfig
//
//  Created by Yi Hu on 10/11/20.
//

#ifndef dumpreduce_hpp
#define dumpreduce_hpp

#include "dumpconfig.hpp"
#include "../TMat/tmatreduce.hpp"

namespace DumpMat{

/// \brief x-ANNNI model reduced transfer matrix dump configurations
class XAM2DDumpConfig: public SymDumpConfig{
public:
  XAM2DDumpConfig(const MParameter &mpara, const std::string &vecfname, const int extend_cond=0);
  double kernel(SPIN_TYPE a, SPIN_TYPE b) const;
  double selfintr(SPIN_TYPE ra) const; // self interaction of a unit
  double crossintr(SPIN_TYPE ra, SPIN_TYPE rb) const; // cross interaction of a unit
  void dump(const int N, const int length, const std::string &fprefix);
  virtual ~XAM2DDumpConfig(){}
private:
  const MParameter mpara;
  double interacton[3];
  const TMat::XAMMatVecProd2D mathandle;
};

class YAM2DDumpConfig: public GenDumpConfig{
public:
  YAM2DDumpConfig(const MParameter &mpara, const std::string &vecrfname, const std::string &veclfname, const int extend_cond=0);
  double kernel(SPIN_TYPE a, SPIN_TYPE b) const;
  double selfintr(SPIN_TYPE ra, SPIN_TYPE rb) const; // self interaction of a unit
  double crossintr(SPIN_TYPE ra, SPIN_TYPE rc) const; // cross interaction of a unit
  void dump(const int N, const int length, const std::string &fprefix);
  virtual ~YAM2DDumpConfig(){}
private:
  const MParameter mpara;
  double interacton[3];
  const TMat::YAMMatVecProd2D mathandle;
};

/// \brief y-full frustrated model reduced transfer matrix dump configurations
class YBM2DDumpConfig: public GenDumpConfig{
public:
  YBM2DDumpConfig(const MParameter &mpara, const std::string &vecrfname, const std::string &veclfname, const int extend_cond=0);
  double kernel(SPIN_TYPE a, SPIN_TYPE b) const;
  double selfintr(SPIN_TYPE ra, SPIN_TYPE rb) const; // self interaction of a unit
  double crossintr(SPIN_TYPE ra, SPIN_TYPE rc) const; // cross interaction of a unit
  void dump(const int N, const int length, const std::string &fprefix);
  virtual ~YBM2DDumpConfig(){}
private:
  const MParameter mpara;
  double interacton[4];
  const TMat::YBMMatVecProd2D mathandle;
};

/// \brief z-full frustrated model reduced transfer matrix dump configurations
class ZBM2DDumpConfig: public GenDumpConfig{
public:
  ZBM2DDumpConfig(const MParameter &mpara, const std::string &vecrfname, const std::string &veclfname, const int extend_cond=0);
  double kernel(SPIN_TYPE ra, SPIN_TYPE rb) const;
  double selfintr(SPIN_TYPE ra) const; // self interaction of a unit
  double crossintr(SPIN_TYPE ra, SPIN_TYPE rb) const; // cross interaction of a unit
  void dump(const int N, const int length, const std::string &fprefix);
  virtual ~ZBM2DDumpConfig(){}
private:
  const MParameter mpara;
  double interacton[4];
  const TMat::ZBMMatVecProd2D mathandle;
};

/// \brief x-DNNI model reduced transfer matrix dump configurations
class XDM2DDumpConfig: public SymDumpConfig{
public:
  XDM2DDumpConfig(const MParameter &mpara, const std::string &vecfname, const int extend_cond=0);
  double kernel(SPIN_TYPE a, SPIN_TYPE b) const;
  double selfintr(SPIN_TYPE ra) const; // self interaction of a unit
  double crossintr(SPIN_TYPE ra, SPIN_TYPE rb) const; // cross interaction of a unit
  void dump(const int N, const int length, const std::string &fprefix);
  virtual ~XDM2DDumpConfig(){}
private:
  const MParameter mpara;
  double interacton[3];
  const TMat::XDMMatVecProd2D mathandle;
};

}

#endif /* dumpreduce_hpp */
