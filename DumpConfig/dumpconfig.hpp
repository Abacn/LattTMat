//
//  dumpconfig.hpp
//  DumpConfig
//
//  Created by Yi Hu on 8/5/20.
//

#ifndef dumpconfig_hpp
#define dumpconfig_hpp

#include <iostream>
#include <fstream>
#include <exception>
#include <cstdint>
#include <string>
#include <vector>
#include <complex>
#include <random>

#include <sys/stat.h>

#include "../TMat/commondef.hpp"

namespace DumpMat{

/// \brief Base interface
class DumpConfig_Base{
public:
  virtual void dump(const int N, const int length, const std::string &fprefix) = 0;
  virtual double kernel(SPIN_TYPE a, SPIN_TYPE b) const = 0;
  virtual ~DumpConfig_Base(){}
  static void output(SPIN_TYPE state, int nspin, std::ofstream &ofs)
  {
    while(--nspin >= 0)
    {
      ofs << (state&1);
      state>>=1;
    }
    ofs << std::endl;//"\n";
  }

};

/// \brief Base class of dump configuration
/// \tparam S scalar typename
template <typename S>
class DumpConfig: public DumpConfig_Base{
public:
  using Scalar = S;
  DumpConfig(const int N, const uint64_t nstate): size(N), nstate(nstate){}
  virtual ~DumpConfig(){}
protected:
  const int size;
  const uint64_t nstate;
  /// \brief Read eigenvector from file
  std::vector<double> readfile(std::string fname, uint64_t fsize=-1)
  {
    struct stat stat_buf;
    int rc = stat(fname.c_str(), &stat_buf);
    bool is_complex = false;
    if(fsize==-1) fsize = nstate;
    if(rc) {throw std::runtime_error("Open \""+fname+"\" failed");}
    if(stat_buf.st_size < sizeof(double)*fsize) {throw std::runtime_error("File \""+fname+"\" size incorrect");}
    else if(!std::is_same<double,Scalar>::value && stat_buf.st_size >= sizeof(Scalar)*fsize) is_complex = true;
    std::ifstream ifs(fname, std::ios::in | std::ios::binary);
    std::vector<double> tmpvec(fsize);
    double tmp;
    for(uint64_t rp=0; rp<fsize; ++rp)
    {
      ifs.read((char*)&tmp, sizeof(double));
      tmpvec[rp] = fabs(tmp);
      if(is_complex) ifs.read((char*)&tmp, sizeof(double));
    }
    if(!ifs){throw std::runtime_error("Read from \""+fname+"\" failed");}
    return tmpvec;
  }
  std::vector<double> qmargin, vec;
};

/// \brief Base class of symmetric matrix
class SymDumpConfig: public DumpConfig<double>{
public:
  SymDumpConfig(const int N, const uint64_t nstate, const std::string &vecfname="");
  /// \brief Initialize qmargin from eigenvectors read from file
  void initq(const std::string &vecfname, uint64_t fsize=-1);
  virtual ~SymDumpConfig(){}
protected:
  
};

/// \brief Base class of generic (non-symmetric) matrix
class GenDumpConfig: public DumpConfig<std::complex<double> >{
public:
  GenDumpConfig(const int N, const uint64_t nstate, const std::string &vecrfname="", const std::string &veclfname="");
  /// \brief Initialize qmargin from eigenvectors read from file
  void initq(const std::string &vecrfname="", const std::string &veclfname="", uint64_t fsize=-1);
  virtual ~GenDumpConfig(){}
protected:
};

/// \brief Dump configuration of X ANNNI model
class XAT2DDumpConfig: public SymDumpConfig{
public:
  XAT2DDumpConfig(const MParameter &mpara, const std::string &vecfname);
  double kernel(SPIN_TYPE a, SPIN_TYPE b) const;
  double selfintr(SPIN_TYPE ra) const; // self interaction of a unit
  double crossintr(SPIN_TYPE ra, SPIN_TYPE rb) const; // cross interaction of a unit
  void dump(const int N, const int length, const std::string &fprefix);
  virtual ~XAT2DDumpConfig(){}
private:
  const MParameter mpara;
  double interacton[3];
};

/// \brief Dump configuration of Y ANNNI model
class YAT2DDumpConfig: public GenDumpConfig{
public:
  YAT2DDumpConfig(const MParameter &mpara, const std::string &vecrfname, const std::string &veclfname);
  double kernel(SPIN_TYPE a, SPIN_TYPE b) const;
  double selfintr(SPIN_TYPE ra, SPIN_TYPE rb) const; // self interaction of a unit
  double crossintr(SPIN_TYPE ra, SPIN_TYPE rc) const; // cross interaction of a unit
  void dump(const int N, const int length, const std::string &fprefix);
  virtual ~YAT2DDumpConfig(){}
private:
  const MParameter mpara;
  double interacton[3];
};

}
#endif /* dumpconfig_hpp */
