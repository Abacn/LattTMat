//
//  tmatreduce.hpp
//  TMat
//
//  Created by Yi Hu on 8/14/20.
//

#ifndef tmatreduce_hpp
#define tmatreduce_hpp

#include <cstdint>
#include <cstring>
#include <vector>
#include <array>

#include "tmat.hpp"

namespace TMat{
  
/**
 Reduced transfer matrix for x-ANNNI
 */
class XAMMatVecProd2D: public XMatVecProd2D
{
public:
  XAMMatVecProd2D(const MParameter &mpara, const int extend_cond=0);
  // y_out = M * x_in
  void perform_op(const double *x_in, double *y_out) const override;
  // y_out = (P*T1) * x_in
  void permute_mult_T1(const double *x_in, double *y_out) const;
  void permute_mult_TL(double *y_tmp) const;
  // override physics
  using SymMatPhys::moments;
  /// \brief Get moments
  virtual int moments(int size, const Scalar *q, double results[TMAT_NMOMENTS]) const override;
  using SymMatPhys::phys_extend;
  virtual int phys_extend(int size, const Scalar *q, std::vector<double> &results) const override;
  /// \brief Get the spin state of irrep state
  std::pair<uint64_t, unsigned char> getIrrepState(SPIN_TYPE idx) const {return { idxlist[idx], duplist[idx] }; }
  /// \brief Get all equivelant spin states of irrep state
  std::vector<SPIN_TYPE> getEqvStates(SPIN_TYPE nowstate) const;
private:
  void genTvec(double normfactorexp=0.0);
  const MParameter mpara;
  // number of bit direct calculation
  const int nleadbit;
  // make use of centrosymmetric or not
  const int censym_flag;
  double bzwj1[2], bzwj2[2];
  Eigen::VectorXd Tvec;
  std::vector<uint64_t> idxlist;
  std::vector<unsigned char> duplist; // number of equivalent spins. Limit: L<64 as 64*4=256
};
  
// forward declaration
class YAMMatVecProd2D;

/**
 Reduced transpose transfer matrix for y-ANNNI
 */
class YAMMatVecProd2D_T: public YMatVecProd2D
{
public:
  YAMMatVecProd2D_T(YAMMatVecProd2D const &tmathandle);
  void perform_op(const double *x_in, double *y_out) const override;
private:
  YAMMatVecProd2D const &tmathandle;
};

/**
 Reduced transfer matrix for y-ANNNI, also for all frustrated
 */
class YAMMatVecProd2D: public YMatVecProd2D
{
public:
  YAMMatVecProd2D(const MParameter &mpara, int extend_cond=0);
  // TODO: implement transpose
  using MatVecProd_T = YAMMatVecProd2D_T;
  MatVecProd_T get_t_instance() const {return MatVecProd_T(*this);}
  // y_out = M * x_in
  void perform_op(const double *x_in, double *y_out) const override;
  // y_out = (P*T1) * x_in
  void permute_mult_T1(const double *x_in, double *y_out) const;
  void permute_mult_TL(double *y_tmp) const;
  // override physics
  virtual int moments(int size, const Scalar *ql, const Scalar *qr, double results[TMAT_NMOMENTS]) const override;
  /// \brief Get the spin state of irrep state
  std::pair<uint64_t, unsigned char> getIrrepState(uint64_t idx) const {return { idxlist[idx], duplist[idx] }; }
  /// \brief Get the irrep state index of a spin state
  uint64_t getIrrepStateIdx(uint64_t nowstate) const;
  /// \brief Get all equivelant spin states of irrep state
  std::vector<SPIN_TYPE> getEqvStates(SPIN_TYPE nowstate) const;
private:
  void genTvec(double normfactorexp=0.0);
  void genTvec1(double normfactorexp=0.0);
  const MParameter mpara;
  const int censym_flag, layertype;
  double bzwj1[2], bzwj2[2], bzwj3[2];
  Eigen::VectorXd Tvec;
  std::vector<uint64_t> idxlist, swapidxlist; // irrep status, irrep status idx while swapping layer
  std::vector<unsigned char> duplist; // number of equivalent spins
  // lead layer period. rightest byte: cyclic period;
  // second rightest byte: reversal identity offset;
  // third rightest byte: flip identity offset
  // leftest byte: flip+reversal identity offset
  std::vector<std::array<unsigned char, 4> > leadprdlist;
  // transpose matrix vector product class
  friend class YAMMatVecProd2D_T;
};
  
/**
 Reduced transfer matrix for y-BNNNI
 */
class YBMMatVecProd2D: public YAMMatVecProd2D
{
public:
  YBMMatVecProd2D(const MParameter &mpara, const int extend_cond = 0);
};

/**
 Reduced transfer matrix for DNNI
 \todo Not implemented yet
 */
class XDMMatVecProd2D: public XMatVecProd2D
{
public:
  XDMMatVecProd2D(const MParameter &mpara, int extend_cond=0);
  // y_out = M * x_in
  void perform_op(const double *x_in, double *y_out) const override;
  // y_out = ((P*T1) * x_in)^(size-nleadbit)
  void permute_mult_TL(double *y_tmp, int outbit) const;
  // override physics
  using SymMatPhys::moments;
  /// \brief Get moments
  virtual int moments(int size, const Scalar *q, double results[TMAT_NMOMENTS]) const override;
  using SymMatPhys::phys_extend;
  virtual int phys_extend(int size, const Scalar *q, std::vector<double> &results) const override;
  /// \brief Get the spin state of irrep state
  std::pair<uint64_t, unsigned char> getIrrepState(SPIN_TYPE idx) const {return { idxlist[idx], duplist[idx] }; }
  /// \brief Get all equivelant spin states of irrep state
  std::vector<SPIN_TYPE> getEqvStates(SPIN_TYPE nowstate) const;
private:
  void genTvec(double normfactorexp=0.0);
  const MParameter mpara;
  // number of bit direct calculation
  const int nleadbit;
  // make use of centrosymmetric or not
  const int censym_flag;
  double bzwj1[2], bzwj2[2];
  Eigen::VectorXd Tvec;
  std::vector<uint64_t> idxlist;
  std::vector<unsigned char> duplist; // number of equivalent spins. Limit: L<64 as 64*4=256
};


/** 2D z-BNNNI model reduced
 */
class ZBMMatVecProd2D;

// transpose matrix vector multiplication T^T x
class ZBMMatVecProd2D_T: public ZMatVecProd2D
{
public:
  ZBMMatVecProd2D_T(const ZBMMatVecProd2D &tmathandle);
  // y_out = M * x_in
  void perform_op(const double *x_in, double *y_out) const override;
private:
  const ZBMMatVecProd2D &tmathandle;
};


class ZBMMatVecProd2D: public ZMatVecProd2D
{
public:
  ZBMMatVecProd2D(const MParameter &mpara, int extend_cond=0);
  // y_out = M * x_in
  void perform_op(const double *x_in, double *y_out) const override;
  // opt direction: 0-M v; 1: M^T v
  void perform_op_internal(const double *x_in, double *y_out, int optdirection) const;
  // operating two spin in y_out = ((P*T1) * x_in)^(size-nleadbit)
  void permute_mult_T1(double *x_in, double *y_out, int optdirection) const;
  void permute_mult_T1_last(double *x_in, double *y_out, SPIN_TYPE lastbit, int optdirection) const;
  // get leading boltzmann factor
  double leadingbz(double *y_tmp, SPIN_TYPE nowstate, SPIN_TYPE leadbits, SPIN_TYPE lastbit, int optdirection) const;
  // override physics
  using GenMatPhys::moments;
  /// \brief Get moments
  virtual int moments(int size, const Scalar *ql, const Scalar *qr, double results[TMAT_NMOMENTS]) const override;
  using GenMatPhys::phys_extend;
  virtual int phys_extend(int size, const Scalar *ql, const Scalar *qr, std::vector<double> &results) const override;
  /// \brief Get the spin state of irrep state
  std::pair<uint64_t, unsigned char> getIrrepState(SPIN_TYPE idx) const {return { idxlist[idx], duplist[idx] }; }
  /// \brief Get all equivelant spin states of irrep state
  std::vector<SPIN_TYPE> getEqvStates(SPIN_TYPE nowstate) const;
  using MatVecProd_T = ZBMMatVecProd2D_T;
  MatVecProd_T get_t_instance() const {return ZBMMatVecProd2D_T(*this);}
private:
  void genTvec(double normfactorexp=0.0);
  const MParameter mpara;
  // number of bit direct calculation
  const int nleadbit;
  // make use of centrosymmetric or not
  const int censym_flag;
  double bzwj1[2], bzwj2[2], bzwj3[2], bzwj2b[4], bzwj2c[16];
  Eigen::VectorXd Tvec;
  std::vector<uint64_t> idxlist;
  std::vector<unsigned char> duplist; // number of equivalent spins. Limit: L<64 as 64*4=256
  friend class ZBMMatVecProd2D_T;
};

}
#endif /* tmatreduce_hpp */
