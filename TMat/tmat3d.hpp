//
//  tmat3d.hpp
//  TMat
//
//  Created by Yi Hu on 1/4/21.
//

#ifndef tmat3d_hpp
#define tmat3d_hpp

#include "tmat.hpp"

namespace TMat{
/**
 Reduced transfer matrix for 3D Ising
 */
class XAMMatVecProd3D: public XMatVecProd3D
{
public:
  XAMMatVecProd3D(const MParameter &mpara, const int extend_cond=0);
  // y_out = M * x_in
  void perform_op(const double *x_in, double *y_out) const override;
  void permute_mult_TL(double *y_tmp) const;
  // override physics
  using SymMatPhys::moments;
  /// \brief Get moments
  virtual int moments(int size, const Scalar *q, double results[TMAT_NMOMENTS]) const override;
  using SymMatPhys::phys_extend;
  virtual int phys_extend(int size, const Scalar *q, std::vector<double> &results) const override;
  /// \brief Get the spin state of irrep state
  std::pair<uint64_t, unsigned short> getIrrepState(SPIN_TYPE idx) const {return { idxlist[idx], duplist[idx] }; }
  /// \brief Get all equivelant spin states of irrep state
  std::vector<SPIN_TYPE> getEqvStates(SPIN_TYPE nowstate) const;
private:
  void genTvec(double normfactorexp=0.0);
  const MParameter mpara;
  // number of bit direct calculation
  const int nleadbit;
  // make use of centrosymmetric or not
  const int censym_flag, ising_flag, rotate_flag;
  double bzwj1[2], bzwj2[2];
  Eigen::VectorXd Tvec;
  std::vector<uint64_t> idxlist;
  std::vector<unsigned short> duplist; // number of equivalent spins. Limit: L<64 as 64*4=256
};

}
#endif /* tmat3d_hpp */
