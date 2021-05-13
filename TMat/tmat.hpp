//
//  xtmat.hpp
//  XTMat
//
//  Created by Yi Hu on 1/20/20.
//

#ifndef tmat_hpp
#define tmat_hpp

#ifndef DEBUG
#define EIGEN_NO_DEBUG
#endif

#include <complex>
#include <Eigen/Core>
#include <vector>

#include "commondef.hpp"
#include "utility.hpp"

/** transfer matrix namespace
 */
namespace TMat
{
  
// functions
/** generate unsymmetrized full x-axial NNNI transfer matrix
 */
Eigen::MatrixXd genXATmat2D(const MParameter &mpara);

/** generate full y-axial NNNI transfer matrix
 */
Eigen::MatrixXd genYATmat2D(const MParameter &mpara);

/** generate unsymmetrized full DNNI transfer matrix
 */
Eigen::MatrixXd genXDTmat2D(const MParameter &mpara);

/** Physical quantity calculator base class
 */
template <typename S>
class MatPhys
{
public:
  using Scalar = S;
  MatPhys(): boundary_cond(0) {}
  MatPhys(int boundary_cond): boundary_cond(boundary_cond) {}
  /// \brief Calculate first TMAT_NMOMENTS moments
  virtual int moments(int size, const Scalar *ql, const Scalar *qr, double results[TMAT_NMOMENTS]) const = 0;
  /// \brief Other physical properties
  virtual int phys_extend(int size, const Scalar *ql, const Scalar *qr, std::vector<double> &results) const = 0;
protected:
  // 0-periodic boundary; 1-open boundary
  const int boundary_cond;
};

/** Physical quantity calculator. Vector size 2^size
 */
class SymMatPhys: public MatPhys<double>
{
public:
  SymMatPhys(int boundary_cond=0);
  virtual int moments(int size, const Scalar *ql, const Scalar *qr, double results[TMAT_NMOMENTS]) const;
  virtual int moments(int size, const Scalar *q, double results[TMAT_NMOMENTS]) const;
  virtual int phys_extend(int size, const Scalar *ql, const Scalar *qr, std::vector<double> &results) const;
  virtual int phys_extend(int size, const Scalar *q, std::vector<double> &results) const;
};

/** Physical quantity calculator. Vector size 4^size
 */
class GenMatPhys: public MatPhys<std::complex<double> >
{
public:
  GenMatPhys(int boundary_cond=0);
  virtual int moments(int size, const Scalar *ql, const Scalar *qr, double results[TMAT_NMOMENTS]) const;
  virtual int phys_extend(int size, const Scalar *ql, const Scalar *qr, std::vector<double> &results) const;
};

/** transfer matrix base class
 */
template<typename S>
class MatVecProd
{
public:
  virtual Eigen::Index rows() const {return ncolsize;}
  virtual Eigen::Index cols() const {return ncolsize;}
  virtual int suggested_nev(int neig) const {return neig>4? 2*neig+1: 2*neig+3;}
  virtual void perform_op(const double *x_in, double *y_out) const = 0;
  virtual Eigen::MatrixXd raw_Tmat() const
  {
    Eigen::MatrixXd mat(ncolsize, ncolsize);
    double *x_in = new double[ncolsize];
    double *y_out = new double[ncolsize];
    for(uint64_t rp=0; rp<ncolsize; ++rp)
    {
      for(uint64_t rq=0; rq<ncolsize; ++rq)
      {
        if(rp==rq) x_in[rq] = 1.0;
        else x_in[rq] = 0.;
      }
      perform_op(x_in, y_out);
      for(uint64_t rq=0; rq<ncolsize; ++rq)
      {
        mat(rq, rp) = y_out[rq]/normfactor;
      }
    }
    delete [] x_in;
    delete [] y_out;
    return mat;
  }
  MatVecProd(SPIN_TYPE nstate): normfactor(1.0), nstate(nstate), ncolsize(nstate)
  {
    OpCounter::reset();
  }
  inline double getnormfactor() const {return normfactor;}
protected:
  double normfactor;
  SPIN_TYPE nstate, ncolsize;
};

/** 2D spin model
 */
template<typename S>
class MatVecProd2D: public MatVecProd<S>
{
public:
  MatVecProd2D(int size, SPIN_TYPE nstate): MatVecProd<S>(nstate), size(size){}
  int getLayerSize(){ return size;}
protected:
  int size;
};

/** 3D spin model
 */
template<typename S>
class MatVecProd3D: public MatVecProd<S>
{
public:
  MatVecProd3D(int size1, int size2, SPIN_TYPE nstate): MatVecProd<S>(nstate), size1(size1), size2(size2), size(size1*size2){}
  int getLayerSize(){ return size;}
protected:
  int size1, size2, size;
};

/** 2D spin model propagates single layer, perpendicular to NNN direction
 */
class XMatVecProd2D: public SymMatPhys, public MatVecProd2D<SymMatPhys::Scalar>
{
public:
  XMatVecProd2D(int size, SPIN_TYPE nstate, int boundary_cond=0);
  static uint64_t getnstate(int size){return 1ull<<size;}
protected:
};

/** 2D spin model propagates along NNN direction
 */
class YMatVecProd2D: public GenMatPhys, public MatVecProd2D<GenMatPhys::Scalar>
{
public:
  YMatVecProd2D(int size, SPIN_TYPE nstate, int boundary_cond=0);
  static uint64_t getnstate(int size){return 1ull<<(2*size);}
protected:
};


/** 2D spin model propagates single layer along diagonal direction
 */
class ZMatVecProd2D: public GenMatPhys, public MatVecProd2D<GenMatPhys::Scalar>
{
public:
  ZMatVecProd2D(int size, SPIN_TYPE nstate, int boundary_cond=0);
  virtual int moments(int size, const Scalar *ql, const Scalar *qr, double results[TMAT_NMOMENTS]) const override;
  static uint64_t getnstate(int size){return 1ull<<size;}
protected:
};

/** 3D Ising model
 */
class XMatVecProd3D: public SymMatPhys, public MatVecProd3D<SymMatPhys::Scalar>
{
public:
  XMatVecProd3D(int size, SPIN_TYPE nstate, int boundary_cond=0);
  static uint64_t getnstate(int sizeL){return 1ull<<sizeL;}
};

/** 3D spin model
 */
class YMatVecProd3D: public GenMatPhys, public MatVecProd3D<GenMatPhys::Scalar>
{
public:
  YMatVecProd3D(int size, SPIN_TYPE nstate, int boundary_cond=0);
  static uint64_t getnstate(int sizeL){return 1ull<<(2*sizeL);}
};

/** 2D x-ANNNI model
 */
class XATMatVecProd2D: public XMatVecProd2D
{
public:
  XATMatVecProd2D(const MParameter &mpara, int boundary_cond=0);
  // y_out = M * x_in
  void perform_op(const double *x_in, double *y_out) const override;
  // y_out = (P*T1) * x_in
  void permute_mult_T1(const double *x_in, double *y_out) const;
  /** generate layer vector
   */
  static Eigen::VectorXd genTvec(const MParameter &mpara, const int boundary_cond=0)
  {
    SPIN_TYPE nstate = getnstate(mpara.size), mask=nstate-1, shift1state, shift2state;
    int n1anti, n2anti, n0anti;
    // H = J1 <s1 s2> + J2 [s1 s2] - hJ1 s1 = -J <s1 s2> + kappa*J [s1 s2] + h*J*s1
    double J1 = -mpara.betaJ, J2 = mpara.betaJ*mpara.kappa, hJ = -mpara.betaJ*mpara.h;
    Eigen::VectorXd V(nstate);
    switch(boundary_cond)
    {
      // periodic boundary condition
      case(0):
        for(uint64_t nowstate = 0LL; nowstate < nstate; ++nowstate)
        {
          shift1state = Spin::offsetState(nowstate, mpara.size, 1);
          shift2state = Spin::offsetState(nowstate, mpara.size, 2);
          n1anti = Spin::nAntiSpin(nowstate, shift1state);
          n2anti = Spin::nAntiSpin(nowstate, shift2state);
          n0anti = Spin::nAntiSpin(nowstate);
          double bz1 = Spin::spinBzExpn(mpara.size, n1anti, J1);
          double bz2 = Spin::spinBzExpn(mpara.size, n2anti, J2);
          double bz0 = Spin::spinBzExpn(mpara.size, n0anti, hJ);
          V(nowstate) = exp(bz1+bz2+bz0);
        }
        break;
      // open boundary condition
      case(1):
        for(uint64_t nowstate = 0LL; nowstate < nstate; ++nowstate)
        {
          // open boundary
          shift1state = nowstate >> 1;
          shift2state = nowstate >> 2;
          n1anti = Spin::nAntiSpin(nowstate & mask>>1, shift1state);
          n2anti = Spin::nAntiSpin(nowstate & mask>>2, shift2state);
          n0anti = Spin::nAntiSpin(nowstate);
          double bz1 = Spin::spinBzExpn(mpara.size-1, n1anti, J1);
          double bz2 = Spin::spinBzExpn(mpara.size-2, n2anti, J2);
          double bz0 = Spin::spinBzExpn(mpara.size, n0anti, hJ);
          V(nowstate) = exp(bz1+bz2+bz0);
        }
        break;
      // fix first and last spin to 0
      case(5):
        for(uint64_t nowstate = 0LL; nowstate < nstate; ++nowstate)
        {
          shift1state = nowstate << 1; // extend state
          shift2state = nowstate >> 1;
          n1anti = Spin::nAntiSpin(nowstate, shift1state);
          n2anti = Spin::nAntiSpin(shift1state & mask, shift2state);
          n0anti = Spin::nAntiSpin(nowstate);
          double bz1 = Spin::spinBzExpn(mpara.size+1, n1anti, J1);
          double bz2 = Spin::spinBzExpn(mpara.size, n2anti, J2);
          double bz0 = Spin::spinBzExpn(mpara.size, n0anti, hJ);
          V(nowstate) = exp(bz1+bz2+bz0);
        }
        break;
      // fix first spin to 0 and last to 1
      case(9):
        for(uint64_t nowstate = 0LL; nowstate < nstate; ++nowstate)
        {
          shift1state = nowstate << 1 | 1; // extend state
          shift2state = nowstate >> 1;
          n1anti = Spin::nAntiSpin(nowstate, shift1state);
          n2anti = Spin::nAntiSpin(shift1state & mask, shift2state);
          n0anti = Spin::nAntiSpin(nowstate);
          double bz1 = Spin::spinBzExpn(mpara.size+1, n1anti, J1);
          double bz2 = Spin::spinBzExpn(mpara.size, n2anti, J2);
          double bz0 = Spin::spinBzExpn(mpara.size, n0anti, hJ);
          V(nowstate) = exp(bz1+bz2+bz0);
        }
        break;
      // first two spin to 0 and last two to 0
      case(13):
        for(uint64_t nowstate = 0LL; nowstate < nstate; ++nowstate)
        {
          shift1state = nowstate << 1; // extend state
          shift2state = nowstate << 2;
          n1anti = Spin::nAntiSpin(nowstate, shift1state);
          n2anti = Spin::nAntiSpin(nowstate, shift2state);
          n0anti = Spin::nAntiSpin(nowstate);
          double bz1 = Spin::spinBzExpn(mpara.size+1, n1anti, J1);
          double bz2 = Spin::spinBzExpn(mpara.size+2, n2anti, J2);
          double bz0 = Spin::spinBzExpn(mpara.size, n0anti, hJ);
          V(nowstate) = exp(bz1+bz2+bz0);
        }
        break;
        // first two spin to 0 and last two to 0
      case(17):
        for(uint64_t nowstate = 0LL; nowstate < nstate; ++nowstate)
        {
          shift1state = nowstate << 1 | 1; // extend state
          shift2state = nowstate << 2 | 3;
          n1anti = Spin::nAntiSpin(nowstate, shift1state);
          n2anti = Spin::nAntiSpin(nowstate, shift2state);
          n0anti = Spin::nAntiSpin(nowstate);
          double bz1 = Spin::spinBzExpn(mpara.size+1, n1anti, J1);
          double bz2 = Spin::spinBzExpn(mpara.size+2, n2anti, J2);
          double bz0 = Spin::spinBzExpn(mpara.size, n0anti, hJ);
          V(nowstate) = exp(bz1+bz2+bz0);
        }
        break;
    }
    return V;
  }
private:
  const MParameter mpara;
  const int boundary_cond;
  double bzwj1[2], bzwj2[2];
  Eigen::VectorXd Tvec;
};


/** 2D y-ANNNI transpose
 */
class YATMatVecProd2D_T: public YMatVecProd2D
{
public:
  YATMatVecProd2D_T(const MParameter &mpara);
  // y_out = M * x_in
  void perform_op(const double *x_in, double *y_out) const override;
  // y_out = (P*T1) * x_in
  void permute_mult_T1(const double *x_in, double *y_out) const;
  /** generate layer vector
   */
  static Eigen::VectorXd genTvec(const MParameter &mpara)
  {
    // H = J1 <s1 s2> + J2 [s1 s2] - hJ1 s1 = -J <s1 s2> + kappa*J [s1 s2] + h*J*s1
    Eigen::MatrixXd Tmat2D = genXATmat2D({mpara.betaJ, 0., 0., mpara.h, mpara.size});
    return Eigen::VectorXd(Eigen::Map<Eigen::VectorXd>(Tmat2D.data(), Tmat2D.cols()*Tmat2D.rows()));
  }
private:
  const MParameter mpara;
  double bzwj1[2], bzwj2[2];
  Eigen::VectorXd Tvec;
};

/** 2D y-ANNNI model
 * I I
 *   I I
 */
class YATMatVecProd2D: public YMatVecProd2D
{
public:
  YATMatVecProd2D(const MParameter &mpara);
  // y_out = M * x_in
  void perform_op(const double *x_in, double *y_out) const override;
  // y_out = (P*T1) * x_in
  /** get transpose matrix vector multiplication instance
   */
  using MatVecProd_T = YATMatVecProd2D_T;
  MatVecProd_T get_t_instance() const {return YATMatVecProd2D_T(mpara);}
  void permute_mult_T1(const double *x_in, double *y_out) const;
  /** generate layer vector
   */
  static Eigen::VectorXd genTvec(const MParameter &mpara)
  {
    // H = J1 <s1 s2> + J2 [s1 s2] - hJ1 s1 = -J <s1 s2> + kappa*J [s1 s2] + h*J*s1
    Eigen::MatrixXd Tmat2D = genXATmat2D({mpara.betaJ, 0., 0., mpara.h, mpara.size});
    return Eigen::VectorXd(Eigen::Map<Eigen::VectorXd>(Tmat2D.data(), Tmat2D.cols()*Tmat2D.rows()));
  }
private:
  const MParameter mpara;
  double bzwj1[2], bzwj2[2];
  Eigen::VectorXd Tvec;
};

/** 2D y2-ANNNI model
 * I I
 *     I I
 */
class Y2ATMatVecProd2D: public YMatVecProd2D
{
public:
  Y2ATMatVecProd2D(const MParameter &mpara);
  
  // TODO: implement transpose
  using MatVecProd_T = Y2ATMatVecProd2D;
  MatVecProd_T get_t_instance() const {return *this;}
  
  // y_out = M * x_in
  void perform_op(const double *x_in, double *y_out) const override;
  // y_out = (P*T1) * x_in
  void permute_mult_T1(const double *x_in, double *y_out) const;
  virtual int moments(int size, const Scalar *ql, const Scalar *qr, double results[TMAT_NMOMENTS]) const override;
  /** overide getnstate
   */
  static uint64_t getnstate(int sizeL){return 1ull<<(2*sizeL);}
  /** generate layer vector
   */
  static Eigen::VectorXd genTvec(const MParameter &mpara)
  {
    uint64_t nstate = getnstate(mpara.size), mask = (1ull<<mpara.size)-1;
    // H = J1 <s1 s2> + J2 [s1 s2] - hJ1 s1 = -J <s1 s2> + kappa*J [s1 s2] + h*J*s1
    double J1 = -mpara.betaJ, hJ = -mpara.betaJ*mpara.h;
    Eigen::VectorXd V(nstate);
    for(uint64_t nowstate = 0LL; nowstate < nstate; ++nowstate)
    {
      uint64_t leftlayer = nowstate&mask,
      rightlayer = nowstate>>mpara.size,
      leftshift = Spin::offsetState(leftlayer, mpara.size, 1),
      rightshift = Spin::offsetState(rightlayer, mpara.size, 1);
      int n1anti = Spin::nAntiSpin(leftlayer, leftshift),
      n2anti = Spin::nAntiSpin(rightlayer, rightshift),
      n0anti = Spin::nAntiSpin(leftlayer, rightlayer);
      int factorA = 2*(n1anti+n2anti+n0anti)-3*mpara.size,
      factorB = 2*(Spin::nAntiSpin(nowstate)-mpara.size);
      V(nowstate) = exp(J1*factorA + hJ*factorB);
    }
    return V;
  }
private:
  const MParameter mpara;
  double bzwj1[2], bzwj2[2];
  Eigen::VectorXd Tvec;
};

/** 2D x-DNNI model
 */
class XDTMatVecProd2D: public XMatVecProd2D
{
public:
  XDTMatVecProd2D(const MParameter &mpara, int boundary_cond=0);
  // y_out = M * x_in
  void perform_op(const double *x_in, double *y_out) const override;
  // y_out = (P*T1) * x_in
  void permute_mult_T1(const double *x_in, double *y_out) const;
  // the last run of (P*T1)*x_in needs information for the auxiliary bit of original spin configuration of a1
  void permute_mult_T1_last(const double *x_in, double *y_out, SPIN_TYPE lastbit) const;
  /** generate layer vector
   */
  static Eigen::VectorXd genTvec(const MParameter &mpara, const int boundary_cond=0)
  {
    uint64_t nstate = getnstate(mpara.size);
    // H = J1 <s1 s2> + J2 [s1 s2] - hJ1 s1 = -J <s1 s2> + kappa*J [s1 s2] + h*J*s1
    double J1 = -mpara.betaJ, hJ = -mpara.betaJ*mpara.h;
    Eigen::VectorXd V(nstate);
    switch(boundary_cond)
    {
      case(0): // periodic boundary condition
      {
#pragma omp parallel for
        for(uint64_t nowstate = 0LL; nowstate < nstate; ++nowstate)
        {
          uint64_t shift1state = Spin::offsetState(nowstate, mpara.size, 1);
          int n1anti = Spin::nAntiSpin(nowstate, shift1state);
          int n0anti = Spin::nAntiSpin(nowstate);
          double bz1 = Spin::spinBzExpn(mpara.size, n1anti, J1);
          double bz0 = Spin::spinBzExpn(mpara.size, n0anti, hJ);
          V(nowstate) = exp(bz1+bz0);
        }
      }
      break;
      // open boundary condition
      case(1):
      {
        SPIN_TYPE maskb = (nstate-1)>>1;
#pragma omp parallel for
        for(uint64_t nowstate = 0LL; nowstate < nstate; ++nowstate)
        {
          int n1anti = Spin::nAntiSpin(nowstate & maskb, nowstate >> 1);
          int n0anti = Spin::nAntiSpin(nowstate);
          double bz1 = Spin::spinBzExpn(mpara.size-1, n1anti, J1);
          double bz0 = Spin::spinBzExpn(mpara.size, n0anti, hJ);
          V(nowstate) = exp(bz1+bz0);
        }
      }
      break;
      // fix first and last spin to 0
      case(5):
      {
#pragma omp parallel for
        for(uint64_t nowstate = 0LL; nowstate < nstate; ++nowstate)
        {
          int n1anti = Spin::nAntiSpin(nowstate << 1, nowstate);
          int n0anti = Spin::nAntiSpin(nowstate);
          double bz1 = Spin::spinBzExpn(mpara.size+1, n1anti, J1);
          double bz0 = Spin::spinBzExpn(mpara.size, n0anti, hJ);
          V(nowstate) = exp(bz1+bz0);
        }
      }
      break;
      // fix first spin to 0 and last to 1
      case(9):
      {
#pragma omp parallel for
        for(uint64_t nowstate = 0LL; nowstate < nstate; ++nowstate)
        {
          int n1anti = Spin::nAntiSpin(nowstate << 1 | 1, nowstate);
          int n0anti = Spin::nAntiSpin(nowstate);
          double bz1 = Spin::spinBzExpn(mpara.size+1, n1anti, J1);
          double bz0 = Spin::spinBzExpn(mpara.size, n0anti, hJ);
          V(nowstate) = exp(bz1+bz0);
        }
      }
      break;
    }
    return V;
  }
private:
  const MParameter mpara;
  const int boundary_cond;
  double bzwj1[2], bzwj2[2];
  Eigen::VectorXd Tvec;
};

/** 2D z-BNNNI model
 */
class ZBTMatVecProd2D;

// transpose matrix vector multiplication T^T x
class ZBTMatVecProd2D_T: public ZMatVecProd2D
{
public:
  ZBTMatVecProd2D_T(const ZBTMatVecProd2D &tmathandle);
  // y_out = M * x_in
  void perform_op(const double *x_in, double *y_out) const override;
private:
  const ZBTMatVecProd2D &tmathandle;
};

// matrix-vector multiplication
class ZBTMatVecProd2D: public ZMatVecProd2D
{
public:
  ZBTMatVecProd2D(const MParameter &mpara);
  // y_out = M * x_in
  void perform_op(const double *x_in, double *y_out) const override;
  using MatVecProd_T = ZBTMatVecProd2D_T;
  MatVecProd_T get_t_instance() const {return *this;}
  
  static Eigen::VectorXd genTvec(const MParameter &mpara)
  {
    uint64_t nstate = getnstate(mpara.size);
    // H = J1 <s1 s2> + J2 [s1 s2] + J3 {s1 s2} - hJ1 s1 = -J <s1 s2> + kappa1*J [s1 s2] + kappa2*J {s1 s2} + h*J*s1
    double J1 = -mpara.betaJ, J3 = mpara.kappa2*mpara.betaJ, hJ = -mpara.betaJ*mpara.h;
    Eigen::VectorXd V(nstate);
    for(uint64_t nowstate = 0LL; nowstate < nstate; ++nowstate)
    {
      uint64_t shift1state = Spin::offsetState(nowstate, mpara.size, 1);
      uint64_t shift2state = Spin::offsetState(nowstate, mpara.size, 2);
      int n1anti = Spin::nAntiSpin(nowstate, shift1state);
      int n2anti = Spin::nAntiSpin(nowstate, shift2state);
      int n0anti = Spin::nAntiSpin(nowstate);
      double bz1 = Spin::spinBzExpn(mpara.size, n1anti, J1);
      double bz2 = Spin::spinBzExpn(mpara.size, n2anti, J3);
      double bz0 = Spin::spinBzExpn(mpara.size, n0anti, hJ);
      V(nowstate) = exp(bz1+bz2+bz0);
    }
    return V;
  }
private:
  // optdirection 0: y_out = M * x_in, 1 : y_out = M^T * x_in
  void perform_op_internal(const double *x_in, double *y_out, int optdirection) const;
  // y_out = (P*T1) * x_in
  void permute_mult_T1(const double *x_in, double *y_out, int optdirection) const;
  // the last run of (P*T1)*x_in needs information for the auxiliary bit of original spin configuration of a1
  void permute_mult_T1_last(const double *x_in, double *y_out, SPIN_TYPE lastbit, int optdirection) const;
  const MParameter mpara;
  double bzwj1[2], bzwj2[2], bzwj3[2];
  Eigen::VectorXd Tvec;
  friend class ZBTMatVecProd2D_T;
};


/** count occupancy
 */
class equivbits{
public:
  equivbits(int size)
  {
    uint64_t nchunk;
    if(size<=6)
    {
      nchunk = 1;
    }
    else
    {
      nchunk = 1ull<<(size-6);
    }
    buffer = new uint64_t[nchunk];
    memset(buffer, 0, sizeof(uint64_t)*nchunk);
   
  }
  
  bool getbit(uint64_t state)
  {
    return buffer[state>>6] & 1ull<<(state&63);
  }
  
  void setbit(uint64_t state)
  {
    buffer[state>>6] |= 1ull<<(state&63);
  }
  
  ~equivbits()
  {
    delete [] buffer;
  }
private:
  uint64_t* buffer;
};



}
#endif /* xtmat_hpp */
