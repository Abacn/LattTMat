//
//  utility.hpp
//  TMat
//
//  Created by Yi Hu on 1/20/20.
//

#ifndef utility_hpp
#define utility_hpp

#include "commondef.hpp"
#include <cmath>

#ifdef _MSC_VER
#  include <intrin.h>
#  define __builtin_popcountll(A) _mm_popcnt_u64(A)
#endif

#ifndef __POPCNT__
#pragma warning( POPCNT instruction not enabled)
#endif

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wparentheses"

namespace Spin
{
  
/**
 * inline double int small pow
 * @param a array of size 2, (base, 1/base)
 * @param b the exponent
 * @return a[0]^(-b)
 * complexity in O(log b). To avoid using exp(...) which is slower.
 */
inline double powdi(const double *a, int b)
{
  double result = 1.;
  if(b<0)
  {
    b = -b;
  }
  else
  {
    ++a;
  }
  double dtmp = *a;
  while(b)
  {
    if(b&1) result *= dtmp;
    dtmp *= dtmp;
    b >>= 1;
  }
  return result;
}
/**
 * Count number of negative spins of a state
 */
inline int nAntiSpin(SPIN_TYPE state)
{
  return __builtin_popcountll(state);
}

/**
 * Count number of opposite direction spins between two state
 */
inline int nAntiSpin(SPIN_TYPE state1, SPIN_TYPE state2)
{
  return __builtin_popcountll(state1 ^ state2);
}

/**
 * integer factor of a state in the exponential term of Boltzmann weight
 */
inline int spinBzexp(SPIN_TYPE state, int nspin)
{
  return 2*__builtin_popcountll(state) - nspin;
}

/**
 * integer factor of two state interaction in the exponential term of Boltzmann weight
 */
inline int spinBzexp(SPIN_TYPE state1, SPIN_TYPE state2, int nspin)
{
  return 2*__builtin_popcountll(state1 ^ state2) - nspin;
}

/**
 * @return a state with spin right shifted by neidist
 */
inline SPIN_TYPE offsetState(SPIN_TYPE state, int nspin, int neidist)
{
  return state >> neidist  | (state & (1ULL<<neidist)-1) << (nspin-neidist);
}

/**
 * @return a state with spin ledt shifted by neidist
 */
inline SPIN_TYPE lffsetState(SPIN_TYPE state, int nspin, int neidist)
{
  return ((state << neidist) & ((1ULL<<nspin)-1)) | (state>>(nspin-neidist));
}

/**
 * Find number of offset that state nowstate recover targetstate
 * if not matched, return nspin
 */
inline int findcyclic(SPIN_TYPE nowstate, SPIN_TYPE targetstate, int nspin)
{
  SPIN_TYPE nextstate=nowstate;
  int rp;
  for(rp=0; rp<nspin; ++rp)
  {
    if(nextstate == targetstate) break;
    nextstate = nextstate>>1 | (nextstate&1)<<(nspin-1);
  }
  return rp;
}
/**
 * @return the exponent of the Boltzmann factor (-J \sum s1 s2)
 * @param nspin number of spin pairs
 * @param nantispin number of spin pairs has opposite sign
 * @param J interaction strength
 */
inline double spinBzExpn(int nspin, int nantispin, double J)
{
  return J*(2.0*nantispin - nspin);
}

/**
 * Calculate momemts of magnitization given a eigenvector
 * Transfer matrix propagates in one layer (XTMat)
 * @param nspin number of spins
 * @param q eigenvector
 * @return success-failure status
 */
int XTmoments(int nspin, const double *q, double results[4]);

/**
 * Reverse bits: count bits backwards
 * @param nowstate now state
 * @param nbit number of bits
 * @return reversed state
 */
uint64_t bitreversal(uint64_t nowstate, int nbit);

/**
 * Calculate number of spin turnover on layer
 */
int nnode(SPIN_TYPE nowstate, int nspin, int boundary_cond);

/**
 * n choose k
 */
inline uint64_t nchoosek(uint64_t n, uint64_t k)
{
  uint64_t x=1, result=1;
  if(k*2>n) k=n-k;
  while(x <= k)
  {
    result *= n--;
    result /= x++;
  }
  return result;
}

// 3D utilities
/**
 * Right offset spin
 */
SPIN_TYPE rffsetState_3D(SPIN_TYPE state, int size1, int size2, int neidist);

/**
 * Down offset spin
 */
SPIN_TYPE dffsetState_3D(SPIN_TYPE state, int size1, int size2, int neidist);

/**
 * Rotate spin 90 degree
 */
SPIN_TYPE rotateState_3D(SPIN_TYPE state, int size);

/**
 * reflect spin layer
 */
SPIN_TYPE reflectState_3D(SPIN_TYPE state, int size1, int size2, int neidist);

}

// global counter
class OpCounter
{
public:
  static void reset(){ counts = 0;}
  static void inc(){++counts;}
  static int get(){return counts;}
private:
  static int counts;
};

#pragma GCC diagnostic pop

#endif /* utility_hpp */
