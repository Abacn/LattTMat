//
//  commondef.hpp
//  ANNNITMat
//
//  Created by Yi Hu on 2/11/20.


#ifndef commondef_h
#define commondef_h

#include <stdint.h>
#include <complex>

typedef uint64_t SPIN_TYPE;

/** Model parameter struct
 */
typedef struct{
  double betaJ;
  double kappa;
  double kappa2;
  double h;
  int size;
} MParameter;

// number of moments calcualted
#define TMAT_NMOMENTS 4

#endif /* commondef_h */
