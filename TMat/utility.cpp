//
//  utility.cpp
//  TMat
//
//  Created by Yi Hu on 1/20/20.
//
#include <cmath>
#include <iostream>
#include "utility.hpp"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wparentheses"

namespace Spin
{

SPIN_TYPE bitreversal(SPIN_TYPE nowstate, int nbit)
{
#define SHIFT1  0x5555555555555555ull
#define SHIFT2  0x3333333333333333ull
#define SHIFT4  0x0f0f0f0f0f0f0f0full
#define SHIFT8  0x00ff00ff00ff00ffull
#define SHIFT16 0x0000ffff0000ffffull
  nowstate = ((nowstate >> 1) & SHIFT1) | ((nowstate & SHIFT1) << 1);
  nowstate = ((nowstate >> 2) & SHIFT2) | ((nowstate & SHIFT2) << 2);
  nowstate = ((nowstate >> 4) & SHIFT4) | ((nowstate & SHIFT4) << 4);
  nowstate = ((nowstate >> 8) & SHIFT8) | ((nowstate & SHIFT8) << 8);
  nowstate = ((nowstate >> 16) & SHIFT16) | ((nowstate & SHIFT16) << 16);
  nowstate = (nowstate >> 32 ) | (nowstate << 32);
  return nowstate >> (64-nbit);
}

int nnode(SPIN_TYPE nowstate, int nspin, int boundary_cond)
{
  SPIN_TYPE shiftstate = nowstate >> 1;
  if(!boundary_cond)
  {
    // periodic boundary
    shiftstate |= (nowstate&1) << (nspin-1);
  }
  else
  {
    // open boundary
    shiftstate |= (nowstate & nspin-1) << 1;
  }
  return __builtin_popcountll(nowstate ^ shiftstate);
}

/**
 * Right offset spin. shift size2 -- O(size1)
 */
SPIN_TYPE rffsetState_3D(SPIN_TYPE state, int size1, int size2, int neidist)
{
  SPIN_TYPE mask = (1ull << size2) - 1, result = 0;
  for(int rp=0; rp<size1; ++rp)
  {
    SPIN_TYPE thisstate = offsetState(state & mask, size2, neidist);
    thisstate <<= rp*size2;
    result |= thisstate;
    state >>= size2;
  }
  return result;
}

/**
 * Down offset spin -- O(1)
 */
SPIN_TYPE dffsetState_3D(SPIN_TYPE state, int size1, int size2, int neidist)
{
  SPIN_TYPE mask = (1ull << size2*neidist) - 1, result = state >> size2*neidist;
  result |= (state & mask) << (size1-neidist)*size2;
  return result;
}

/**
 * Rotate spin 90 degree -- O(size1*size2)
 */
SPIN_TYPE rotateState_3D(SPIN_TYPE state, int size)
{
  SPIN_TYPE result = 0;
  for(int rp=0; rp<size; ++rp)
  {
    for(int rq=0; rq<size; ++rq)
    {
      int loc1 = rp*size+rq;
      if(state & (1ull << loc1))
        result |= 1ull << ((rq+1)*size - (rp+1));
    }
  }
  return result;
}

/**
 * reflect spin layer -- O(size1)
 *  \param neidist 0: left-right reflect; 1: up-down reflect
 */
SPIN_TYPE reflectState_3D(SPIN_TYPE state, int size1, int size2, int neidist)
{
  SPIN_TYPE mask = (1ull << size2) - 1, result = 0;
  if(!neidist)
  {
    for(int rp=0; rp<size1; ++rp)
    {
      SPIN_TYPE thisstate = bitreversal(state & mask, size2);
      thisstate <<= rp*size2;
      result |= thisstate;
      state >>= size2;
    }
  }
  else
  {
    for(int rp=0; rp<size1; ++rp)
    {
      result <<= size2;
      result |= state & mask;
      state >>= size2;
    }
  }
  return result;
}


}

// static variable
int OpCounter::counts = 0;

#pragma GCC diagnostic pop
