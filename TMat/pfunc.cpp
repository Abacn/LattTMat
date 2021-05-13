//
//  pfunc.cpp
//  TMat
//
//  Created by Yi Hu on 2/12/20.
//  Calculate partition function
//

#include <iostream>

#include "tmat.hpp"
#include "tmatreduce.hpp"
#include "tmat3d.hpp"
#include "pfunc.hpp"

void pfunc(const Read_Input &input)
{
  // check model specification: model, dim, axial-direction
  // 0, 1, 2, 3: w, x, y, z
  // +4: reduced matrix; +8: open boundary
  const MParameter mpara = get_MParameter(input);
  const int multieigs = input.neigs>1? 2: 0;
  if('A' == input.model) // ANNNI model
  {
    if(2 == input.dim)
    {
      switch(input.axialdir)
      {
        // y-double-layer-ANNNI
        case(0): pfunc_inner(input, TMat::Y2ATMatVecProd2D(mpara));
          break;
        // x-ANNNI
        case(1): pfunc_inner(input, TMat::XATMatVecProd2D(mpara));
          break;
        // y-ANNNI
        case(2): pfunc_inner(input, TMat::YATMatVecProd2D(mpara));
          break;
        // x-reduced ANNNI
        case(5): pfunc_inner(input, TMat::XAMMatVecProd2D(mpara, multieigs));
          break;
        // y-reduced ANNNI
        case(6): pfunc_inner(input, TMat::YAMMatVecProd2D(mpara, multieigs));
          break;
        // x-ANNNI, open boundary
        case(9): pfunc_inner(input, TMat::XATMatVecProd2D(mpara, 1));
          break;
        //x-reduced-ANNNI, open boundary
        case(13): pfunc_inner(input, TMat::XAMMatVecProd2D(mpara, multieigs|1));
          break;
        //open boundary with confined spins
        case(17): case(18): pfunc_inner(input, TMat::XATMatVecProd2D(mpara, input.axialdir&1 ? 5 : 9));
          break;
        //open boundary with two confined spins in each side
        case(19): case(20): pfunc_inner(input, TMat::XATMatVecProd2D(mpara, input.axialdir&1 ? 13 : 17));
          break;
        default:
        {
          std::cout << "Function yet implemented" << std::endl;
          assert(false);
        }
      }
    }
    else if(3 == input.dim)
    {
      if(1 == input.axialdir)
      {
        pfunc_inner(input, TMat::XAMMatVecProd3D(mpara));
      }
      else
      {
        std::cout << "Function yet implemented" << std::endl;
        assert(false);
      }
    }
  }
  else if('B' == input.model) // BNNNI model
  {
    switch(input.axialdir)
    {
      case(3): pfunc_inner(input, TMat::ZBTMatVecProd2D(mpara)); break;
      case(6): pfunc_inner(input, TMat::YBMMatVecProd2D(mpara, multieigs)); break;
      case(7): pfunc_inner(input, TMat::ZBMMatVecProd2D(mpara, multieigs)); break;
      default:
      {
        std::cout << "Function yet implemented" << std::endl;
        assert(false);
      }
    }
  }
  else if('D' == input.model) // DNNI (J1-J2) model
  {
    switch(input.axialdir)
    {
      // DNNI TMat is equivalent in x and y direction
      case(1): case(2): pfunc_inner(input, TMat::XDTMatVecProd2D(mpara)); break;
      // reduced matrix
      case(5): case(6): pfunc_inner(input, TMat::XDMMatVecProd2D(mpara, multieigs)); break;
      //open boundary with confined spins
      case(17): case(18): pfunc_inner(input, TMat::XDTMatVecProd2D(mpara, input.axialdir&1 ? 5 : 9)); break;
      default:
        // diagonal direction
        std::cout << "Function yet implemented" << std::endl;
        assert(false);
    }
  }
}

MParameter get_MParameter(const Read_Input &input)
{
  return {input.betaJ, input.kappa, input.kappa2, input.h, input.size};
}
