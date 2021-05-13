//
//  utility.cpp
//  DumpConfig
//
//  Created by Yi Hu on 8/5/20.
//

#include <random>
#include <ctime>

#include "utilityb.hpp"

namespace{
  std::mt19937 rg;
  std::uniform_real_distribution <> rrand;
};

namespace MyUtility{
void setseed(void)
{
  rg.seed((unsigned int)time(NULL));
}

void setseed(unsigned int seed)
{
  rg.seed(seed);
}

double rand()
{
  return rrand(rg);
}
};
