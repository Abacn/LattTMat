//
//  read_input.hpp
//  TMat
//
//  Created by Yi Hu on 2/13/20.
//

#ifndef read_input_hpp
#define read_input_hpp

#include <string>

/** input parameters
 */
class Read_Input{
public:
  int read(const char *inputf);
  char model;
  int dim;
  int axialdir;
  int neigs;
  int savevec[2];
  std::string outprefix;
  int size;
  double betaJ;
  double kappa, kappa2;
  double h;
};


#endif /* read_input_hpp */
