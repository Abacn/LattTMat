//
//  read_input.cpp
//  TMat
//
//  Created by Yi Hu on 2/13/20.
//

#include <iostream>
#include <sstream>
#include <fstream>

#include "read_input.hpp"

int Read_Input::read(const char *inputf)
{
  int error = 0, size1 = 0, size2 = 0;
  std::ifstream infile(inputf);
  if(!infile)
  {
    std::cout << "Can't open " << inputf << " for input." << std::endl;
    error = 2;
    return error;
  }
  else
  {
    std::cout << "Reading input from file " << inputf << std::endl;
  }
  char buf[300],c;
  infile.get(buf,300,'='); infile.get(c); infile >> model;
  infile.get(buf,300,'='); infile.get(c); infile >> dim;
  infile.get(buf,300,'='); infile.get(c); infile >> axialdir;
  infile.get(buf,300,'='); infile.get(c); infile >> neigs;
  infile.get(buf,300,'='); infile.get(c); infile >> savevec[0] >> savevec[1];
  infile.get(buf,300,'='); infile.get(c); infile >> outprefix;
  infile.get(buf,300,'='); infile.get(c);
  if(2==dim)
    infile >> size;
  else if(3==dim)
  {
    infile >> size1 >> size2;
    size = size1 << 16 | size2; // embedded
  }
  else
  {
    std::cout << "Invalid dim " << dim << std::endl;
    error = 4;
    return error;
  }
  infile.get(buf,300,'='); infile.get(c); infile >> betaJ;
  // deal with kappa line
  infile.get(buf,300,'=');
  std::string kappaline;
  std::getline(infile, kappaline);
  std::istringstream kappass(kappaline);
  kappass.get(c); kappass >> kappa; kappass >> kappa2;
  if(!kappass) // error occurred
  {
    kappa2 = 0.;
  }
  infile.get(buf,300,'='); infile.get(c); infile >> h;
  if(infile.fail())
  {
    std::cout << "Error reading input file " << inputf << std::endl;
    error = 3;
  }
  
  std::cout << "     MODEL : " << dim << "D_" << model << axialdir << "\n";
  std::cout << "Save & flg : " << outprefix << " " << savevec[0] << " " << savevec[1] << "\n";
  if(2==dim)
    std::cout << "      size : " << size << "\n";
  else
    std::cout << "      size : " << size1 << " " << size2 << "\n";
  std::cout << "     betaJ : " << betaJ << "\n";
  std::cout << "     kappa : " << kappa << std::endl;
  if('B' == model)
  {
    std::cout << "    kappa2 : " << kappa2 << std::endl;
  }
  std::cout << "         h : " << h << std::endl;
  // check valid
  if(model != 'A' && model != 'B' && model != 'D')
  {
    std::cout << "Invalid model " << model << std::endl;
    error = 4;
  }
  else if(size<4)
  {
    std::cout << "Invalid size (<4) " << size << std::endl;
    error = 4;
  }
  return error;
}
