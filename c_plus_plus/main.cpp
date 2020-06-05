#include <iostream>
#include <fstream>
#include "GradedMatrix.h"

int main(int argc, char** argv) {

  std::ifstream ifstr("firep.txt");

  phat::GradedMatrix<> GM1,GM2;

  create_matrix_from_firep(ifstr,GM1,GM2);
  
  GM1.print();
  GM2.print();

  phat::GradedMatrix<> MG;

  
  min_gens(GM1,MG);

  std::cout << "After" << std::endl;
  GM1.print();

  std::cout << "Min Gens" << std::endl;
  MG.print();

  
  phat::GradedMatrix<> Ker;
  ker_basis(GM2,Ker);

  std::cout << "Ker basis" << std::endl;
  Ker.print(false);

  phat::GradedMatrix<> semi_min_rep;
  reparameterize(MG,Ker,semi_min_rep);

  std::cout << "Semi-min" << std::endl;
  semi_min_rep.print();

  std::cout << "Minimize" << std::endl;
  minimize(semi_min_rep);

  return 0;
  
}
