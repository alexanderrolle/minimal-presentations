#define SMART_REDUCTION 1

#define LAZY_MINIMIZATION 1

#define TIMERS 1

#define SWAP_GRADE 0

#include <iostream>
#include <fstream>
#include "GradedMatrix.h"

typedef phat::GradedMatrix<phat::vector_vector> GrMat;

#if TIMERS
#include <boost/timer/timer.hpp>

boost::timer::cpu_timer overall_timer,firep_timer, mingens_timer, kerbasis_timer, 
  reparam_timer,minimize_timer,
  test_timer1, test_timer2, test_timer3, test_timer4;

void initialize_timers() {
  overall_timer.start();
  overall_timer.stop();
  firep_timer.start();
  firep_timer.stop();
  mingens_timer.start();
  mingens_timer.stop();
  kerbasis_timer.start();
  kerbasis_timer.stop();
  reparam_timer.start();
  reparam_timer.stop();
  minimize_timer.start();
  minimize_timer.stop();

  test_timer1.start();
  test_timer1.stop();
  test_timer2.start();
  test_timer2.stop();
  test_timer3.start();
  test_timer3.stop();
  test_timer4.start();
  test_timer4.stop();
}

void print_timers() {
  std::cout << "Overall timer: " << double(overall_timer.elapsed().wall)/std::pow(10,9) << std::endl;
  std::cout << "Firep timer: " << double(firep_timer.elapsed().wall)/std::pow(10,9) << "      ( " << double(firep_timer.elapsed().wall)/double(overall_timer.elapsed().wall)*100 << "% )" << std::endl;
  std::cout << "Mingens timer: " << double(mingens_timer.elapsed().wall)/std::pow(10,9) << "     ( " << double(mingens_timer.elapsed().wall)/double(overall_timer.elapsed().wall)*100 << "% )" << std::endl;
  std::cout << "Kerbasis timer: " << double(kerbasis_timer.elapsed().wall)/std::pow(10,9) << "    ( " << double(kerbasis_timer.elapsed().wall)/double(overall_timer.elapsed().wall)*100 << "% )" << std::endl;
  std::cout << "Reparam timer: " << double(reparam_timer.elapsed().wall)/std::pow(10,9) << "    ( " << double(reparam_timer.elapsed().wall)/double(overall_timer.elapsed().wall)*100 << "% )" << std::endl;
  std::cout << "Minimize timer: " << double(minimize_timer.elapsed().wall)/std::pow(10,9) << "    ( " << double(minimize_timer.elapsed().wall)/double(overall_timer.elapsed().wall)*100 << "% )" << std::endl;
  std::cout << "Test timer1: " << double(test_timer1.elapsed().wall)/std::pow(10,9) << std::endl;
  std::cout << "Test timer2: " << double(test_timer2.elapsed().wall)/std::pow(10,9) << std::endl;
  std::cout << "Test timer3: " << double(test_timer3.elapsed().wall)/std::pow(10,9) << std::endl;
  std::cout << "Test timer4: " << double(test_timer4.elapsed().wall)/std::pow(10,9) << std::endl;
}

#endif

int main(int argc, char** argv) {

  if(argc==1) {
    std::cerr << "Input file?" << std::endl;
    std::exit(1);
  }

#if TIMERS
  initialize_timers();
  overall_timer.start();
#endif

  std::ifstream ifstr(argv[1]);

  GrMat GM1,GM2;

#if TIMERS
  firep_timer.start();
#endif
  create_matrix_from_firep(ifstr,GM1,GM2);

#if TIMERS
  firep_timer.stop();
#endif

  //GM1.print();
  //GM2.print();

  GrMat MG;

#if TIMERS
  mingens_timer.start();
#endif
  std::cout << "Min Gens..." << std::flush;
  min_gens(GM1,MG);
  std::cout << "done" << std::endl;

#if TIMERS
  mingens_timer.stop();
#endif

  
  //MG.print();

  
  GrMat Ker;

#if TIMERS
  kerbasis_timer.start();
#endif
  std::cout << "Ker basis..." << std::flush;
  ker_basis(GM2,Ker);
   std::cout << "done" << std::endl;
#if TIMERS
  kerbasis_timer.stop();
#endif

  //Ker.print(false);

  GrMat semi_min_rep;

#if TIMERS
  reparam_timer.start();
#endif
  std::cout << "Reparameterize..." << std::flush;
  reparameterize(MG,Ker,semi_min_rep);
  std::cout << "done" << std::endl;
#if TIMERS
  reparam_timer.stop();
#endif

  //semi_min_rep.print(true,false);

  GrMat min_rep;
#if TIMERS
  minimize_timer.start();
#endif
  std::cout << "Minimize..." << std::flush;
  minimize(semi_min_rep,min_rep);
  std::cout << "done" << std::endl;
#if TIMERS
  minimize_timer.stop();
#endif

  std::cout << "Resulting minimal presentation has " << min_rep.get_num_cols() << " columns and " << min_rep.num_rows << " rows" << std::endl;
  //min_rep.print(true,false);

#if TIMERS
  print_timers();
#endif

  return 0;
  
}
