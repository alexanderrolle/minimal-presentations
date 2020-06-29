#ifndef SMART_REDUCTION
#define SMART_REDUCTION 1
#endif

#ifndef LAZY_MINIMIZATION
#define LAZY_MINIMIZATION 1
#endif

#ifndef CHUNK_PREPROCESSING
#define CHUNK_PREPROCESSING 1
#endif

#ifndef CLEARING
#define CLEARING 0
#endif

#define TIMERS 1

#define SWAP_GRADE 0

#include <iostream>
#include <fstream>
#include <cmath>

#if TIMERS
#include <boost/timer/timer.hpp>

boost::timer::cpu_timer overall_timer, chunk_timer,firep_timer, mingens_timer, kerbasis_timer, 
  reparam_timer,minimize_timer,
  test_timer1, test_timer2, test_timer3, test_timer4;

void initialize_timers() {
  overall_timer.start();
  overall_timer.stop();
  chunk_timer.start();
  chunk_timer.stop();
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
  std::cout << "Chunk timer: " << double(chunk_timer.elapsed().wall)/std::pow(10,9) << "      ( " << double(chunk_timer.elapsed().wall)/double(overall_timer.elapsed().wall)*100 << "% )" << std::endl;
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


#include "GradedMatrix.h"

typedef phat::GradedMatrix<phat::vector_vector> GrMat;


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

  GrMat* preGM1 = new GrMat();
  GrMat* preGM2 = new GrMat();

#if TIMERS
  firep_timer.start();
#endif
  create_matrix_from_firep(ifstr,*preGM1,*preGM2);

#if TIMERS
  firep_timer.stop();
#endif

  /*
  std::cout << "Number of threads: " << omp_get_num_threads() << std::endl;
  std::cout << "Max Number of threads: " << omp_get_max_threads() << std::endl;
  omp_set_num_threads(omp_get_max_threads());
  std::cout << "Number of threads: " << omp_get_num_threads() << std::endl;
  */


#if CHUNK_PREPROCESSING

#if TIMERS
  chunk_timer.start();
#endif

  std::cout << "Chunk preprocessing..." << std::flush;
  GrMat *GM1=new GrMat();
  GrMat *GM2=new GrMat();
  chunk_preprocessing(*preGM1,*preGM2,*GM1,*GM2);
  std::cout << "done" << std::endl;
  delete preGM1;
  delete preGM2;

#if TIMERS
  chunk_timer.stop();
#endif

  //GM1.print();
  //GM2.print();

#else
  GrMat* GM1=preGM1;
  GrMat* GM2=preGM2;
#endif

  GrMat* MG = new GrMat();

#if TIMERS
  mingens_timer.start();
#endif
  std::cout << "Min Gens..." << std::flush;
  min_gens(*GM1,*MG);
  std::cout << "done" << std::endl;

  delete GM1;

#if TIMERS
  mingens_timer.stop();
#endif

  
  //MG.print();

  
  GrMat* Ker = new GrMat();

#if TIMERS
  kerbasis_timer.start();
#endif
  std::cout << "Ker basis..." << std::flush;
  ker_basis(*GM2,*Ker,*MG);
   std::cout << "done" << std::endl;
   delete GM2;
#if TIMERS
  kerbasis_timer.stop();
#endif

  //Ker.print(false);

  GrMat* semi_min_rep = new GrMat();

#if TIMERS
  reparam_timer.start();
#endif
  std::cout << "Reparameterize..." << std::flush;
  reparameterize(*MG,*Ker,*semi_min_rep);
  std::cout << "done" << std::endl;
  delete MG;
  delete Ker;
#if TIMERS
  reparam_timer.stop();
#endif

  //semi_min_rep.print(true,false);

  GrMat *min_rep=new GrMat();
#if TIMERS
  minimize_timer.start();
#endif
  std::cout << "Minimize..." << std::flush;
  minimize(*semi_min_rep,*min_rep);
  std::cout << "done" << std::endl;
  delete semi_min_rep;
#if TIMERS
  minimize_timer.stop();
#endif

  std::cout << "Resulting minimal presentation has " << min_rep->get_num_cols() << " columns and " << min_rep->num_rows << " rows" << std::endl;
  //min_rep.print(true,false);

  if(argc>=3) {
    std::cout << "Writing to file \"" << argv[2] << "\"..." << std::flush;
    std::ofstream ofstr(argv[2]);
    print_in_rivet_format(*min_rep,ofstr);
    ofstr.close();
    std::cout << "done" << std::endl;
  }

  delete min_rep;

#if TIMERS
  print_timers();
#endif

  return 0;
  
}
