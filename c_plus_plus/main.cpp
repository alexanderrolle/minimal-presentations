#ifndef CHUNK_PREPROCESSING
#define CHUNK_PREPROCESSING 1
#endif

#ifndef SMART_REDUCTION
#define SMART_REDUCTION 1
#endif

#ifndef SPARSE_GRID_TRAVERSAL
#define SPARSE_GRID_TRAVERSAL 1
#endif

#ifndef CLEARING
#define CLEARING 0
#endif

#ifndef MIN_GENS_AND_KER_BASIS_IN_PARALLEL
#define MIN_GENS_AND_KER_BASIS_IN_PARALLEL 0
#endif


#ifndef LAZY_MINIMIZATION
#define LAZY_MINIMIZATION 0
#endif


#ifndef PARALLEL_FOR_LOOPS
#define PARALLEL_FOR_LOOPS 1
#endif

// The reading from input uses strtok_r as a thread-save
// version of strtok. If this function is not enabled,
// that part of code is not parallelized, which results
// in a slight performance penalty
//
// Disabled, because it causes problems on server... 
#ifndef STRINGTOK_R_AVAILABLE
#define STRINGTOK_R_AVAILABLE 0
#endif

#if CLEARING && MIN_GENS_AND_KER_BASIS_IN_PARALLEL
#error CLEARING and  MIN_GENS_AND_KER_BASIS_IN_PARALLEL cannot both be enabled
#endif

#if !SMART_REDUCTION && SPARSE_GRID_TRAVERSAL
#error SPARSE_GRID_TRAVERSAL requires SMART_REDUCTION
#endif

int gl_no_column_additions=0;

#define TIMERS 1

#ifndef SWAP_GRADE
#define SWAP_GRADE 0
#endif

#ifndef PERTURB
#define PERTURB 0
#endif

#if PERTURB
#include<random>
#endif

#include <iostream>
#include <fstream>
#include <cmath>

#if TIMERS
#include <boost/timer/timer.hpp>

boost::timer::cpu_timer overall_timer, chunk_timer,firep_timer, mingens_timer, kerbasis_timer, 
  reparam_timer,minimize_timer,
  test_timer1, test_timer2, test_timer3, test_timer4,test_timer5;

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
  test_timer5.start();
  test_timer5.stop();
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
  std::cout << "Test timer5: " << double(test_timer5.elapsed().wall)/std::pow(10,9) << std::endl;
}

#endif


#include "GradedMatrix.h"

typedef phat::GradedMatrix<phat::vector_vector> GrMat;


int main(int argc, char** argv) {

  if(argc==1) {
    std::cerr << "Input file?" << std::endl;
    std::exit(1);
  }

  if(std::string(argv[1])=="-c") {
    std::cerr << "Running in check only mode" << std::endl;
    if(argc<3) {
      std::cerr << "Input file?" << std::endl;
      std::exit(1);
    }
    GrMat GM1,GM2;
    create_matrix_from_firep(argv[2],GM1,GM2);
    check_grade_sanity(GM1);
    check_boundaries(GM1,"first");
    check_boundaries(GM2,"second");
    std::cerr << "Input is valid, exiting..." << std::endl;
    std::exit(0);
  }

#if PERTURBED
  std::cout << "Warning: PERTURBED is enabled, writing to output file is supressed" << std::endl;
#endif

#if TIMERS
  initialize_timers();
  overall_timer.start();
#endif

#if PARALLEL_FOR_LOOPS
  std::cout << "Execution parralized, max Number of threads: " << omp_get_max_threads() << std::endl;
#endif

#if MIN_GENS_AND_KER_BASIS_IN_PARALLEL
  if(omp_get_max_threads()<2) {
    std::cerr << "Option MIN_GENS_AND_KER_BASIS_IN_PARALLEL not possible with one thread" << std::endl;
    std::exit(1);
  }
#endif



  GrMat preGM1,preGM2;

#if TIMERS
  firep_timer.start();
#endif
  create_matrix_from_firep(argv[1],preGM1,preGM2);

#if TIMERS
  firep_timer.stop();
#endif



#if CHUNK_PREPROCESSING

#if TIMERS
  chunk_timer.start();
#endif

  std::cout << "Chunk preprocessing..." << std::flush;
  GrMat GM1;
  GrMat GM2;
  //test_timer1.start();
  chunk_preprocessing(preGM1,preGM2,GM1,GM2);
  //test_timer1.stop();
  std::cout << "done" << std::endl;

  preGM1=GrMat();
  preGM2=GrMat();
  

#if TIMERS
  chunk_timer.stop();
#endif

  //GM1.print();
  //GM2.print();

#else
  GrMat& GM1=preGM1;
  GrMat& GM2=preGM2;
#endif




  GrMat MG,Ker;

#if MIN_GENS_AND_KER_BASIS_IN_PARALLEL
  std::cout << "Min Gens and ker_basis in parallel..." << std::flush;
  kerbasis_timer.start();
#pragma omp parallel num_threads(2)
  {
    
    int tid = omp_get_thread_num();
    if(tid==0) {
#endif

#if !MIN_GENS_AND_KER_BASIS_IN_PARALLEL
#if TIMERS
  mingens_timer.start();
#endif
  std::cout << "Min Gens..." << std::flush;
#endif
  min_gens(GM1,MG);
#if !MIN_GENS_AND_KER_BASIS_IN_PARALLEL
  std::cout << "done" << std::endl;
#if TIMERS
  mingens_timer.stop();
#endif
#endif

#if MIN_GENS_AND_KER_BASIS_IN_PARALLEL
    }
    if(tid==1) {
#endif

#if !MIN_GENS_AND_KER_BASIS_IN_PARALLEL
#if TIMERS
  kerbasis_timer.start();
#endif

  std::cout << "Ker basis..." << std::flush;
#endif
  ker_basis(GM2,Ker,MG);
#if !MIN_GENS_AND_KER_BASIS_IN_PARALLEL
  std::cout << "done" << std::endl;
#if TIMERS
  kerbasis_timer.stop();
#endif
#endif

#if MIN_GENS_AND_KER_BASIS_IN_PARALLEL
    }
  }
  kerbasis_timer.stop();
  std::cout << "done" << std::endl;
#endif

  GM1=GrMat();
  GM2=GrMat();

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

  std::cout << "Resulting semi-minimal presentation has " << semi_min_rep.get_num_cols() << " columns and " << semi_min_rep.num_rows << " rows" << std::endl;

  MG=GrMat();

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

#if !PERTURBED
  if(argc>=3) {
    std::cout << "Writing to file \"" << argv[2] << "\"..." << std::flush;
    std::ofstream ofstr(argv[2]);
    print_in_rivet_format(min_rep,ofstr);
    ofstr.close();
    std::cout << "done" << std::endl;
  }
#endif

#if TIMERS
  print_timers();
#endif
  std::cout << "Total number of column additions: " << gl_no_column_additions << std::endl;

  return 0;
  
}
