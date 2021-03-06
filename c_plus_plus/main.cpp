#ifndef CHUNK_PREPROCESSING
#define CHUNK_PREPROCESSING 0
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
#define PARALLEL_FOR_LOOPS 0
#endif

#ifndef COLUMN_REPRESENTATION
#define COLUMN_REPRESENTATION vector_vector
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

// If parallelized, the counting does not work (not thread-safe)
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

boost::timer::cpu_timer overall_timer, chunk_timer,io_timer, mingens_timer, kerbasis_timer, 
  reparam_timer,minimize_timer,
  test_timer1, test_timer2, test_timer3, test_timer4,test_timer5;

void initialize_timers() {
  overall_timer.start();
  overall_timer.stop();
  chunk_timer.start();
  chunk_timer.stop();
  io_timer.start();
  io_timer.stop();
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
  std::cout << "IO timer: " << double(io_timer.elapsed().wall)/std::pow(10,9) << "      ( " << double(io_timer.elapsed().wall)/double(overall_timer.elapsed().wall)*100 << "% )" << std::endl;
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

typedef phat::GradedMatrix<phat::COLUMN_REPRESENTATION> GrMat;


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

  if(std::string(argv[1])=="-o" || std::string(argv[1])=="-op") {
    std::cerr << "Running in write-only mode" << std::endl;
    if(argc<4) {
      std::cerr << "Input file and output file needed" << std::endl;
      std::exit(1);
    }
    GrMat GM1,GM2;
    create_matrix_from_firep(argv[2],GM1,GM2);
    if(std::string(argv[1])=="-o") {
      write_matrix_to_firep(argv[3],GM1,GM2);
    } else {
      write_matrix_to_phat(argv[3],GM1);
    }
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


#if TIMERS
  io_timer.start();
#endif

  GrMat preGM1,preGM2;

  create_matrix_from_firep(argv[1],preGM1,preGM2);

#if TIMERS
  io_timer.stop();
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
  GM1.grid_scheduler=phat::Grid_scheduler(GM1);
  min_gens(GM1,MG);
#if !MIN_GENS_AND_KER_BASIS_IN_PARALLEL
  std::cout << "done, size is " << MG.num_rows << "x" << MG.get_num_cols() << std::endl;
  GM1=GrMat();
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
  GM2.grid_scheduler=phat::Grid_scheduler(GM2);
  ker_basis(GM2,Ker,MG);
#if !MIN_GENS_AND_KER_BASIS_IN_PARALLEL
  std::cout << "done, size is " << Ker.num_rows << "x" << Ker.get_num_cols() << std::endl;
  GM2=GrMat(); 

#if TIMERS
  kerbasis_timer.stop();
#endif
#endif

#if MIN_GENS_AND_KER_BASIS_IN_PARALLEL
    }
  }
  kerbasis_timer.stop();
  std::cout << "done" << std::endl;
  std::cout << "size of MG is " << MG.num_rows << "x" << MG.get_num_cols() << std::endl; 
  std::cout << "size of KER is " << Ker.num_rows << "x" << Ker.get_num_cols() << std::endl; 
  GM1=GrMat();
  GM2=GrMat();
#endif


#if TIMERS
  reparam_timer.start();
#endif
  GrMat semi_min_rep;
  std::cout << "Reparameterize..." << std::flush;
  reparameterize(MG,Ker,semi_min_rep);
  std::cout << "done" << std::endl;
  std::cout << "Resulting semi-minimal presentation has " << semi_min_rep.get_num_cols() << " columns and " << semi_min_rep.num_rows << " rows" << std::endl;
  MG=GrMat();
  Ker=GrMat();
  
#if TIMERS
  reparam_timer.stop();
#endif

  //semi_min_rep.print(true,false);


#if TIMERS
  minimize_timer.start();
#endif
  GrMat min_rep;
  std::cout << "Minimize..." << std::flush;
  minimize(semi_min_rep,min_rep);
  std::cout << "done" << std::endl;
  std::cout << "Resulting minimal presentation has " << min_rep.get_num_cols() << " columns and " << min_rep.num_rows << " rows" << std::endl;
  semi_min_rep=GrMat(); 
#if TIMERS
  minimize_timer.stop();
#endif


#if !PERTURBED
  if(argc>=3) {
#if TIMERS
    io_timer.resume();
#endif

    std::cout << "Writing to file \"" << argv[2] << "\"..." << std::flush;
    std::ofstream ofstr(argv[2]);
    print_in_rivet_format(min_rep,ofstr);
    ofstr.close();
    std::cout << "done" << std::endl;
#if TIMERS
    io_timer.stop();
#endif

  }
#endif

#if TIMERS
  print_timers();
#endif
#if !PARALLEL_FOR_LOOPS && !MIN_GENS_AND_KER_BASIS_IN_PARALLEL
  std::cout << "Total number of column additions: " << gl_no_column_additions << std::endl;
#endif
  return 0;
  
}
