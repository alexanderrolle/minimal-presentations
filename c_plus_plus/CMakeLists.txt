project( Code_ )

cmake_minimum_required(VERSION 2.8.10)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_VERBOSE_MAKEFILE ON)

set(CMAKE_BUILD_TYPE Release)

#set(CMAKE_CXX_COMPILER g++-8)

FIND_PACKAGE(OpenMP)

include_directories(phat_mod/include)

find_package(Boost REQUIRED timer)

set (CMAKE_CXX_FLAGS "-fopenmp")

add_executable(main "main.cpp")
target_link_libraries(main "-lboost_timer -lgmp -lgmpxx")


add_executable(main_rivet "main.cpp")
target_compile_options(main_rivet PUBLIC -DCHUNK_PREPROCESSING=0 -DSMART_REDUCTION=0 -DSPARSE_GRID_TRAVERSAL=0 -DCLEARING=0 -DMIN_GENS_AND_KER_BASIS_IN_PARALLEL=0 -DLAZY_MINIMIZATION=0 -DPARALLEL_FOR_LOOPS=0 -DSWAP_GRADE=0 -DPERTURB=0)
target_link_libraries(main_rivet "-lboost_timer -lgmp -lgmpxx")

add_executable(main_rivet_swapped "main.cpp")
target_compile_options(main_rivet_swapped PUBLIC -DCHUNK_PREPROCESSING=0 -DSMART_REDUCTION=0 -DSPARSE_GRID_TRAVERSAL=0 -DCLEARING=0 -DMIN_GENS_AND_KER_BASIS_IN_PARALLEL=0 -DLAZY_MINIMIZATION=0 -DPARALLEL_FOR_LOOPS=0 -DSWAP_GRADE=1 -DPERTURB=0)
target_link_libraries(main_rivet_swapped "-lboost_timer -lgmp -lgmpxx")

add_executable(main_lazy "main.cpp")
target_compile_options(main_lazy PUBLIC -DCHUNK_PREPROCESSING=0 -DSMART_REDUCTION=0 -DSPARSE_GRID_TRAVERSAL=0 -DCLEARING=0 -DMIN_GENS_AND_KER_BASIS_IN_PARALLEL=0 -DLAZY_MINIMIZATION=1 -DPARALLEL_FOR_LOOPS=0 -DSWAP_GRADE=0 -DPERTURB=0)
target_link_libraries(main_lazy "-lboost_timer -lgmp -lgmpxx")

add_executable(main_smart_sparse "main.cpp")
target_compile_options(main_smart_sparse PUBLIC -DCHUNK_PREPROCESSING=0 -DSMART_REDUCTION=1 -DSPARSE_GRID_TRAVERSAL=1 -DCLEARING=0 -DMIN_GENS_AND_KER_BASIS_IN_PARALLEL=0 -DLAZY_MINIMIZATION=0 -DPARALLEL_FOR_LOOPS=0 -DSWAP_GRADE=0 -DPERTURB=0)
target_link_libraries(main_smart_sparse "-lboost_timer -lgmp -lgmpxx")

add_executable(main_lazy_smart "main.cpp")
target_compile_options(main_lazy_smart PUBLIC -DCHUNK_PREPROCESSING=0 -DSMART_REDUCTION=1 -DSPARSE_GRID_TRAVERSAL=0 -DCLEARING=0 -DMIN_GENS_AND_KER_BASIS_IN_PARALLEL=0 -DLAZY_MINIMIZATION=1 -DPARALLEL_FOR_LOOPS=0 -DSWAP_GRADE=0 -DPERTURB=0)
target_link_libraries(main_lazy_smart "-lboost_timer -lgmp -lgmpxx")

add_executable(main_lazy_smart_sparse "main.cpp")
target_compile_options(main_lazy_smart_sparse PUBLIC -DCHUNK_PREPROCESSING=0 -DSMART_REDUCTION=1 -DSPARSE_GRID_TRAVERSAL=1 -DCLEARING=0 -DMIN_GENS_AND_KER_BASIS_IN_PARALLEL=0 -DLAZY_MINIMIZATION=1 -DPARALLEL_FOR_LOOPS=0 -DSWAP_GRADE=0 -DPERTURB=0)
target_link_libraries(main_lazy_smart_sparse "-lboost_timer -lgmp -lgmpxx")

add_executable(main_lazy_smart_sparse_chunk "main.cpp")
target_compile_options(main_lazy_smart_sparse_chunk PUBLIC -DCHUNK_PREPROCESSING=1 -DSMART_REDUCTION=1 -DSPARSE_GRID_TRAVERSAL=1 -DCLEARING=0 -DMIN_GENS_AND_KER_BASIS_IN_PARALLEL=0 -DLAZY_MINIMIZATION=1 -DPARALLEL_FOR_LOOPS=0 -DSWAP_GRADE=0 -DPERTURB=0)
target_link_libraries(main_lazy_smart_sparse_chunk "-lboost_timer -lgmp -lgmpxx")

add_executable(main_lazy_smart_sparse_chunk_parfor "main.cpp")
target_compile_options(main_lazy_smart_sparse_chunk_parfor PUBLIC -DCHUNK_PREPROCESSING=1 -DSMART_REDUCTION=1 -DSPARSE_GRID_TRAVERSAL=1 -DCLEARING=0 -DMIN_GENS_AND_KER_BASIS_IN_PARALLEL=0 -DLAZY_MINIMIZATION=1 -DPARALLEL_FOR_LOOPS=1 -DSWAP_GRADE=0 -DPERTURB=0)
target_link_libraries(main_lazy_smart_sparse_chunk_parfor "-lboost_timer -lgmp -lgmpxx")

add_executable(main_lazy_smart_sparse_chunk_parfor_clearing "main.cpp")
target_compile_options(main_lazy_smart_sparse_chunk_parfor_clearing PUBLIC -DCHUNK_PREPROCESSING=1 -DSMART_REDUCTION=1 -DSPARSE_GRID_TRAVERSAL=1 -DCLEARING=1 -DMIN_GENS_AND_KER_BASIS_IN_PARALLEL=0 -DLAZY_MINIMIZATION=1 -DPARALLEL_FOR_LOOPS=1 -DSWAP_GRADE=0 -DPERTURB=0)
target_link_libraries(main_lazy_smart_sparse_chunk_parfor_clearing "-lboost_timer -lgmp -lgmpxx")

add_executable(main_lazy_smart_sparse_chunk_parfor_parmgkb "main.cpp")
target_compile_options(main_lazy_smart_sparse_chunk_parfor_parmgkb PUBLIC -DCHUNK_PREPROCESSING=1 -DSMART_REDUCTION=1 -DSPARSE_GRID_TRAVERSAL=1 -DCLEARING=0 -DMIN_GENS_AND_KER_BASIS_IN_PARALLEL=1 -DLAZY_MINIMIZATION=1 -DPARALLEL_FOR_LOOPS=1 -DSWAP_GRADE=0 -DPERTURB=0)
target_link_libraries(main_lazy_smart_sparse_chunk_parfor_parmgkb "-lboost_timer -lgmp -lgmpxx")


add_executable(main_vector_vector "main.cpp")
target_compile_options(main_vector_vector PUBLIC -DCHUNK_PREPROCESSING=1 -DSMART_REDUCTION=1 -DSPARSE_GRID_TRAVERSAL=1 -DCLEARING=0 -DMIN_GENS_AND_KER_BASIS_IN_PARALLEL=0 -DLAZY_MINIMIZATION=1 -DPARALLEL_FOR_LOOPS=0 -DSWAP_GRADE=0 -DPERTURB=0 -DCOLUMN_REPRESENTATION=vector_vector)
target_link_libraries(main_vector_vector "-lboost_timer -lgmp -lgmpxx")

add_executable(main_vector_vector_parfor "main.cpp")
target_compile_options(main_vector_vector_parfor PUBLIC -DCHUNK_PREPROCESSING=1 -DSMART_REDUCTION=1 -DSPARSE_GRID_TRAVERSAL=1 -DCLEARING=0 -DMIN_GENS_AND_KER_BASIS_IN_PARALLEL=0 -DLAZY_MINIMIZATION=1 -DPARALLEL_FOR_LOOPS=1 -DSWAP_GRADE=0 -DPERTURB=0 -DCOLUMN_REPRESENTATION=vector_vector)
target_link_libraries(main_vector_vector_parfor "-lboost_timer -lgmp -lgmpxx")

add_executable(main_vector_heap "main.cpp")
target_compile_options(main_vector_heap PUBLIC -DCHUNK_PREPROCESSING=1 -DSMART_REDUCTION=1 -DSPARSE_GRID_TRAVERSAL=1 -DCLEARING=0 -DMIN_GENS_AND_KER_BASIS_IN_PARALLEL=0 -DLAZY_MINIMIZATION=1 -DPARALLEL_FOR_LOOPS=0 -DSWAP_GRADE=0 -DPERTURB=0 -DCOLUMN_REPRESENTATION=vector_heap)
target_link_libraries(main_vector_heap "-lboost_timer -lgmp -lgmpxx")

add_executable(main_vector_heap_parfor "main.cpp")
target_compile_options(main_vector_heap_parfor PUBLIC -DCHUNK_PREPROCESSING=1 -DSMART_REDUCTION=1 -DSPARSE_GRID_TRAVERSAL=1 -DCLEARING=0 -DMIN_GENS_AND_KER_BASIS_IN_PARALLEL=0 -DLAZY_MINIMIZATION=1 -DPARALLEL_FOR_LOOPS=1 -DSWAP_GRADE=0 -DPERTURB=0 -DCOLUMN_REPRESENTATION=vector_heap)
target_link_libraries(main_vector_heap_parfor "-lboost_timer -lgmp -lgmpxx")

add_executable(main_vector_list "main.cpp")
target_compile_options(main_vector_list PUBLIC -DCHUNK_PREPROCESSING=1 -DSMART_REDUCTION=1 -DSPARSE_GRID_TRAVERSAL=1 -DCLEARING=0 -DMIN_GENS_AND_KER_BASIS_IN_PARALLEL=0 -DLAZY_MINIMIZATION=1 -DPARALLEL_FOR_LOOPS=0 -DSWAP_GRADE=0 -DPERTURB=0 -DCOLUMN_REPRESENTATION=vector_list)
target_link_libraries(main_vector_list "-lboost_timer -lgmp -lgmpxx")

add_executable(main_vector_list_parfor "main.cpp")
target_compile_options(main_vector_list_parfor PUBLIC -DCHUNK_PREPROCESSING=1 -DSMART_REDUCTION=1 -DSPARSE_GRID_TRAVERSAL=1 -DCLEARING=0 -DMIN_GENS_AND_KER_BASIS_IN_PARALLEL=0 -DLAZY_MINIMIZATION=1 -DPARALLEL_FOR_LOOPS=1 -DSWAP_GRADE=0 -DPERTURB=0 -DCOLUMN_REPRESENTATION=vector_list)
target_link_libraries(main_vector_list_parfor "-lboost_timer -lgmp -lgmpxx")

add_executable(main_vector_set "main.cpp")
target_compile_options(main_vector_set PUBLIC -DCHUNK_PREPROCESSING=1 -DSMART_REDUCTION=1 -DSPARSE_GRID_TRAVERSAL=1 -DCLEARING=0 -DMIN_GENS_AND_KER_BASIS_IN_PARALLEL=0 -DLAZY_MINIMIZATION=1 -DPARALLEL_FOR_LOOPS=0 -DSWAP_GRADE=0 -DPERTURB=0 -DCOLUMN_REPRESENTATION=vector_set)
target_link_libraries(main_vector_set "-lboost_timer -lgmp -lgmpxx")

add_executable(main_vector_set_parfor "main.cpp")
target_compile_options(main_vector_set_parfor PUBLIC -DCHUNK_PREPROCESSING=1 -DSMART_REDUCTION=1 -DSPARSE_GRID_TRAVERSAL=1 -DCLEARING=0 -DMIN_GENS_AND_KER_BASIS_IN_PARALLEL=0 -DLAZY_MINIMIZATION=1 -DPARALLEL_FOR_LOOPS=1 -DSWAP_GRADE=0 -DPERTURB=0 -DCOLUMN_REPRESENTATION=vector_set)
target_link_libraries(main_vector_set_parfor "-lboost_timer -lgmp -lgmpxx")

add_executable(main_full_pivot_column "main.cpp")
target_compile_options(main_full_pivot_column PUBLIC -DCHUNK_PREPROCESSING=1 -DSMART_REDUCTION=1 -DSPARSE_GRID_TRAVERSAL=1 -DCLEARING=0 -DMIN_GENS_AND_KER_BASIS_IN_PARALLEL=0 -DLAZY_MINIMIZATION=1 -DPARALLEL_FOR_LOOPS=0 -DSWAP_GRADE=0 -DPERTURB=0 -DCOLUMN_REPRESENTATION=full_pivot_column)
target_link_libraries(main_full_pivot_column "-lboost_timer -lgmp -lgmpxx")

add_executable(main_full_pivot_column_parfor "main.cpp")
target_compile_options(main_full_pivot_column_parfor PUBLIC -DCHUNK_PREPROCESSING=1 -DSMART_REDUCTION=1 -DSPARSE_GRID_TRAVERSAL=1 -DCLEARING=0 -DMIN_GENS_AND_KER_BASIS_IN_PARALLEL=0 -DLAZY_MINIMIZATION=1 -DPARALLEL_FOR_LOOPS=1 -DSWAP_GRADE=0 -DPERTURB=0 -DCOLUMN_REPRESENTATION=full_pivot_column)
target_link_libraries(main_full_pivot_column_parfor "-lboost_timer -lgmp -lgmpxx")

add_executable(main_heap_pivot_column "main.cpp")
target_compile_options(main_heap_pivot_column PUBLIC -DCHUNK_PREPROCESSING=1 -DSMART_REDUCTION=1 -DSPARSE_GRID_TRAVERSAL=1 -DCLEARING=0 -DMIN_GENS_AND_KER_BASIS_IN_PARALLEL=0 -DLAZY_MINIMIZATION=1 -DPARALLEL_FOR_LOOPS=0 -DSWAP_GRADE=0 -DPERTURB=0 -DCOLUMN_REPRESENTATION=heap_pivot_column)
target_link_libraries(main_heap_pivot_column "-lboost_timer -lgmp -lgmpxx")

add_executable(main_heap_pivot_column_parfor "main.cpp")
target_compile_options(main_heap_pivot_column_parfor PUBLIC -DCHUNK_PREPROCESSING=1 -DSMART_REDUCTION=1 -DSPARSE_GRID_TRAVERSAL=1 -DCLEARING=0 -DMIN_GENS_AND_KER_BASIS_IN_PARALLEL=0 -DLAZY_MINIMIZATION=1 -DPARALLEL_FOR_LOOPS=1 -DSWAP_GRADE=0 -DPERTURB=0 -DCOLUMN_REPRESENTATION=heap_pivot_column)
target_link_libraries(main_heap_pivot_column_parfor "-lboost_timer -lgmp -lgmpxx")

add_executable(main_sparse_pivot_column "main.cpp")
target_compile_options(main_sparse_pivot_column PUBLIC -DCHUNK_PREPROCESSING=1 -DSMART_REDUCTION=1 -DSPARSE_GRID_TRAVERSAL=1 -DCLEARING=0 -DMIN_GENS_AND_KER_BASIS_IN_PARALLEL=0 -DLAZY_MINIMIZATION=1 -DPARALLEL_FOR_LOOPS=0 -DSWAP_GRADE=0 -DPERTURB=0 -DCOLUMN_REPRESENTATION=sparse_pivot_column)
target_link_libraries(main_sparse_pivot_column "-lboost_timer -lgmp -lgmpxx")

add_executable(main_sparse_pivot_column_parfor "main.cpp")
target_compile_options(main_sparse_pivot_column_parfor PUBLIC -DCHUNK_PREPROCESSING=1 -DSMART_REDUCTION=1 -DSPARSE_GRID_TRAVERSAL=1 -DCLEARING=0 -DMIN_GENS_AND_KER_BASIS_IN_PARALLEL=0 -DLAZY_MINIMIZATION=1 -DPARALLEL_FOR_LOOPS=1 -DSWAP_GRADE=0 -DPERTURB=0 -DCOLUMN_REPRESENTATION=sparse_pivot_column)
target_link_libraries(main_sparse_pivot_column_parfor "-lboost_timer -lgmp -lgmpxx")

add_executable(main_bit_tree_pivot_column "main.cpp")
target_compile_options(main_bit_tree_pivot_column PUBLIC -DCHUNK_PREPROCESSING=1 -DSMART_REDUCTION=1 -DSPARSE_GRID_TRAVERSAL=1 -DCLEARING=0 -DMIN_GENS_AND_KER_BASIS_IN_PARALLEL=0 -DLAZY_MINIMIZATION=1 -DPARALLEL_FOR_LOOPS=0 -DSWAP_GRADE=0 -DPERTURB=0 -DCOLUMN_REPRESENTATION=bit_tree_pivot_column)
target_link_libraries(main_bit_tree_pivot_column "-lboost_timer -lgmp -lgmpxx")

add_executable(main_bit_tree_pivot_column_parfor "main.cpp")
target_compile_options(main_bit_tree_pivot_column_parfor PUBLIC -DCHUNK_PREPROCESSING=1 -DSMART_REDUCTION=1 -DSPARSE_GRID_TRAVERSAL=1 -DCLEARING=0 -DMIN_GENS_AND_KER_BASIS_IN_PARALLEL=0 -DLAZY_MINIMIZATION=1 -DPARALLEL_FOR_LOOPS=1 -DSWAP_GRADE=0 -DPERTURB=0 -DCOLUMN_REPRESENTATION=bit_tree_pivot_column)
target_link_libraries(main_bit_tree_pivot_column_parfor "-lboost_timer -lgmp -lgmpxx")




