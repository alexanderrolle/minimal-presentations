#Benchmarks for random delaunay, small data
for dim in 1 2
  do
  for x in 5000 10000 20000 40000
    do
      for i in 1 2 3 4 5
	do
          echo BENCHMARK FOR random_del_firep_${x}_${i}_${dim}.txt
	  echo RUNNING main_lazy_smart_sparse_chunk_parfor
	  ./instance.sh ./main_lazy_smart_sparse_chunk_parfor random_del_fireps/instance_${x}_${i}_dim_${dim}.firep reference.out
	  echo RUNNING rivet
	  ./rivet_instance.sh random_del_fireps/instance_${x}_${i}_dim_${dim}.firep compare.out
	  diff -sq reference.out compare.out
	  for prog in ./main_rivet ./main_smart_sparse ./main_lazy_smart_sparse ./main_lazy_smart_sparse_chunk ./main_lazy_smart_sparse_chunk_parfor
	  do
	    echo RUNNING $prog
	    ./instance.sh $prog random_del_fireps/instance_${x}_${i}_dim_${dim}.firep reference.out
	    diff -sq reference.out compare.out
	  done
	done
      done
  done

#Benchmark for random delaunay, large data

for dim in 1 2
  do
  for x in 80000 160000 320000 640000
    do
      for i in 1 2 3 4 5
	do
          echo BENCHMARK FOR random_del_firep_${x}_${i}_dim_${dim}.txt
	  echo RUNNING main_lazy_smart_sparse
	  ./instance.sh ./main_lazy_smart_sparse random_del_fireps/instance_${x}_${i}_dim_${dim}.firep reference.out
	  for prog in ./main_lazy_smart_sparse_chunk ./main_lazy_smart_sparse_chunk_parfor ./main_lazy_smart_sparse_chunk_parfor_parmgkb
	  do
	    echo RUNNING $prog
	    ./instance.sh $prog random_del_fireps/instance_${x}_${i}_dim_${dim}.firep compare.out
	    diff -sq reference.out compare.out
	  done
	  echo RUNNING main_lazy_smart_sparse_chunk_parfor_clearing
	  ./instance.sh ./main_lazy_smart_sparse_chunk_parfor_clearing random_del_fireps/instance_${x}_${i}_dim_${dim}.firep compare.out
	  wc -l reference.out compare.out
	  #diff -sq reference.out compare.out
	done
      done
  done
