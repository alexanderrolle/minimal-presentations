#Benchmarks for points on sphere, large data
for x in 200000 400000 800000 1600000
  do
    for i in 1 2 3 4 5
      do
        echo BENCHMARK FOR points_on_sphere_firep_${x}_${i}.txt
	echo RUNNING main_lazy_smart_sparse
        ./instance.sh ./main_lazy_smart_sparse points_on_sphere_fireps/instance_${x}_${i}.firep reference.out
	for prog in ./main_lazy_smart_sparse_chunk ./main_lazy_smart_sparse_chunk_parfor ./main_lazy_smart_sparse_chunk_parfor_parmgkb
	do
	  echo RUNNING $prog
	  ./instance.sh $prog points_on_sphere_fireps/instance_${x}_${i}.firep compare.out
	  diff -sq reference.out compare.out
	done
	echo RUNNING main_lazy_smart_sparse_chunk_parfor_clearing
        ./instance.sh ./main_lazy_smart_sparse_chunk_parfor_clearing points_on_sphere_fireps/instance_${x}_${i}.firep compare.out
	wc -l reference.out compare.out
      done
  done
