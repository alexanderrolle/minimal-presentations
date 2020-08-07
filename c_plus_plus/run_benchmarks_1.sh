#Benchmarks for points on sphere, small data
for x in 7500 15000 30000 60000
  do
    for i in 1 2 3 4 5
      do
        echo BENCHMARK FOR points_on_sphere_firep_${x}_${i}.txt
	echo RUNNING rivet
	./rivet_instance.sh points_on_sphere_fireps/instance_${x}_${i}.firep compare.out
	echo RUNNING main_rivet
        ./instance.sh ./main_rivet points_on_sphere_fireps/instance_${x}_${i}.firep reference.out
        diff -sq reference.out compare.out
	for prog in ./main_smart_sparse ./main_lazy_smart_sparse ./main_lazy_smart_sparse_chunk ./main_lazy_smart_sparse_chunk_parfor
	do
	  echo RUNNING $prog
	  ./instance.sh $prog points_on_sphere_fireps/instance_${x}_${i}.firep compare.out
	  diff -sq reference.out compare.out	  
	done
      done
  done
