# Benchmarks for functional rips
for x in 9 10 11 12 13
  do
    for i in 0 1 2 3 4
      do
        echo BENCHMARK FOR noisy_circle_firep_${x}_${i}.txt
	echo RUNNING rivet
	./rivet_instance.sh noisy_circle_fireps/noisy_circle_firep_${x}_${i}.txt compare.out
	echo RUNNING main_rivet
        ./instance.sh ./main_rivet noisy_circle_fireps/noisy_circle_firep_${x}_${i}.txt reference.out
	diff -sq reference.out compare.out 
	for prog in ./main_lazy_smart_sparse ./main_lazy_smart_sparse_chunk ./main_lazy_smart_sparse_chunk_parfor
	do
	  echo RUNNING $prog
	  ./instance.sh $prog noisy_circle_fireps/noisy_circle_firep_${x}_${i}.txt compare.out
	  diff -sq reference.out compare.out
	done
	echo RUNNING main_lazy_smart_sparse_chunk_parfor_clearing
        ./instance.sh ./main_lazy_smart_sparse_chunk_parfor_clearing noisy_circle_fireps/noisy_circle_firep_${x}_${i}.txt compare.out
	#echo Diff of output files: 
        wc -l reference.out compare.out
      done
  done

for x in 8 9 10 11 12
  do
    for i in 0 1 2 3 4
      do
        echo BENCHMARK FOR noisy_sphere_2_firep_${x}_${i}.txt
	echo RUNNING rivet
	./rivet_instance.sh noisy_2_sphere_fireps/noisy_2_sphere_firep_${x}_${i}.txt compare.out
	echo RUNNING main_rivet
        ./instance.sh ./main_rivet noisy_2_sphere_fireps/noisy_2_sphere_firep_${x}_${i}.txt reference.out
	diff -sq reference.out compare.out
	for prog in ./main_lazy_smart_sparse ./main_lazy_smart_sparse_chunk ./main_lazy_smart_sparse_chunk_parfor
	do
	  echo RUNNING $prog
	  ./instance.sh $prog noisy_2_sphere_fireps/noisy_2_sphere_firep_${x}_${i}.txt compare.out
	  diff -sq reference.out compare.out
	done
	echo RUNNING main_lazy_smart_sparse_chunk_parfor_clearing
        ./instance.sh ./main_lazy_smart_sparse_chunk_parfor_clearing noisy_2_sphere_fireps/noisy_2_sphere_firep_${x}_${i}.txt compare.out
	#echo Diff of output files: 
        wc -l reference.out compare.out
      done
  done
