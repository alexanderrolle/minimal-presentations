# Benchmarks for functional rips
for x in 9 10 11 12 13
  do
    for i in 0 1 2 3 4
      do
        echo BENCHMARK FOR noisy_circle_firep_${x}_${i}.txt
	echo RUNNING main_rivet
        /usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ./main_rivet noisy_circle_fireps/noisy_circle_firep_${x}_${i}.txt reference.out
	echo RUNNING main_lazy_smart_sparse
        /usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ./main_lazy_smart_sparse noisy_circle_fireps/noisy_circle_firep_${x}_${i}.txt compare.out
	echo Diff of output files: 
        diff -sq reference.out compare.out
	echo RUNNING main_lazy_smart_sparse_chunk
        /usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ./main_lazy_smart_sparse_chunk noisy_circle_fireps/noisy_circle_firep_${x}_${i}.txt compare.out
	echo Diff of output files: 
        diff -sq reference.out compare.out
	echo RUNNING main_lazy_smart_sparse_chunk_parfor
        /usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ./main_lazy_smart_sparse_chunk_parfor noisy_circle_fireps/noisy_circle_firep_${x}_${i}.txt compare.out
	echo Diff of output files: 
        diff -sq reference.out compare.out
	echo RUNNING main_lazy_smart_sparse_chunk_parfor_clearing
        /usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ./main_lazy_smart_sparse_chunk_parfor_clearing noisy_circle_fireps/noisy_circle_firep_${x}_${i}.txt compare.out
	#echo Diff of output files: 
        wc -l reference.out compare.out
      done
  done

for x in 9 10 11 12 13
  do
    for i in 0 1 2 3 4
      do
        echo BENCHMARK FOR noisy_sphere_2_firep_${x}_${i}.txt
	echo RUNNING main_rivet
        /usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ./main_rivet noisy_2_sphere_fireps/noisy_2_sphere_firep_${x}_${i}.txt reference.out
	echo RUNNING main_lazy_smart_sparse
        /usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ./main_lazy_smart_sparse noisy_2_sphere_fireps/noisy_2_sphere_firep_${x}_${i}.txt compare.out
	echo Diff of output files: 
        diff -sq reference.out compare.out
	echo RUNNING main_lazy_smart_sparse_chunk
        /usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ./main_lazy_smart_sparse_chunk noisy_2_sphere_fireps/noisy_2_sphere_firep_${x}_${i}.txt compare.out
	echo Diff of output files: 
        diff -sq reference.out compare.out
	echo RUNNING main_lazy_smart_sparse_chunk_parfor
        /usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ./main_lazy_smart_sparse_chunk_parfor noisy_2_sphere_fireps/noisy_2_sphere_firep_${x}_${i}.txt compare.out
	echo Diff of output files: 
        diff -sq reference.out compare.out
	echo RUNNING main_lazy_smart_sparse_chunk_parfor_clearing
        /usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ./main_lazy_smart_sparse_chunk_parfor_clearing noisy_2_sphere_fireps/noisy_2_sphere_firep_${x}_${i}.txt compare.out
	#echo Diff of output files: 
        wc -l reference.out compare.out
      done
  done
