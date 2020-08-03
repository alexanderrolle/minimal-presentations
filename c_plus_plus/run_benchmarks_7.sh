#Benchmarks for points on sphere, small data
for x in 7500 15000 30000 60000
  do
    for i in 1 2 3 4 5
      do
        echo BENCHMARK FOR points_on_sphere_firep_${x}_${i}.txt
	echo RUNNING rivet
	timeout 1800 /usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ~/rivet/rivet_console points_on_sphere_fireps/instance_${x}_${i}.firep --minpres --num_threads=1 > compare.out
	echo RUNNING main_rivet
        timeout 1800 /usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ./main_rivet points_on_sphere_fireps/instance_${x}_${i}.firep reference.out
	echo Diff of output files: 
        diff -sq reference.out compare.out
	echo RUNNING main_smart_sparse
        timeout 1800 /usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ./main_smart_sparse points_on_sphere_fireps/instance_${x}_${i}.firep compare.out
	echo Diff of output files: 
        diff -sq reference.out compare.out
	echo RUNNING main_lazy_smart_sparse
        timeout 1800 /usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ./main_lazy_smart_sparse points_on_sphere_fireps/instance_${x}_${i}.firep compare.out
	echo Diff of output files: 
        diff -sq reference.out compare.out
	echo RUNNING main_lazy_smart_sparse_chunk
        timeout 1800 /usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ./main_lazy_smart_sparse_chunk points_on_sphere_fireps/instance_${x}_${i}.firep compare.out
	echo Diff of output files: 
        diff -sq reference.out compare.out
	echo RUNNING main_lazy_smart_sparse_chunk_parfor
        timeout 1800 /usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ./main_lazy_smart_sparse_chunk_parfor points_on_sphere_fireps/instance_${x}_${i}.firep compare.out
	echo Diff of output files: 
        diff -sq reference.out compare.out
      done
  done