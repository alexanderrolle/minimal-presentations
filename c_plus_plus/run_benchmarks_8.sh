#Benchmarks for points on sphere, large data
for x in 200000 400000 800000 1600000
  do
    for i in 1 2 3 4 5
      do
        echo BENCHMARK FOR points_on_sphere_firep_${x}_${i}.txt
	echo RUNNING main_lazy_smart_sparse
        /usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ./main_lazy_smart_sparse points_on_sphere_fireps/instance_${x}_${i}.firep reference.out
	echo RUNNING main_lazy_smart_sparse_chunk
        /usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ./main_lazy_smart_sparse_chunk points_on_sphere_fireps/instance_${x}_${i}.firep compare.out
	echo Diff of output files: 
        diff -sq reference.out compare.out
	echo RUNNING main_lazy_smart_sparse_chunk_parfor
        /usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ./main_lazy_smart_sparse_chunk_parfor points_on_sphere_fireps/instance_${x}_${i}.firep compare.out
	echo Diff of output files: 
        diff -sq reference.out compare.out
	echo RUNNING main_lazy_smart_sparse_chunk_parfor_parmgkb
        /usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ./main_lazy_smart_sparse_chunk_parfor_parmgkb points_on_sphere_fireps/instance_${x}_${i}.firep compare.out
	echo Diff of output files: 
        diff -sq reference.out compare.out
      done
  done