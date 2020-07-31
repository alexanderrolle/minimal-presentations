# Benchmarks for multi-cover, small instances
for x in 96 175 330 640
  do
    for i in 1 2 3 4 5
      do
        echo BENCHMARK FOR k_fold_${x}_10_${i}.txt
	echo RUNNING RIVET
	/usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ~/rivet/rivet_console ~/minimal-presentations/c_plus_plus/multi_cover_data/k_fold_${x}_10_${i}.txt --minpres --num_threads=1 > reference.out
	echo RUNNING main_rivet
        /usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ./main_rivet multi_cover_data/k_fold_${x}_10_${i}.txt compare.out
	diff -sq reference.out compare.out
	echo RUNNING main_lazy
        /usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ./main_lazy multi_cover_data/k_fold_${x}_10_${i}.txt compare.out
	diff -sq reference.out compare.out
	echo RUNNING main_lazy_smart
        /usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ./main_lazy_smart multi_cover_data/k_fold_${x}_10_${i}.txt compare.out
	diff -sq reference.out compare.out
	echo RUNNING main_lazy_smart_sparse
        /usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ./main_lazy_smart_sparse multi_cover_data/k_fold_${x}_10_${i}.txt compare.out
	diff -sq reference.out compare.out
      done
  done