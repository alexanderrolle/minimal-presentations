# Benchmarks for multi-cover, large instances
for x in 640 1250 2455 4990
  do
    for i in 1 2 3 4 5
      do
        echo BENCHMARK FOR k_fold_${x}_10_${i}.txt
	echo RUNNING main_lazy_smart_sparse
        /usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ./main_lazy_smart_sparse multi_cover_data/k_fold_${x}_10_${i}.txt reference.out
	echo RUNNING main_lazy_smart_sparse_chunk
        /usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ./main_lazy_smart_sparse_chunk multi_cover_data/k_fold_${x}_10_${i}.txt compare.out
	diff -sq reference.out compare.out
	echo RUNNING main_lazy_smart_sparse_chunk_parfor
        /usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ./main_lazy_smart_sparse_chunk_parfor multi_cover_data/k_fold_${x}_10_${i}.txt compare.out
	diff -sq reference.out compare.out
	echo RUNNING main_lazy_smart_sparse_chunk_parfor_clearing
        /usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ./main_lazy_smart_sparse_chunk_parfor_clearing multi_cover_data/k_fold_${x}_10_${i}.txt compare.out
	echo Compare file sizes:
	wc -l reference.out compare.out
      done
  done