# Benchmarks for multi-cover, perturbed instances (obsolete)
for x in 20 29 46
  do
    for i in 1 2 3 4 5
      do
        echo BENCHMARK FOR k_fold_${x}_10_${i}.txt
	echo RUNNING main_lazy_smart_perturbed
        /usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ./main_lazy_smart_perturbed multi_cover_data/k_fold_${x}_10_${i}.txt 
	echo RUNNING main_lazy_smart_sparse_perturbed
        /usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ./main_lazy_smart_sparse_perturbed multi_cover_data/k_fold_${x}_10_${i}.txt
      done
  done