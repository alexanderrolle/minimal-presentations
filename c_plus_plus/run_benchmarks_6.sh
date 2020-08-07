# Benchmarks for multi-cover, small instances
for x in 96 175 330 640
  do
    for i in 1 2 3 4 5
      do
        echo BENCHMARK FOR k_fold_${x}_10_${i}.txt
	echo RUNNING RIVET
	./rivet_instance.sh  ~/minimal-presentations/c_plus_plus/multi_cover_data/k_fold_${x}_10_${i}.txt reference.out
	for prog in ./main_rivet ./main_smart_sparse ./main_lazy_smart_sparse ./main_lazy_smart_sparse_chunk ./main_lazy_smart_sparse_chunk_parfor 
	do
	  echo RUNNING $prog
	  ./instance.sh $prog multi_cover_data/k_fold_${x}_10_${i}.txt compare.out
	  diff -sq reference.out compare.out
	done
	echo RUNNING main_lazy_smart_sparse_chunk_parfor_clearing
        /usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ./main_lazy_smart_sparse_chunk_parfor_clearing multi_cover_data/k_fold_${x}_10_${i}.txt compare.out
	echo Compare file sizes:
	wc -l reference.out compare.out
      done
  done

# Benchmarks for multi-cover, large instances
for x in 1250 2455 4990
  do
    for i in 1 2 3 4 5
      do
        echo BENCHMARK FOR k_fold_${x}_10_${i}.txt
	echo RUNNING main_lazy_smart_sparse
        ./instance.sh ./main_lazy_smart_sparse multi_cover_data/k_fold_${x}_10_${i}.txt reference.out
	for prog in ./main_lazy_smart_sparse_chunk ./main_lazy_smart_sparse_chunk_parfor ./main_lazy_smart_sparse_chunk_parfor_parmgkb
	do
	  echo RUNNING $prog
	  ./instance.sh $prog multi_cover_data/k_fold_${x}_10_${i}.txt compare.out
	  diff -sq reference.out compare.out
	done
	echo RUNNING main_lazy_smart_sparse_chunk_parfor_clearing
        ./instance.sh ./main_lazy_smart_sparse_chunk_parfor_clearing multi_cover_data/k_fold_${x}_10_${i}.txt compare.out
	#echo Compare file sizes:
	wc -l reference.out compare.out
      done
  done
