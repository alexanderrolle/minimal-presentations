cmake .
make

for x in 1250 2455 4990
  do
    for i in 1 2 3 4 5
      do
        echo BENCHMARK FOR k_fold_${x}_10_${i}.txt
	echo RUNNING main_lazy_smart_sparse
        /usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ./main_lazy_smart_sparse multi_cover_data/k_fold_${x}_10_${i}.txt k_fold_${x}_10_${i}_result_main.txt
	echo RUNNING main_lazy_smart_sparse_chunk
        /usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ./main_lazy_smart_sparse_chunk multi_cover_data/k_fold_${x}_10_${i}.txt k_fold_${x}_10_${i}_result_tmp.txt
	echo Diff of output files: 
        diff k_fold_${x}_10_${i}_result_main.txt k_fold_${x}_10_${i}_result_tmp.txt
	echo RUNNING main_lazy_smart_sparse_chunk_parfor
        /usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ./main_lazy_smart_sparse_chunk_parfor multi_cover_data/k_fold_${x}_10_${i}.txt k_fold_${x}_10_${i}_result_tmp.txt
	echo Diff of output files: 
        diff k_fold_${x}_10_${i}_result_main.txt k_fold_${x}_10_${i}_result_tmp.txt
	echo RUNNING main_lazy_smart_sparse_chunk_parfor_clearing
        /usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ./main_lazy_smart_sparse_chunk_parfor_clearing multi_cover_data/k_fold_${x}_10_${i}.txt k_fold_${x}_10_${i}_result_tmp.txt
	echo Compare file sizes:
	wc -l  k_fold_${x}_10_${i}_result_main.txt k_fold_${x}_10_${i}_result_tmp.txt
        #diff k_fold_${x}_10_${i}_result_main.txt k_fold_${x}_10_${i}_result_tmp.txt
      done
  done