# Benchmarks for mesh data from aim@shape
for x in hand eros
  do
        echo BENCHMARK FOR ${x}.off.firep
	echo RUNNING main_lazy_smart_sparse_chunk_parfor
        /usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ./main_lazy_smart_sparse_chunk_parfor from_off_fireps/${x}.off.firep reference.out
	echo RUNNING main_rivet
        /usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ./main_rivet from_off_fireps/${x}.off.firep compare.out
	echo Diff of output files: 
        diff -sq reference.out compare.out
	echo RUNNING main_lazy_smart_sparse
        /usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ./main_lazy_smart_sparse from_off_fireps/${x}.off.firep compare.out
	echo Diff of output files: 
        diff -sq reference.out compare.out
	echo RUNNING main_lazy_smart_sparse_chunk
        /usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ./main_lazy_smart_sparse_chunk from_off_fireps/${x}.off.firep compare.out
	echo Diff of output files: 
        diff -sq reference.out compare.out
	echo RUNNING main_lazy_smart_sparse_chunk_parfor_clearing
        /usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ./main_lazy_smart_sparse_chunk_parfor_clearing from_off_fireps/${x}.off.firep compare.out
	#echo Diff of output files: 
        #diff -sq reference.out compare.out
	wc -l reference.out compare.out
  done

for x in dragon raptor
  do
        echo BENCHMARK FOR ${x}.off.firep
	echo RUNNING main_lazy_smart_sparse_chunk_parfor
        /usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ./main_lazy_smart_sparse_chunk_parfor from_off_fireps/${x}.off.firep reference.out
	echo RUNNING main_lazy_smart_sparse
        /usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ./main_lazy_smart_sparse from_off_fireps/${x}.off.firep compare.out
	echo Diff of output files: 
        diff -sq reference.out compare.out
	echo RUNNING main_lazy_smart_sparse_chunk
        /usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ./main_lazy_smart_sparse_chunk from_off_fireps/${x}.off.firep compare.out
	echo Diff of output files: 
        diff -sq reference.out compare.out
	echo RUNNING main_lazy_smart_sparse_chunk_parfor_clearing
        /usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ./main_lazy_smart_sparse_chunk_parfor_clearing from_off_fireps/${x}.off.firep compare.out
	#echo Diff of output files: 
        #diff -sq reference.out compare.out
	wc -l reference.out compare.out
  done