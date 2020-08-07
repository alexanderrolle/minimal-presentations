# Benchmarks for mesh data from aim@shape
for x in hand eros
  do
        echo BENCHMARK FOR ${x}.off.firep
	echo RUNNING rivet
	./rivet_instance.sh from_off_fireps/${x}.off.firep compare.out
	echo RUNNING main_lazy_smart_sparse_chunk_parfor
        ./instance.sh ./main_lazy_smart_sparse_chunk_parfor from_off_fireps/${x}.off.firep reference.out
	diff -sq reference.out compare.out
	for prog in ./main_rivet ./main_lazy_smart_sparse ./main_lazy_smart_sparse_chunk
	do
	  echo RUNNING $prog
	  ./instance.sh $prog from_off_fireps/${x}.off.firep compare.out
	  diff -sq reference.out compare.out
	done
	echo RUNNING main_lazy_smart_sparse_chunk_parfor_clearing
        ./instance.sh ./main_lazy_smart_sparse_chunk_parfor_clearing from_off_fireps/${x}.off.firep compare.out
	#echo Diff of output files: 
        #diff -sq reference.out compare.out
	wc -l reference.out compare.out
  done

for x in dragon raptor
  do
        echo BENCHMARK FOR ${x}.off.firep
	echo RUNNING main_lazy_smart_sparse_chunk_parfor
        ./instance.sh ./main_lazy_smart_sparse_chunk_parfor from_off_fireps/${x}.off.firep reference.out
	for prog in ./main_lazy_smart_sparse ./main_lazy_smart_sparse_chunk
	do
	  echo RUNNING $prog
	  ./instance.sh $prog from_off_fireps/${x}.off.firep compare.out
	  diff -sq reference.out compare.out
	done
	echo RUNNING main_lazy_smart_sparse_chunk_parfor_clearing
        ./instance.sh ./main_lazy_smart_sparse_chunk_parfor_clearing from_off_fireps/${x}.off.firep compare.out
	#echo Diff of output files: 
        #diff -sq reference.out compare.out
	wc -l reference.out compare.out
  done