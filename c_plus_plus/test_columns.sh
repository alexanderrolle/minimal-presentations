echo RUNNING main_vector_vector
./instance.sh ./main_vector_vector $1 reference.out

echo RUNNING main_vector_vector_parfor 
/usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ./main_vector_vector_parfor $1 compare.out 
diff -sq reference.out compare.out

for col in vector_heap vector_list vector_set full_pivot_column sparse_pivot_column heap_pivot_column bit_tree_pivot_column
do
  echo RUNNING main_$col
  ./instance.sh ./main_$col $1 compare.out
  diff -sq reference.out compare.out
  echo RUNNING main_${col}_parfor 
  ./instance.sh ./main_${col}_parfor $1 compare.out 
  diff -sq reference.out compare.out
done 






