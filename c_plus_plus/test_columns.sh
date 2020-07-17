echo RUNNING main_vector_vector
/usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ./main_vector_vector $1 reference.out
echo RUNNING main_vector_heap
/usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ./main_vector_heap $1 result.out
if  cmp -s reference.out result.out
then
    echo "OUTPUTS DIFFER"
    exit 1
fi
echo RUNNING main_vector_list
/usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ./main_vector_list $1 result.out
diff -s reference.out result.out
echo RUNNING main_vector_set
/usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ./main_vector_set $1 result.out
diff -s reference.out result.out
echo RUNNING main_full_pivot_column
/usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ./main_full_pivot_column $1 result.out
diff -s reference.out result.out
echo RUNNING main_sparse_pivot_column
/usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ./main_sparse_pivot_column $1 result.out
diff -s reference.out result.out
echo RUNNING main_heap_pivot_column
/usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ./main_heap_pivot_column $1 result.out
diff -s reference.out result.out
echo RUNNING main_bit_tree_pivot_column
/usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" ./main_bit_tree_pivot_column $1 result.out
diff -s reference.out result.out





