/usr/bin/time --format "RESULTS %C\nTime: %e\nMemory: %M\nMemorySwaps: %W" timeout 3600 ~/rivet/rivet_console $1 --minpres --num_threads=1  > $2