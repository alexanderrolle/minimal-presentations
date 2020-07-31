# Benchmarks for column types
for x in multi_cover_data/k_fold_1250_10_1.txt multi_cover_data/k_fold_4990_10_1.txt noisy_circle_fireps/noisy_circle_firep_10_1.txt noisy_circle_fireps/noisy_circle_firep_13_1.txt from_off_fireps/hand.off.firep from_off_fireps/dragon.off.firep
  do
        echo BENCHMARK FOR ${x}
	./test_columns.sh $x
  done