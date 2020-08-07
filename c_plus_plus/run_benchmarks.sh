echo "Starting:"
date

cmake .
make

for i in 1 2 3 4 5 6 7
do
  echo EXECUTING BENCHMARK $i
  date
  ./run_benchmarks_$i.sh 2>&1 
done

echo "Ending:"
date
