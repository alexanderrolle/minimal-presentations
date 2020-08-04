echo "Starting:"
date

cmake .
make

./run_benchmarks_1.sh 2>&1 #| tee -a $1
#./run_benchmarks_2.sh 2>&1 #| tee -a $1
./run_benchmarks_3.sh 2>&1 #| tee -a $1
./run_benchmarks_4.sh 2>&1 #| tee -a $1
./run_benchmarks_5.sh 2>&1 #| tee -a $1
./run_benchmarks_6.sh 2>&1 #| tee -a $1
./run_benchmarks_7.sh 2>&1 #| tee -a $1
./run_benchmarks_8.sh 2>&1 #| tee -a $1
./run_benchmarks_9.sh 2>&1 #| tee -a $1

echo "Ending:"
date
