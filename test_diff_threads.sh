#!/bin/bash
# Number of threads to test on
n_threads=(1 2 4 6 8 12)

for t in "${n_threads[@]}"; do
    echo "Using ${t} threads"
    # Loop over the array and run the program for each pair
    bash test_all.sh
    echo
done
