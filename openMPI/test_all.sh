#!/bin/bash

# Define the list of input pairs as an array
input_pairs=(ex1:20 ex2:40 gen10:20 gen15:25 gen19:23 gen20:5000 gen22:25000 gen24:50000 gen26:50000 gen30:5000)

# Loop over the array and run the program for each pair
for (( i=0; i<${#input_pairs[@]}; i++ )); do
    pair=${input_pairs[$i]}
    Y=${pair%:*}
    X=${pair#*:}
    echo "Running program for pair ($Y,$X)..."
    ./tsp pub-instances/$Y-$X.in $X
done
