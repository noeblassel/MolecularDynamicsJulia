#!/bin/bash

sim=$1

for dt in 1e-2 6e-3 7e-3 8e-3 9e-3
do
nohup /libre/blasseln/julia-1.7.2/bin/julia test_bias.jl 1.5 0.15 $dt 50.0 1000.0 4 16000 4.0 $sim bias_test_64.out > nohup$dt.out &
done