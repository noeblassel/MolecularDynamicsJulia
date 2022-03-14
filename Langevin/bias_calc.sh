#!/bin/bash

sim=$1

for dt in 6e-3 7e-3 8e-3 9e-3 1e-2
do
nohup /libre/blasseln/julia-1.7.2/bin/julia test_bias.jl 1.25 0.25 $dt 50.0 1000.0 5 2000 4.0 $sim bias_test.out > nohup$dt.out &
done