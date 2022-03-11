#!/bin/bash

sim=$1

for dt in 1e-3 2e-3 3e-3 5e-4 7e-4
do
nohup /libre/blasseln/julia-1.7.2/bin/julia test_bias.jl 1.25 0.25 $dt 50.0 500.0 5 100 4.0 $sim bias_test.out > nohup$dt.out &
done