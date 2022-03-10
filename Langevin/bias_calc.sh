#!/bin/bash

sim=$1

for dt in 1e-4 2e-4 3e-4 4e-4 5e-4
do
nohup /libre/blasseln/julia-1.7.2/bin/julia test_bias.jl 1.25 0.25 $dt 10000 $sim bias_test.out > nohup$dt.out &
done