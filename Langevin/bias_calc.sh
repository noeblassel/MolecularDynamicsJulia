#!/bin/bash

sim=$1

for dt in 1e-1 2e-1 3e-1 4e-1 5e-1
do
nohup /libre/blasseln/julia-1.7.2/bin/julia test_bias.jl 1.5 0.3 $dt 20.0 2000.0 3 40000 4.0 $sim bias_test_64.$sim > nohup$dt.out &
done