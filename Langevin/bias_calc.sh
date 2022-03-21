#!/bin/bash

sim=$1

for dt in 5e-3 6e-3 7e-3 8e-3 9e-3 1e-2 1.5e-2 2e-2
do
nohup /libre/blasseln/julia-1.7.2/bin/julia --threads=1 test_bias.jl 1.25 0.25 $dt 200.0 1000.0 3 100000 4.0 $sim > nohup.out &
done