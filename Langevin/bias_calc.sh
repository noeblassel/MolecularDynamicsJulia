#!/bin/bash

sim=$1

for dt in 2e-3 3e-3 4e-3 5e-3 6e-3 7e-3
do
nohup /libre/blasseln/julia-1.7.2/bin/julia --threads=1 test_bias.jl 1.25 0.25 $dt 200.0 1000.0 3 100000 $sim > nohup.out &
done