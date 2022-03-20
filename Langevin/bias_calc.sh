#!/bin/bash

sim=$1

for dt in 1e-1 2e-1 3e-1 4e-1 5e-1 1e-2 2e-2 3e-2 4e-2 5e-2 6e-2 7e-2 8e-2 9e-2
do
nohup /libre/blasseln/julia-1.7.2/bin/julia --threads=1 test_bias.jl 1.5 0.3 $dt 20.0 2000.0 3 100 4.0 $sim > nohup.out &
done