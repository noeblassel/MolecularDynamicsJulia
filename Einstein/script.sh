#!/bin/bash
N=$1

for i in {1..$N}
do
nohup /libre/blasseln/julia-1.7.2/bin/julia --threads=4 compute_msd.jl 2.5 0.7 5e-3 1.0 100.0 1000000.0 1000.0 10 BAOAB 3.0 > nohup$i.out &
done

