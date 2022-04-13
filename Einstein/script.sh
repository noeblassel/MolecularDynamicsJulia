#!/bin/bash
N=$1

for i in $( eval echo {0..$N} )
do
nohup /libre/blasseln/julia-1.7.2/bin/julia --threads=4 compute_msd.jl 2.5 0.7 5e-3 1.0 100.0 1000000.0 1000.0 10 BAOAB 3.0 msd_history$i.out > nohup$i.out &
done

