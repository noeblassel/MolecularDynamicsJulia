#!/bin/bash

sim=$1

for dt in 1e-4 2e-4 3e-4
do
nohup /libre/blasseln/julia-1.7.2/bin/julia test_bias.jl 1.25 0.25 $dt 500 $sim bias_test.out > nohup$dt.out &
rm nohup$dt.out
done

for dt in 4e-4 5e-4 6e-4
do
nohup /libre/blasseln/julia-1.7.2/bin/julia test_bias.jl 1.25 0.25 $dt 500 $sim bias_test.out > nohup$dt.out &
rm nohup$dt.out
done

for dt in 7e-4 8e-4 9e-4
do
nohup /libre/blasseln/julia-1.7.2/bin/julia test_bias.jl 1.25 0.25 $dt 500 $sim bias_test.out > nohup$dt.out &
rm nohup$dt.out
done

for dt in 1e-3 2e-3 3e-3
do
nohup /libre/blasseln/julia-1.7.2/bin/julia test_bias.jl 1.25 0.25 $dt 500 $sim bias_test.out > nohup$dt.out &
rm nohup$dt.out
done

for dt in 4e-3 5e-3
do
nohup /libre/blasseln/julia-1.7.2/bin/julia test_bias.jl 1.25 0.25 $dt 500 $sim bias_test.out > nohup$dt.out &
rm nohup$dt.out
done