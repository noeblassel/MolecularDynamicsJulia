#!/bin/bash

simulator=BABO

for dt in 1e-6 2e-6 3e-6 4e-6 5e-6 7e-6 8e-6 9e-6 1e-5 2e-5 3e-5 4e-5 5e-5
do
nohup /libre/blasseln/julia-1.7.2/bin/julia test_bias.jl 1.25 0.25 $dt 500 $sim bias_test.out > nohup.out &
done