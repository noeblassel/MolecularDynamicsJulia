method=$1

for eta in 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1
do
nohup /libre/blasseln/julia-1.7.2/bin/julia --threads=8 estimate_mobility.jl 1.25 0.6 1e-3 1.0 $eta $method 100.0 5000000.0 10 BAOAB 2.5 > nohup.out &
done

