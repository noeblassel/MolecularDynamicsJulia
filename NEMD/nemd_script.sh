method=$1

for eta in 0.01 0.02 0.03 0.04 0.05
do
nohup /libre/blasseln/julia-1.7.2/bin/julia --threads=8 estimate_mobility.jl 2.5 0.7 5e-3 1.0 $eta $method 100.0 5000000.0 10 BAOAB 3.0 > nohup.out &
done

