method=$1

for eta in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
do
nohup /libre/blasseln/julia-1.7.2/bin/julia --threads=8 estimate_mobility.jl 1.25 0.6 1e-3 1.0 $eta $method 100.0 5000000.0 10 BAOAB 2.5 > nohup.out &
done

