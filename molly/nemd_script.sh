method=$1

for eta in 0.1 0.2 0.3 0.4 0.5
do
nohup /libre/blasseln/julia-1.7.2/bin/julia --threads=3 estimate_mobility.jl 2.5 0.7 5e-3 1.0 $eta $1 100.0 5000.0 10 BAOAB 3.0 > nohup.out &
done