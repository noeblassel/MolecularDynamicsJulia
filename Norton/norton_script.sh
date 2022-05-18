method=$1

v_min=$2
v_max=$3
v_inc=$4

for v in `LANG=en_US seq $v_min $v_inc $v_max`
do
nohup /libre/blasseln/julia-1.7.2/bin/julia --threads=8 estimate_mobility.jl 1.25 0.6 1e-3 1.0 $v $method 100.0 5000000.0 10 BAOAB 2.5 > nohup.out &
done

