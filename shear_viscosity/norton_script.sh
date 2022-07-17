method=$1

eta_min=$2
eta_max=$3
eta_inc=$4

for eta in `LANG=en_US seq $eta_min $eta_inc $eta_max`
do
nohup /libre/blasseln/julia-1.7.2/bin/julia --threads=8 shear_viscosity.jl 0.8 0.7 1e-3 1.0 $eta $method 100.0 5000000 10 2.5 > nohup.out &
done

