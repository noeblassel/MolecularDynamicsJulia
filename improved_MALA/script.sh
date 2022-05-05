rule=$1
proposal=$2

for lg_dt in -5.0 -4.9 -4.8 -4.7 -4.6 -4.5 -4.4 -4.3 -4.2 -4.1 -4.0
do nohup /libre/blasseln/julia-1.7.2/bin/julia --threads=8 compute_diffusion.jl 4 $lg_dt 0.4 1000000 0.4 $rule $proposal>nohup.out&
done
