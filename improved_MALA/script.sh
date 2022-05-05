rule=$1
proposal=$2

for lg_dt in -5.0 -4.8 -4.6 -4.4 -4.2
do nohup /libre/blasseln/julia-1.7.2/bin/julia --threads=8 compute_diffusion.jl 4 $lg_dt 0.4 1000000 0.75 $rule $proposal>nohup.out&
done
