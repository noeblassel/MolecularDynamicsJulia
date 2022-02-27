declare -A rhomins
declare -A rhomaxs

rhomins[17]=0.45
rhomins[18]=0.5
rhomins[23]=0.55
rhomins[24]=0.6
rhomis[25]=0.65

rhomaxs[17]=0.5
rhomaxs[18]=0.55
rhomaxs[23]=0.6
rhomaxs[24]=0.65
rhomaxs[25]=0.7

for n in 17 18 23 24 25
do
ssh clustern$n "cd /libre/blasseln/MolecularDynamicsJulia/Langevin; git pull; nohup /libre/blasseln/julia-1.7.2/bin/julia test_nist.jl ${rhomins[$n]} ${rhomaxs[$n]} 12 >nohup.out&"
done
