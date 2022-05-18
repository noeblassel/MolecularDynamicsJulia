#!/bin/bash
ssh clustern14 -X  'cd /libre/blasseln/MolecularDynamicsJulia/Norton/results; /libre/blasseln/julia-1.7.2/bin/julia plot_mobility_estimates.jl; exit'
scp clustern14:/libre/blasseln/MolecularDynamicsJulia/Norton/results/*.pdf .
