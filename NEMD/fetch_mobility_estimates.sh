#!/bin/bash
#15 : color drift
#16: one drift

declare -A nodeA
declare -A nodeB

nodeA[color]=15
nodeA[single]=16

nodeB[color]=17
nodeB[single]=19

mode=$1 #color or single

scp clustern${nodeA[$mode]}:/libre/blasseln/MolecularDynamicsJulia/NEMD/mobility*.out ./results/
#scp clustern${nodeB[$mode]}:/libre/blasseln/MolecularDynamicsJulia/NEMD/mobility*.out ./mobility_computations_$mode/