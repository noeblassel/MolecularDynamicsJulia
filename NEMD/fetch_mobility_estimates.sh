#!/bin/bash
#16, 19 : single drift
#17, 14 : color drift

declare -A nodeA
declare -A nodeB

nodeA[color]=14
nodeA[single]=16

nodeB[color]=17
nodeB[single]=19

mode=$1 #color or single

scp clustern${nodeA[$mode]}:/libre/blasseln/MolecularDynamicsJulia/NEMD/mobility*.out ./mobility_computations_$mode/
scp clustern${nodeB[$mode]}:/libre/blasseln/MolecularDynamicsJulia/NEMD/mobility*.out ./mobility_computations_$mode/