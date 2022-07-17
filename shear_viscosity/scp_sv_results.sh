#!/bin/bash

for node in 15 16 19 20 24 25 
do
scp clustern$node:/libre/blasseln/MolecularDynamicsJulia/shear_viscosity/*.out .
done