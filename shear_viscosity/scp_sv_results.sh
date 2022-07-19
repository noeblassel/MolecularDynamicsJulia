#!/bin/bash

for node in 15 16 18
do
echo clustern$node
scp clustern$node:/libre/blasseln/MolecularDynamicsJulia/shear_viscosity/*.out .
done