#!/bin/bash

for n in 14 16 17 18
do
scp clustern$n:/libre/blasseln/MolecularDynamicsJulia/BAOA_tests/*.out .
done