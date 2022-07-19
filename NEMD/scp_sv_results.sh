#!/bin/bash
for node in 15 16 18 do 
scp clustern$node:/libre/blasseln/MolecularDynamicsJulia/NEMD/*.out .
done