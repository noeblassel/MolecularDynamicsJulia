#!/bin/bash

for rhomin in 0.4502 0.45771249999999997 0.465225 0.4727375 0.48025 0.4877625 0.495275 0.5027874999999999 0.5103 0.5178125 0.525325 0.5328375 0.54035 0.5478625 0.555375 0.5628875 0.5704
do
nohup julia test_nist.jl $rhomin > nohup.out &

done