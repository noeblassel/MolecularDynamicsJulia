#!/bin/bash

ssh clustern15 "cd /libre/blasseln/MolecularDynamicsJulia/Langevin/ ; cat BAO*>BAO.csv"
scp clustern15:/libre/blasseln/MolecularDynamicsJulia/Langevin/BAO.csv ./bias_dumps/BAO.csv
ssh clustern15 "cd /libre/blasseln/MolecularDynamicsJulia/Langevin/ ; rm BAO.csv"

ssh clustern16 "cd /libre/blasseln/MolecularDynamicsJulia/Langevin/ ; cat BAOAB*>BAOAB.csv"
scp clustern16:/libre/blasseln/MolecularDynamicsJulia/Langevin/BAOAB.csv ./bias_dumps/BAOAB.csv
ssh clustern16 "cd /libre/blasseln/MolecularDynamicsJulia/Langevin/ ; rm BAOAB.csv"

ssh clustern17 "cd /libre/blasseln/MolecularDynamicsJulia/Langevin/ ; cat BABO*>BABO.csv"
scp clustern17:/libre/blasseln/MolecularDynamicsJulia/Langevin/BABO.csv ./bias_dumps/BABO.csv
ssh clustern17 "cd /libre/blasseln/MolecularDynamicsJulia/Langevin/ ; rm BABO.csv"

ssh clustern18 "cd /libre/blasseln/MolecularDynamicsJulia/Langevin/ ; cat BAOA*>BAOA.csv"
scp clustern18:/libre/blasseln/MolecularDynamicsJulia/Langevin/BAOA.csv ./bias_dumps/BAOA.csv
ssh clustern18 "cd /libre/blasseln/MolecularDynamicsJulia/Langevin/ ; rm BAOA.csv"