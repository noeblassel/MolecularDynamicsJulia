#!/bin/bash


#scp clustern14:/libre/blasseln/MolecularDynamicsJulia/Langevin/BABO* clustern17:/libre/blasseln/MolecularDynamicsJulia/Langevin/
#scp clustern14:/libre/blasseln/MolecularDynamicsJulia/Langevin/BAOA* clustern18:/libre/blasseln/MolecularDynamicsJulia/Langevin/
#scp clustern16:/libre/blasseln/MolecularDynamicsJulia/Langevin/BAO* clustern15:/libre/blasseln/MolecularDynamicsJulia/Langevin/

ssh clustern23 "cd /libre/blasseln/MolecularDynamicsJulia/Langevin/ ; cat GHMC0*>GHMC.csv"
scp clustern23:/libre/blasseln/MolecularDynamicsJulia/Langevin/GHMC.csv ./bias_dumps/GHMC.csv
ssh clustern23 "cd /libre/blasseln/MolecularDynamicsJulia/Langevin/ ; rm GHMC.csv"

ssh clustern23 "cd /libre/blasseln/MolecularDynamicsJulia/Langevin/ ; cat BAOAB0*>BAOAB.csv"
scp clustern23:/libre/blasseln/MolecularDynamicsJulia/Langevin/BAOAB.csv ./bias_dumps/BAOAB.csv
ssh clustern23 "cd /libre/blasseln/MolecularDynamicsJulia/Langevin/ ; rm BAOAB.csv"

ssh clustern25 "cd /libre/blasseln/MolecularDynamicsJulia/Langevin/ ; cat BABO0*>BABO.csv"
scp clustern25:/libre/blasseln/MolecularDynamicsJulia/Langevin/BABO.csv ./bias_dumps/BABO.csv
ssh clustern25 "cd /libre/blasseln/MolecularDynamicsJulia/Langevin/ ; rm BABO.csv"

ssh clustern25 "cd /libre/blasseln/MolecularDynamicsJulia/Langevin/ ; cat BAOA0*>BAOA.csv"
scp clustern25:/libre/blasseln/MolecularDynamicsJulia/Langevin/BAOA.csv ./bias_dumps/BAOA.csv
ssh clustern25 "cd /libre/blasseln/MolecularDynamicsJulia/Langevin/ ; rm BAOA.csv"

ssh clustern25 "cd /libre/blasseln/MolecularDynamicsJulia/Langevin/ ; cat BAO0*>BAO.csv"
scp clustern25:/libre/blasseln/MolecularDynamicsJulia/Langevin/BAO.csv ./bias_dumps/BAO.csv
ssh clustern25 "cd /libre/blasseln/MolecularDynamicsJulia/Langevin/ ; rm BAO.csv"