#!/bin/bash
cd /home/theochem/allouche/tmp/CChemI-211215/cchemi/examples/VPT2/QFFOrcaRIB2PLYP

node=`uname -n`

export LD_LIBRARY_PATH=/softs/openmpi-1.4.3-gcc/lib:/home/theochem/allouche/lib:/softs/openmpi-1.4.3-gcc/lib:/home/theochem/allouche/Softwares/ati-stream-sdk-v2.2-lnx32/lib/x86
export PATH=/softs/openmpi-1.4.3-gcc/bin:/home/theochem/allouche/shell:/home/theochem/allouche/Gabedit64:/home/theochem/allouche/Gabedit64:/softs/CompChemPackages/shell:/home/theochem/allouche/Softwares/bin:/home/theochem/allouche/bin:/home/theochem/allouche/shell:/home/theochem/allouche/Gabedit64:/home/theochem/allouche/Gabedit64:/softs/openmpi-1.4.3-gcc/bin:/home/theochem/allouche/Gabedit64:/home/theochem/allouche/Gabedit:/home/theochem/allouche/Gabedit:/softs/CompChemPackages/Turbomole/TURBOMOLE/bin/x86_64-unknown-linux-gnu:/softs/CompChemPackages/Turbomole/TURBOMOLE/scripts:/softs/CompChemPackages/Turbomole/TURBOMOLE/../bin:/usr/local/lam-mpi/bin/:/home/theochem/allouche/bin:/usr/local/bin:/usr/bin:/bin:/usr/games:.:.:/softs/CompChemPackages/shell:/softs/CompChemPackages/bin:/softs/CompChemPackages/Gabedit:/softs/CompChemPackages/Gaussian09D/g09/bsd:/softs/CompChemPackages/Gaussian09D/g09/local:/softs/CompChemPackages/Gaussian09D/g09/extras:/softs/CompChemPackages/Gaussian09D/g09:/home/admin/CompChemPackages/deMonNano-nosrc/bin:/softs/CompChemPackages/DFTBP/bin:/softs/CompChemPackages/nbo6/bin:/softs/CompChemPackages/CP2K:.:/home/theochem/allouche/logiciels/tlse/shell:/home/theochem/allouche/logiciels/tlse/bin:/home/theochem/allouche/shell:/home/theochem/allouche/psi3/bin
xMakeQFF2MRVPT2 > "/home/theochem/allouche/tmp/CChemI-211215/cchemi/examples/VPT2/QFFOrcaRIB2PLYP/xMakeQFF2MRVPT2.txt"

