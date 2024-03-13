#!/bin/bash
export LIBCCHEMIDIR=/home/allouche/MySoftwares/CChemI/CChemISRC/cchemilib
dirName="$PWD/OpenBabel232/lib"
export OBLIBDIR=$dirName
fileName="$dirName/CONFIG"
echo "OBLIBDIR=$dirName" > $fileName
echo "RM = rm -f">> $fileName
echo "MAKE = make">>$fileName
echo "CCPP = g++ -fPIC -I \$(OBLIBDIR)/include">>$fileName
echo "COMMONCFLAGS = -Wall -O2 \$(CLCFLAGS)">>$fileName

make clean
