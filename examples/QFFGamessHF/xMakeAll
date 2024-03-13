#!/bin/bash
echo " Begin QFF2MR-VPT2 calculation Hybrid B3LYP/PM6 ">  $0.log

# SCRIPT=$(readlink -f $0) works only under linux

# Works on both linux and Mac OS X.
function read_link() {
    local path=$1
        if [ -d $path ] ; then
        local abspath=$(cd $path; pwd)
        else
            local prefix=$(cd $(dirname -- $path) ; pwd)
            local suffix=$(basename $path)
            local abspath="$prefix/$suffix"
        fi
    echo $abspath
}
SCRIPT=$(read_link $0)

export  DEFDIR=`dirname $SCRIPT`
echo "DEFDIR=$DEFDIR"
cd $DEFDIR

export IGVPT2DIR=$HOME/iGVPT2
export PATH=$IGVPT2DIR/bin:$PATH
export PATH=$IGVPT2DIR/shell:$PATH

rm $0.err

# Opt + Freq => h2oOptFreq, runGamess
#To show the first derivatives of the dipole, you need to make an input file like this :
#First you ask an optimisation of the geometry
#Then you ask to read the optimized geometry and compute the frequencies. Do a two steps job with --Link1--.
#For the second step, the frequencies, you need to add the iop : iop(7/33=1)
#runGamess h2oOptFreq.inp >$0.log 2>  $0.err
#rm -f $0.err

# generate QFF energies files
/bin/rm -f h2oVQFF_*.ici
/bin/rm -f h2oVQFF_*.gab
/bin/rm -f h2oVQFF_*.log
/bin/rm -f *FirstDerivatives.txt
igvpt2 h2oV.ici >> $0.log 2>>  $0.err

#echo "cd $DEFDIR" > xp
#ls h2oVQFF_* >> xp
ls h2oVQFF_* > xp

sed -i 's/h2oVQFF/runiGVPT2 h2oVQFF/' xp
#cat xp
chmod u+x xp

#parallel

# SLURM
#module purge
#module use /softs/modulefiles
#module load gcc/4.7.3 openmpi/1.7.4/openmpi-gcc47
#mpirun parallel ./xp >>  $0.log 2>> $0.err

#PBS
#export LD_LIBRARY_PATH=/softs/openmpi-1.4.3-gcc/lib:$LD_LIBRARY_PATH
#export PATH=/softs/openmpi-1.4.3-gcc/bin:$PATH
#export PATH=$HOME/bin:$PATH
#mpiexec -machinefile $PBS_NODEFILE parallel xp >>  $0.log 2>> $0.err

#single 
./xp >>  $0.log 2>> $0.err

sed -i 's/RunType=GenerateQFFnMRFiles/#RunType=GenerateQFFnMRFiles/g' h2oV.ici
sed -i 's/#RunType=ComputeQFFnMRFromFiles/RunType=ComputeQFFnMRFromFiles/g' h2oV.ici

# generate txt data for VPT2 calculation
igvpt2 h2oV.ici  >>  $0.log 2>> $0.err

sed -i 's/#RunType=GenerateQFFnMRFiles/RunType=GenerateQFFnMRFiles/g' h2oV.ici
sed -i 's/RunType=ComputeQFFnMRFromFiles/#RunType=ComputeQFFnMRFromFiles/g' h2oV.ici

# make ici  file for VPT2 calculation
echo 'RunType=VPT2' > h2oVQFF.ici
#cat h2oVFirstDerivatives.txt >> h2oVQFF.ici
cat h2oVQFF.txt >> h2oVQFF.ici

# generate txt data for VPT2 calculation
igvpt2 h2oVQFF.ici  >  $0.out 2>> $0.err


#/bin/rm -f xp
/bin/rm -f xp_*
#/bin/rm -f h2oVQFF_*.ici
#/bin/rm -f h2oVQFF_*.gab
#/bin/rm -f h2oVQFF_*.log
echo "==========================================================================="
echo "See $0.out , $0.log and $0.err files "
echo "==========================================================================="
