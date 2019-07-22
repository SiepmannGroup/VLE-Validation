#!/bin/bash
# set the number of nodes and processes per node
#SBATCH --nodes=1

# set the number of tasks (processes) per node.
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1

# set max wallclock time
#SBATCH --time=120:00:00

# set name of job
#SBATCH --job-name=C5OH_510_1

#SBATCH --partition=E5-26xx
#SBATCH --output test.out


module purge
module load intel ompi
# "PROG" points to the executable
PROG=~/exe-06-2019/src/topmon
FILELOG=log

WRKDIR=$PWD
cd $WRKDIR

# "FCUR" loops over a sequence of strings representing different phases of the simulation.
# If there is a corresponding fort.4.$FCUR file, it will be copied to fort.4 to be used as
# the input file; otherwise the existing fort.4 file will be used.
for FCUR in prod;do
    if [ -e "run.$FCUR" ]; then
        echo "$FCUR already finished"
        continue
    fi


    # echo "$FCUR: " >> $FILELOG
    [ -e "fort.4.$FCUR" ] && cp -f "fort.4.$FCUR" fort.4

    $PROG
    #/usr/bin/time -ao $FILELOG mpirun -np 4 $PROG || exit -1
    if grep 'ERROR\|error\|fault\|severe' test.out
    then
        echo "$FCUR crashed"
        break
    else
        echo "$FCUR finished"
    fi
    cp -f config1a.dat fort.77
    mv -f run1a.dat "run.$FCUR"
    mv -f config1a.dat "config.$FCUR"
    mv -f movie1a.dat "movie.$FCUR"
    mv -f fort.12 "fort12.$FCUR"
    for j in 1 2 ;do
        mv -f "box${j}config1a.xyz" "box${j}config.$FCUR"
        mv -f "box${j}movie1a.xyz" "box${j}movie.$FCUR"
    done
done

