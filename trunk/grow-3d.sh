#!/bin/bash

# AUTHOR: Steve Sintay 10/10
####################################################################
## OVERVIEW  ##
# Script will facilitate the generation of a digital microstructure
# and texture based on experimental data.
####################################################################

dir=${PWD}

########################################################
## Compile and Install                                ##
########################################################
echo Building MC
make -f Makefile.MC MC

if [ -x MC2vti ] ; then
    echo MC2vti up to date
else
    echo Building MC2vti
    g++ -o MC2vti MC2vti.cpp
fi

echo Building 3d_enum
cd $dir/3d_enumerate
make
cp 3d_enum $dir/

cd $dir

###########################################################
## Read and process commandline arguments                ##
###########################################################

# check for proper commandline arguments 

if [[ $# -lt 2 ]] ; then
    echo or     ./grow-3d.sh Inputfile number_of_iterations Periodic
    exit 1
elif [[ $# -eq 3 ]]; then
    input_dir=`dirname $1`
    input_file=$1
    iterations=$2
    Periodic=$3
fi
echo
echo INPUT_DIR "    "$input_dir
echo INPUT_FILE "   "$input_file
echo ITERATIONS "   "$iterations
echo PERIODIC "     "$Periodic




###########################################################
## Grow grains to provide realism to structure           ##
###########################################################
echo Executing MC growth on structure
if [[ $iterations -lt 10 ]] ; then
    grow_keyword='grw0'$iterations
elif [[  $iterations -ge 10 ]] ; then
    grow_keyword='grw'$iterations
fi
count=1
LIMIT=25

keyword2=${grow_keyword}_v${count}

while [ $count -lt $LIMIT ]
do
 if [ -d $input_dir/${keyword2} ] ; then
   #echo satisfied for count=$count
   keyword2=${grow_keyword}_v${count}
 else
   echo new mc growth direcory will be named ${keyword2} and will be located
   echo in the $input_dir directory.
   break  
 fi
 let "count+=1"
done

grow_dir=$input_dir/${keyword2}
mkdir $input_dir/${keyword2}

if [ ! -d $grow_dir ]; then
    echo ERROR: growth directory not created 
    exit 
fi

$dir/MC $input_file $grow_dir/ellipsoid_$grow_keyword.mc $iterations $Periodic
mv mc_rewrite.ph $grow_dir/ellipsoid_$grow_keyword.ph

cd $grow_dir

## Visualization
Vti=1
echo Generate VTI file? yes=1 no=2 [1]
read vti
if [[ "$vti" > /dev/null ]] ; then
    Vti=$vti
fi
if [[ "$Vti" -eq 1 ]] ; then
    $dir/MC2vti ellipsoid_$grow_keyword.mc
fi


## Stats
Stats=1
echo Run statistics analysis? yes=1 no=2 [1]
read stats

if [[ "$stats" > /dev/null ]] ; then
    Stats=$stats
fi

if [[ "$Stats" -eq 1 ]] ; then
    $dir/stat3d ellipsoid_BaseMC.ph $Periodic
    $dir/awk_stats.sh
    $dir/3d_enum    
fi

