#!/bin/bash

####################################################################
## OVERVIEW  ##
# Script will facilitate the generation of a digital
# microstructure. There are several core sections.
# 1) Setup the environment 
# 2) Compile codes
# 3) Structure definition
# 4) Structure generation
# 5) Post processing
####################################################################

########################################################
## Setup environment                                  ##
########################################################
dir=${PWD}
bin=$dir/bin
src=$dir/src
includes=$dir/includes
tmp_keyword='tmpoutput_'`date '+%s'` # Date in seconds
mkdir /tmp/$tmp_keyword
tmpdir=/tmp/$tmp_keyword
modelfile=$tmpdir/rve.txt
if [ ! -d $tmpdir ]; then
    echo ERROR: tmp directory not created 
    exit
fi
time1=`date`
echo \# >> $modelfile
echo \# Generated: $time1 >> $modelfile
echo \# >> $modelfile


########################################################
## Install and compile                                ##
########################################################
echo Building CA executables
os=`uname -a | grep Darwin`

if [ "$os" > /dev/null ] ; then
    echo MAC os found $os
    make -f Makefile.Darwin
else
    make -f Makefile.gfortran
fi

echo Building MC
make -f Makefile.MC MC

if [ -x MC2vti ] ; then
    echo MC2vti up to date
else
    echo Building MC2vti
    g++ -o MC2vti MC2vti.cpp
fi

if [ -x 3d_enum ] ; then
    echo 3d_enum up to date
elif [ -e 3d_enumerate  ] ; then
    echo Building 3d_enum
    cd $dir/3d_enumerate
    make
    cp 3d_enum $dir/
else 
    echo No 3d_enum available
fi

cd $dir

if [ ! -f bimodal.pl -a bell_curve.pl -a recursiveSampler.pl ]; then 
  echo bimodal.pl or bell_curve.pl or recursiveSampler.pl does not exist in $dir 
  exit 1
fi

if [ ! -d $includes/voxelconvert ]; then
    echo voxelconvert not installed in $includes
    echo downloading voxelconvert
    cd $includes
    pwd
    svn co https://latir.materials.cmu.edu/svn/voxelconvert voxelconvert
fi

cd $includes/voxelconvert
make

cd $dir


#######################################################################
## Structure Definition                                              ##
#######################################################################

###########################################################
## Read and process commandline arguments                ##
###########################################################

# check for proper commandline arguments 

# if [[ $# -gt 1 ]] ; then
#     echo or     ./grow-3d.sh Inputfile number_of_iterations Periodic
#     exit 1
# elif [[ $# -eq 3 ]]; then
#     input_dir=`dirname $1`
#     input_file=$1
#     iterations=$2
#     Periodic=$3
# fi
# echo
# echo INPUT_DIR "    "$input_dir
# echo INPUT_FILE "   "$input_file
# echo ITERATIONS "   "$iterations
# echo PERIODIC "     "$Periodic

###########################################################
## Interactive                                           ##
###########################################################

echo
echo
echo "======Structure definition============================="
Eta=1
echo What is the scale of the model \(microns/voxel\)? [$Eta]
read eta
if [ "$eta" > /dev/null ] ; then
    Eta=$eta
fi
## Define dimension of the RVE in voxels
Vx=100
Vy=100
Vz=100
echo What are the dimensions of the RVE \(voxels\)? [$Vx $Vy $Vz]
read vx vy vz
if [ "$vx" > /dev/null ] ; then
Vx=$vx
fi
if [ "$vy" > /dev/null ] ; then
Vy=$vy
fi
if [ "$vz" > /dev/null ] ; then
Vz=$vz
fi

## Comput the size of the RVE in microns
Rx=`echo $Vx $eta|awk '{print $1*$2}'`
Ry=`echo $Vy $eta|awk '{print $1*$2}'`
Rz=`echo $Vz $eta|awk '{print $1*$2}'`

## Compute the unit_RVE
Vmax=$Vx
if [ "$Vy" -gt "$Vmax" ]; then 
Vmax=$Vy
fi
if [ "$Vz" -gt "$Vmax" ]; then 
Vmax=$Vz
fi

Rmax=`echo $Vmax $eta|awk '{print $1*$2}'`

Ux=`echo $Vx $Vmax|awk '{print $1/$2}'`
Uy=`echo $Vy $Vmax|awk '{print $1/$2}'`
Uz=`echo $Vz $Vmax|awk '{print $1/$2}'`

echo "Unit"
echo $Ux $Uy $Uz
echo
## Define RVE boundary conditions
Periodic=true
echo Are the boundary conditions periodic? true=1 false=2 [1]
read periodic

if [ "$periodic" > /dev/null ] ; then
    if [[ "$periodic" -eq 1 ]] ; then
	Preiodic=true
    elif [[ "$periodic" -ne 1 ]] ; then
	Periodic=false
    fi
fi

#######################################################################
## recursiveSampler.pl 
#
# Description: resursiveSampler.pl is a code that implements a
# recusive strategy to pack ellipsoids in to the unit RVE. It executes
# a perl script to generate a discrete set of ellipsoids and place
# them in the unit RVE. The script defines the number and size
# distribution of grains in the RVE. It monitors if ellipsoids overlap
# and recursively divides the unit RVE into smaller and smaller
# sections to guide placement of ellipsoids into unfilled regions.
#
# What does Fraction Mean? 
#
# (taken from recursiveSampler.pl) recursiveSampler.pl samples the
# computational space (either specified by -bbox, or the unit cube)
# randomly, and calls ellipsoidFunction at each of the sample points.
# it keeps track of the smallest dimension of the ellipsoid, and
# compares that with the dimension of it\'s bounding box.  if the
# bounding box is greater than the -fraction parameter (0.333 by
# default) times the smallest dimension of the ellipsoid, then
# recursiveSampler resamples the box with recursively subdivided
# boxes.
#######################################################################
##  Read ellipsoid overlap fraction
fraction=1.05
echo Please enter fraction of ellipsoid overlap. [$fraction]
read frac

if [ "$frac" > /dev/null ] ; then
fraction=$frac
fi

## Read ellipsoid distribution information and exectue recussiveSampler.pl
distfile="List"
ListFileName="ellipsoids.txt"
elscript="randomEntry.pl $ListFileName $Rx $Ry $Rz"
solution=1 
echo Ellipsoid Distribution file? ellipsoidList=1 bell_curve=2 bimodal=3 [$solution] 
read sol

if [ "$sol" > /dev/null ] ; then
solution=$sol
fi

if [ "$solution" -eq "3" ] ; then
    distfile="bimodal.pl"
    elscript="bimodal.pl"
    cp $elscript $tmpdir
elif [ "$solution" -eq "2" ] ; then
    distfile="bell_curve.pl"
    elscript="bell_curve.pl"
    cp $elscript $tmpdir
elif [ "$solution" -eq "1" ] ; then
    distfile="List"
    echo Enter the name of the ellipsoid list file. [$ListFileName]
    read listFileName
    if [ "$listFileName" > /dev/null ]; then
	ListFileName=$listFileName
    fi
    elscript="randomEntry.pl $ListFileName $Rx $Ry $Rz"
    cp $ListFileName $tmpdir
else
    echo ERROR: invalid response
    exit 1
fi

## Print the results of Structure Definition screen
echo
echo
echo "======Structure Summary============================="
echo Feature size: $a microns
echo Voxels for feature: $Va
echo RVE scale \(microns/voxel\): $eta
#echo RVE dimensions \(voxels\): $Vx $Vy $Vz
echo RVE dimensions \(microns\): $Rx $Ry $Rz
echo Unit RVE ratio:  $Ux $Uy $Uz
echo Periodic $Periodic
echo Fraction: $fraction
echo Ellipsoids: $distfile 
if [ "$solution" -eq "1" ] ; then
    echo EllipsoidList: $ListFileName
fi
echo 

## Print the results of Structure definition to model file
echo Scale\(microns/voxel\)  $Eta >> $modelfile
echo Voxels $Vx $Vy $Vz >> $modelfile
echo Periodic $Periodic >> $modelfile
echo Fraction $fraction >> $modelfile
echo Ellipsoids $distfile >> $modelfile
if [ "$solution" -eq "1" ] ; then
    echo EllipsoidList $ListFileName >> $modelfile
fi



#######################################################################
## Structure Generation                                              ##
#######################################################################

echo
echo
echo "======Structure generation============================="


###########################################################
## Execute recursiveSampler.pl to pack RVE               ##
###########################################################
echo =======executing recursiveSampler.pl with $distfile

if [ "$solution" -eq "3" ] ; then
    perl recursiveSampler.pl -bbox 0 $Ux 0 $Uy 0 $Uz -fraction $fraction 'perl bimodal.pl' > $tmpdir/test.input
elif [ "$solution" -eq "2" ] ; then
    perl recursiveSampler.pl -bbox 0 $Ux 0 $Uy 0 $Uz -fraction $fraction "perl bell_curve.pl $Rmax" > $tmpdir/test.input
elif [ "$solution" -eq "1" ] ; then
    perl recursiveSampler.pl -bbox 0 $Rx 0 $Ry 0 $Rz -fraction $fraction -list $ListFileName > $tmpdir/test.input
fi 

###########################################################
## Execute ellipticalFoam to choose a set of axtive      ##
## ellipsoids in the RVE to fill space and minimize      ##
## overlap.                                              ##
###########################################################
echo =======executing ellipticalFoam

if [ -f cell.ctrl ] ; then rm cell.ctrl ; fi
cd $tmpdir
ulimit
$dir/ellipticalFoam test.input
rm paddedPoints pointKey metaData annealing.log

###########################################################
## Grow ellipsoids using Cellular Automaton              ##                 
###########################################################
echo =======executing CA
if [ "$solution" -eq "3" ] ; then
    $dir/myCA  $Ux $Uy $Uz $Vx $Vy $Vz $Periodic
elif [ "$solution" -eq "2" ] ; then
    $dir/myCA  $Ux $Uy $Uz $Vx $Vy $Vz $Periodic
elif [ "$solution" -eq "1" ] ; then
    $dir/myCA  $Rx $Ry $Rz $Vx $Vy $Vz $Periodic
fi 



#######################################################################
## Post processing                                                   ##
#######################################################################
echo
echo
echo "======Post processing============================="

###########################################################
## Renumber as unique grains and set grain size threshold
###########################################################
$dir/includes/voxelconvert/voxelconvert ellipsoid_BaseCA.ph ellipsoid_BaseCA.mc

Threshold=0
echo How many voxels for the smallest grain? [$Threshold]
read threshold
if [[ "threshold" > /dev/null ]] ; then
    Threshold=$threshold
fi

$dir/includes/voxelconvert/ug2 ellipsoid_BaseCA.mc grains.mc $Threshold $Periodic
rm ellipsoid_BaseCA.mc
rm ellispoid_BaseCA.ph

###########################################################
## Output XML format for texturelist                     ##
###########################################################
$dir/includes/voxelconvert/voxel2XML grains.mc

###########################################################
## Visualization                                         ##
###########################################################
$dir/MC2vti grains.mc

###########################################################
## Stats                                                 ##
###########################################################
$dir/awk_stats.sh

###########################################################
## Grow grains to provide realism to structure           ##
###########################################################
iterations=0
echo
echo How many steps of MC grain growth? [$iterations]
read iters

if [ "$iters" > /dev/null ]; then
    iterations=$iters
fi

if [[ $iterations -gt 0 ]]; then

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
      if [ -d $tmpdir/${keyword2} ] ; then
	  keyword2=${grow_keyword}_v${count}
      else
	  break  
      fi
      let "count+=1"
    done

    grow_dir=${tmpdir}/${keyword2}
    mkdir ${grow_dir}

    if [ ! -d $grow_dir ]; then
	echo ERROR: growth directory not created 
	exit 
    fi
    $dir/MC ellipsoid_BaseMC.ph $grow_dir/grains_$grow_keyword.mc $iterations $Periodic
    mv $tmpdir/ellipsoid_BaseMC.ph $grow_dir/grains_$grow_keyword.ph

    cd $grow_dir

    ## Visualization
    Vti=1
    echo Generate VTI file? yes=1 no=2 [$Vti]
    read vti
    if [[ "$vti" > /dev/null ]] ; then
	Vti=$vti
    fi
    if [[ "$Vti" -eq 1 ]] ; then
	$dir/MC2vti ellipsoid_$grow_keyword.mc
    fi


    ## Stats
    Stats=1
    echo Run statistics analysis? yes=1 no=2 [$Stats]
    read stats

    if [[ "$stats" > /dev/null ]] ; then
	Stats=$stats
    fi

    if [[ "$Stats" -eq 1 ]] ; then
	$dir/stat3d ellipsoid_$grow_keyword.ph $Periodic
	$dir/awk_stats.sh
	$dir/3d_enum    
    fi
fi

##############################################################
## Move data files to unique folder to prevent overwrite    ##
##############################################################
#Safety check -- ensure directory with same filename does not exist
echo Moving data to OUTPUT_FILES directory

count=1
LIMIT=25

keyword='newoutput_'`date '+%m%d%y'`
keyword2=${keyword}_v${count}

while [ $count -lt $LIMIT ]
do
 if [ -d $dir/OUTPUT_FILES/${keyword2} ] ; then
   #echo satisfied for count=$count
   keyword2=${keyword}_v${count}
 else
   break  
 fi
 let "count+=1"
done

cd $dir
outputdir=OUTPUT_FILES/${keyword2}
mv $tmpdir $outputdir


echo Program Finished
echo Program started at $time1
echo Program ended at `date`
echo results are located in $dir/OUTPUT_FILES/${keyword2}
exit
