#!/bin/bash

# AUTHOR: Chris Roberts 4/07
####################################################################
## OVERVIEW  ##
# Script will facilitate the generation of a digital microstructure
# and texture based on experimental data.
####################################################################

dir=${PWD}
time1=`date`
########################################################
## FILE CHECK                                         ##
########################################################
if [ -x annealfinal -a -x RodToE -a -x stat3d -a -x myCA \
     -a -x odfextract_final -a -x mdfextract_final -a -x ellipticalFoam \
     -a -x wts2ang -a -x MC_mdf_check ]; then
  echo All necessary executables have been built.
else
  echo you are missing an executable from the package
  echo type " make " upon exiting this script
  exit 1
fi
 

###########################################################
## Calculate ellipsoid statistical features              ##
## Create an XML file for texture-fitting programs       ##
###########################################################
echo Extracting ellipsoid statistics from file
echo What is filename?
read inputfile
echo $inputfile will be renamed ellipsoid.yourCA.ph
cp $inputfile ellipsoid.yourCA.ph
#./stat3d ellipsoid.yourCA.ph 
./stat3d2 ellipsoid.yourCA.ph 
#rm size*yourCA aspect*yourCA
mv cellIdealization2.xml cellIdealization.xml

awk '{if($1 ~/id/)id=$2}{if($1~/volume/)print id,$2^(0.333333333)}' \
cellIdealization.xml > gsd.txt


#####################################################################
## Determine ODF and MDF based on ANG and GF#1 files (experiment)  ##
#####################################################################

cp symop.txt odfextract_final mdfextract_final $dir/ANG_FILE_STORAGE
if [ $? -ne 0 ] ; then  echo Problem with odf \& mdf  programs; exit 1 ; fi
cd $dir/ANG_FILE_STORAGE

## Directory CLEANUP
if [ -f mdfoutputlist ] ; then rm mdfoutputlist ; fi
if [ -f odfoutputlist ] ; then rm odfoutputlist ; fi
if [ -f mdfoutname.txt ] ; then rm mdfoutname.txt ; fi
if [ -f odfoutname.txt ] ; then rm odfoutname.txt ; fi
if [ -f ctrlfile.txt ] ; then  rm ctrlfile.txt ; fi

# Extract MDF in Homochoric space from each data file listed in inputlist

ls *.gf1.txt  *.ang > texlist

wc -l texlist | awk '{print $1}'  > ctrlfile.txt

while read fname
do
   echo $fname >> ctrlfile.txt
   echo 1.0 0.0 0.0 >> ctrlfile.txt
   echo 0.0 1.0 0.0 >> ctrlfile.txt
   echo 0.0 0.0 1.0 >> ctrlfile.txt
done < texlist

rm texlist

./mdfextract_final
./odfextract_final

cp evodf.txt evmdf.txt $dir
cd $dir

#######################################################################
## FIT ODF and MDF to digital microstucture using Stochastic Anneal  ##
#######################################################################
echo Fitting texture data to digital microstructure
echo To view progress, open extra window
echo cd $dir and type 'tail -n 50 -f anneal.log'
./annealfinal

###################################################
## Orientation Conversion and Output             ##
###################################################
echo Converting texture file from Rodrigues space to Euler space
./RodToE                # Rodrigues to Euler space
./RodToE_wts            # creates wts file
./rod2eul               # creates wts file for rex3d program
cp fridyTexin texin1

# generate MDF and WTS file from digital MS and texture
#./MC_mdf_check  ellipsoid.yourCA.ph
./xml_texture cellIdealization.xml 
# converts WTS file to an ANG file for TSL software
./wts2ang  ellipsoid.yourCA.odf ellipsoid.yourCA.ang

##############################################################
## Move data files to unique folder to prevent overwrite    ##
##############################################################
#Safety check -- ensure directory with same filename does not exist
echo Moving data to OUTPUT_FILES directory

count=2
LIMIT=25

keyword='newoutput_'`date '+%m%d%y'`
keyword2=$keyword

while [ $count -lt $LIMIT ]
do
 if [ -d $dir/OUTPUT_FILES/${keyword2} ] ; then
   #echo satisfied for count=$count
   keyword2=${keyword}_v${count}
 else
   echo new directory will be named ${keyword2} and will be located
   echo in the $dir/OUTPUT_FILES directory.
   break  
 fi
 let "count+=1"
done

mkdir OUTPUT_FILES/${keyword2}

mv evodf.txt evmdf.txt ellipsoid.yourCA*  EAorts.txt anneal.log \
 orts.txt orts.wts fridyTexin rmdf.txt rodf.txt texin1  annealing.log \
 cellIdealization.xml  pointKey paddedPoints anneal.log active.list \
 metaData test.input elliptical.cells size*yourCA aspect*yourCA \
 OUTPUT_FILES/${keyword2}

echo Program Finished
echo All files were moved to OUTPUT_FILES/${keyword2}

echo Program started at $time1
echo Program ended at `date`
exit
