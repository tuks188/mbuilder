#!/bin/bash

echo number of arguments $#
dir=`pwd`
outputdir="work"
for ((i=1;i<=$#;i+=1)); do
var=${!i}
#echo $var
f=$(basename $var)
cp $dir/$f/Dist_Facets_by_size_ellipsoid_.txt $dir/$output/Dist_Facets_by_size_ellipsoid_${f}.txt
cp $dir/$f/Dist_freq.ellipsoid_.txt $dir/$output/Dist_freq.ellipsoid_${f}.txt
cp $dir/$f/Dist_PDF.ellipsoid_.txt $dir/$output/Dist_PDF.ellipsoid_${f}.txt
cp $dir/$f/NN.ellipsoid_*.ph.txt $dir/$output/NN.ellipsoid_xx.ph_${f}.txt
done


