#!/bin/bash

awk '{if($1 ~/\<nregions\>/)grains=$2} \
{if($1 ~/\<id\>/)id=$2}\
{if($1 ~/\<id\>/)ids[id]=$2} \
{if($1~/volume/) vols[id] = $2} \
{if($1~/volume/) rads[id] = $2^(0.33333333)} \
{if($1~/volume/) vol_sum+=$2} \
{if($1~/volume/) rad_sum+=$2^(0.33333333)} \
{if($1~/npatches/) neighs[id] = $2} \
{if($1~/npatches/) neigh_sum += $2} \
END { \
print "ID \t Vol/<Vol> \t Neigh/<Neigh> \t Rad/<Rad> \t Vol \t Neigh \t Rad"
    for (x = 0; x <= grains-1; x++) \
    print ids[x] "\t" \
          vols[x]/vol_sum "\t" \
          neighs[x]/neigh_sum "\t" \
          rads[x]/rad_sum "\t" \
          vols[x] "\t" \
          neighs[x] "\t" \
          rads[x]  \
}' grains.xml > grain_stats.txt




