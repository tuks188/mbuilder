# Makefile
# Minimal GNU makefile for programs using MMSP
# Questions/comments to ssintay@andrew.cmu.edu (Stephen Sintay)

# compiler/flags
compiler = g++
flags = -g -ansi 

# dependencies
core = MMSP.data.hpp MMSP.grid.hpp MCgrid.hpp voxelio.hpp tokenize.hpp

# single processor program
all : voxel2XML voxelconvert

voxelconvert : voxelconvert.cpp $(core)
	$(compiler) $(flags) $< -o $@

voxel2XML : voxel2XML.cpp $(core)
	$(compiler) $(flags) $< -o $@

clean : 
	rm voxel2XML voxelconvert 

archive : 
	tar -cvf voxelconvertPackage.tar Makefile $(core) voxel2XML.cpp voxelconvert.cpp

