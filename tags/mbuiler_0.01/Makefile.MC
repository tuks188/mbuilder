# Makefile.MC
# Minimal GNU makefile for Monte Carlo programs using MMSP
# Questions/comments to jgruber@andrew.cmu.edu (Jason Gruber)

# compiler/flags
compiler = g++
flags = -O3 -ansi -g
pcompiler = mpic++
pflags = -O3 -ansi

# dependencies
core = MMSP.data.hpp MMSP.grid.hpp MCgrid.hpp
single = $(core) MCgrid.update.hpp 
parallel = $(core) MCgrid.update.parallel.hpp

# single processor program
MC : MCgrid.cpp $(single)
	$(compiler) $(flags) $< -o $@

# parallel processor program
MC.parallel : MCgrid.parallel.cpp $(parallel)
	$(pcompiler) $(pflags) $< -o $@

