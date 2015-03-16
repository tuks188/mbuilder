# CMU Microstructure Builder (mbuilder) #

Microstructure builder or **mbuilder** is a strategy to construct simulated 3D polycrystalline materials. The input is typically grain size and shape data as obtained from orthogonal images (optical or SEM). The output is a 3D voxel structure that matches the size and shape statistics provided at input.

The voxel structures can be used directly as input to Monte Carlo simulations or can be converted to mesh structures for use in FE structural analysis.

## History of **mbuilder** ##

Microstructure Builder started as a collaboration between David Saylor at Carnegie Mellon University (CMU) and Joe Fridy at the Alcoa Technical Center, with help from Tony Rollett, Bassem El-Dasher and Kyung-Jun Kung (all at CMU).  It was supported by the Mesoscale Interface Mapping Project or MIMP under the NSF-supported Materials Research Science and Engineering Center at CMU (mimp.materials.cmu.edu).  Various individuals have contributed to **mbuilder** over the years, including Chris Roberts, Abhijit Brahme, Sukbin Lee and Steve Sintay.  Programs that have supported it include the Computational Materials Science Network (CMSN), DARPA under the SIPS program and the Commonwealth of Pennsylvania.

## Getting started ##

**mbuilder** is a collection of various tools that are managed through the shell script, [generate-3d.sh](http://code.google.com/p/mbuilder/source/browse/trunk/generate-3d.sh). This is an all purpose code that will do the following:
  1. Download additional, necessary tools from a server at CMU using SVN
  1. Compile source
  1. Interactively guide the user in microstructure generation

### Additional tools ###
The voxelconvert code is downloaded from a server at CMU using SVN. This tool allows for file format management between the various steps of the microstructure generation process. It also provides post processing ability that allows the user to specify a size threshold on the smallest grains that are represented.

### Compile source ###
This section will likely need to be modified for each user depending on the compilers available. The code will attempt to determine if a user is running a OSX operating system and will call a different make file in this case.

### Interactive microstructure generation ###
This is a step by step process that guide the user through the following steps
  1. Define the scale of the voxels in the microstructure. (i.e microns/voxel)
  1. Define the number of voxels in the x,y,z directions. When combined with the previous step this determines the overall size and shape of the generated volume element.
  1. Define the boundary conditions as periodic or non-periodic
  1. Define the ellipsoid distribution for grain representation and some parameters for packing the ellipsoids into the generated volume element.
  1. Perform microstructure generation to produce a voxel filled volume with unique regions segmented as grains.
  1. Perform post processing on this microstructure.

<a href='Hidden comment: 
<img height=400 src=http://neon.mems.cmu.edu/rollett/ADR_CMU_organ_19Aug06.jpg>
'></a>