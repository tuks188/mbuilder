// MCgrid.example.hpp
// Isotropic coarsening algorithm for 2D and 3D Monte Carlo methods
// Questions/comments to jgruber@andrew.cmu.edu (Jason Gruber)

#ifndef MCGRID_UPDATE
#define MCGRID_UPDATE
#include"MCgrid.hpp"

namespace MMSP{

void MCgrid2D::update(int steps)
{
	MCgrid2D& grid = *this;
	int n = nx[0]*nx[1];
	const float kT = 0.75;

	for (int step=0; step<steps; step++) {
		for (int h=0; h<n; h++) {
			int x = rand()%nx[0];
			int y = rand()%nx[1];
			int spin1 = grid[x][y];
			
			std::vector<int> neighbors = nonzero(x,y);
			int spin2 = neighbors[rand()%neighbors.size()];

			if (spin1!=spin2) {
				float dE = -1.0;
				for (int i=-1; i<=1; i++)
					for (int j=-1; j<=1; j++) {
						int spin = grid.neighbor(x,y,i,j);
						dE += (spin!=spin2)-(spin!=spin1);
					}
				float num = static_cast<float>(rand())/static_cast<float>(RAND_MAX);
				if (dE<=0.0) grid[x][y] = spin2;
				else if (dE>0.0 && num<exp(-dE/kT)) grid[x][y] = spin2;
			}
		}
	}
}

void MCgrid3D::update(int steps)
{
	MCgrid3D& grid = *this;
	int n = nx[0]*nx[1]*nx[2];
	//const float kT = 1.5;
	const float kT =1.5;
	int bridge_count=0;

	for (int step=0; step<steps; step++) {
		for (int h=0; h<n; h++) {
			int x = rand()%nx[0];
			int y = rand()%nx[1];
			int z = rand()%nx[2];
			int spin1 = grid[x][y][z];
			
			// Choose spin2 from all NN
			std::vector<int> neighbors = nonzero(x,y,z);
			int spin2 = neighbors[rand()%neighbors.size()];

			// Find orthogonal NN for spin1
			std::vector<int> orth_neighs;
			orth_neighs.push_back(grid.neighbor(x,y,z,1,0,0));
			orth_neighs.push_back(grid.neighbor(x,y,z,0,1,0));
			orth_neighs.push_back(grid.neighbor(x,y,z,0,0,1));
			orth_neighs.push_back(grid.neighbor(x,y,z,-1,0,0));
			orth_neighs.push_back(grid.neighbor(x,y,z,0,-1,0));
			orth_neighs.push_back(grid.neighbor(x,y,z,0,0,-1));

			// Chose spin2 from orthogonal NN only
			//int candidate = rand()%orth_neighs.size();
			//int spin2 = orth_neighs[candidate];

			// Bypass if the spins are the same
			if (spin1!=spin2){

			  int like_count=0;
			  std::vector<int> like_nn;
			  bool united=true;

			  
			  // Check to see how many orthogonal NN are the same spin
			  for(int ii=0; ii < orth_neighs.size(); ii++){
				if(orth_neighs[ii] == spin1){
				  like_count++;
				}
			  }

			  // Check to see if we have a case where spin1 is a
			  // bridge between two and only two like voxels
			  if(like_count == 2){
				for(int ii=0; ii < orth_neighs.size(); ii++){
				  if(orth_neighs[ii] == spin1){
					if( ii > 2 ){
					  if( orth_neighs[ii-3] == spin1){
						united = false;
						bridge_count++;
					  }
					}else if( ii < 2 ){
					  if( orth_neighs[ii+3] == spin1){
						united = false;
						bridge_count++;
					  }
					}
				  }
				}
			  }
			  
			  if(united){
				float dE = -1.0;
				for (int i=-1; i<=1; i++)
				  for (int j=-1; j<=1; j++)
					for (int k=-1; k<=1; k++) {
					  int spin = grid.neighbor(x,y,z,i,j,k);
					  dE += (spin!=spin2)-(spin!=spin1);
					}
				float num = static_cast<float>(rand())/static_cast<float>(RAND_MAX);
				if (dE<=0.0) grid[x][y][z] = spin2;
				else if (dE>0.0 && num<exp(-dE/kT)) grid[x][y][z] = spin2;
			  }
			}
		}
	}
	std::cout << "Bridge attempts: " << bridge_count << std::endl;
}

} // namespace MMSP

#endif
