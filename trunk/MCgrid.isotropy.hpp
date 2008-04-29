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
	const float kT = 1.5;

	for (int step=0; step<steps; step++) {
		for (int h=0; h<n; h++) {
			int x = rand()%nx[0];
			int y = rand()%nx[1];
			int z = rand()%nx[2];
			int spin1 = grid[x][y][z];
			
			std::vector<int> neighbors = nonzero(x,y,z);
			int spin2 = neighbors[rand()%neighbors.size()];

			if (spin1!=spin2) {
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

} // namespace MMSP

#endif
