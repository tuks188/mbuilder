// MCgrid.hpp
// Class definitions for rectilinear Monte Carlo grids
// Questions/comments to jgruber@andrew.cmu.edu (Jason Gruber)

#ifndef MCGRID
#define MCGRID
#include"MMSP.grid.hpp"
#include<string.h>

namespace MMSP{

class MCgrid1D : public grid1D<scalar,int> {
public:
        // constructors
        MCgrid1D(int x)
                : grid1D<scalar,int>(x) {}
        MCgrid1D(const char* filename)
                : grid1D<scalar,int>(filename) {}
        MCgrid1D(const char* filename, int id, int np, int ng = 1)
                : grid1D<scalar,int>(filename,id,np,ng) {}

        // the update function
        void update(int steps = 1);
        void update(int steps, int id, int np, int ng = 1);

        // utility functions
        std::vector<int> nonzero(int x) const;
};

std::vector<int> MCgrid1D::nonzero(int x) const
{
        std::vector<int> nonzero;
        for (int i=-1; i<=1; i++)
                //if (abs(i)<=1)
                        nonzero.push_back(neighbor(x,i));
        sort(nonzero.begin(),nonzero.end());
        nonzero.erase(unique(nonzero.begin(),nonzero.end()),nonzero.end());
        return nonzero;
}


class MCgrid2D : public grid2D<scalar,int> {
public:
        // constructors
        MCgrid2D(int x, int y)
                : grid2D<scalar,int>(x,y) {}
        MCgrid2D(const char* filename)
                : grid2D<scalar,int>(filename) {}
        MCgrid2D(const char* filename, int id, int np, int ng = 1)
                : grid2D<scalar,int>(filename,id,np,ng) {}

        // the update function
        void update(int steps = 1);
        void update(int steps, int id, int np, int ng = 1);

        // utility functions
        std::vector<int> nonzero(int x, int y) const;
};

std::vector<int> MCgrid2D::nonzero(int x, int y) const
{
        std::vector<int> nonzero;
        for (int i=-1; i<=1; i++)
                for (int j=-1; j<=1; j++)
                        //if (abs(i)+abs(j)<=1) 
                                nonzero.push_back(neighbor(x,y,i,j));
        sort(nonzero.begin(),nonzero.end());
        nonzero.erase(unique(nonzero.begin(),nonzero.end()),nonzero.end());
        return nonzero;
}


class MCgrid3D : public grid3D<scalar,int> {
public:
        // constructors
        MCgrid3D(int x, int y, int z)
                : grid3D<scalar,int>(x,y,z) {}
        MCgrid3D(const char* filename)
                : grid3D<scalar,int>(filename) {}
        MCgrid3D(const char* filename, int id, int np, int ng = 1)
                : grid3D<scalar,int>(filename,id,np,ng) {}

        // the update function
        void update(int steps = 1);
        void update(int steps, int id, int np, int ng = 1);

        // utility functions
        std::vector<int> nonzero(int x, int y, int z) const;
};

std::vector<int> MCgrid3D::nonzero(int x, int y, int z) const
{
        std::vector<int> nonzero;
        for (int i=-1; i<=1; i++)
                for (int j=-1; j<=1; j++)
                        for (int k=-1; k<=1; k++)
                                //if (abs(i)+abs(j)+abs(k)<=1)
                                        nonzero.push_back(neighbor(x,y,z,i,j,k));
        sort(nonzero.begin(),nonzero.end());
        nonzero.erase(unique(nonzero.begin(),nonzero.end()),nonzero.end());
        return nonzero;
}

} // namespace MMSP

#endif
