// MMSP.grid.hpp
// Class definitions for rectilinear grids (base class)
// Questions/comments to jgruber@andrew.cmu.edu (Jason Gruber)

#ifndef MMSP_GRID
#define MMSP_GRID
#include"MMSP.data.hpp"
#include<iostream>
#include<fstream>
#include<cstdarg>

namespace MMSP{

template <int dimension, typename type>
class grid{
public:
	// constructors
	grid(int nf, ...);
	grid(const char* filename) {this->input(filename);}
	grid(const char* filename, int id, int np, int ng = 1);

	// subscript operator
	type& operator[](int i) {return data[i];}
	const type& operator[](int i) const {return data[i];}

	// buffer I/O
	int buffer_size(int x0 = 0, int sx = 0) const;
	void to_buffer(char* buffer, int x0 = 0, int sx = 0) const;
	void from_buffer(char* buffer, int x0 = 0, int sx = 0);

	// file I/O
	void input(const char* filename);
	void input(const char* filename, int id, int np, int ng = 1);
	void output(const char* filename) const;
	void output(const char* filename, int id, int np, int ng = 1) const;

	// grid attribute functions
	int fields() const {return data[0].fields();}
	int size(int i) const {return nx[i];}
	bool boundary(int i) const {return px[i];}
	bool& boundary(int i) {return px[i];}
	float spacing(int i) const {return dx[i];}
	float& spacing(int i) {return dx[i];}

	// grid operations
	void swap(grid& temp) {data.swap(temp.data);}
  void ghostswap(int id, int np, int ng = 1);
  void set_boundary(int i, bool periodic){px[i] = periodic;};

protected:
	// grid data
	std::vector<type> data;
	int   nx[dimension];
	bool  px[dimension];
	float dx[dimension];
};

template <int dimension, typename type>
grid<dimension,type>::grid(int nf, ...)
{
	va_list list;
	va_start(list,nf);
	for (int i=0; i<dimension; i++) {
		nx[i] = va_arg(list,int);
		px[i] = false;
		dx[i] = 1.0;
	}
	va_end(list);

	if (dimension==1) data.resize(nx[0],type(nf));
	if (dimension==2) data.resize(nx[0],type(nf,nx[1]));
	if (dimension==3) data.resize(nx[0],type(nf,nx[1],nx[2]));
}

template <int dimension, typename type>
int grid<dimension,type>::buffer_size(int x0, int sx) const
{
	int buffer_size = 0;
	if (sx==0) sx = data.size();
	for (int x=x0; x<x0+sx; x++)
		buffer_size += data[x].buffer_size();
	return buffer_size;
}

template <int dimension, typename type>
void grid<dimension,type>::to_buffer(char* buffer, int x0, int sx) const
{
	char* p = buffer;
	if (sx==0) sx = data.size();
	for (int x=x0; x<x0+sx; x++) {
		data[x].to_buffer(p);
		p += data[x].buffer_size();
	}
}

template <int dimension, typename type>
void grid<dimension,type>::from_buffer(char* buffer, int x0, int sx)
{
	char* p = buffer;
	if (sx==0) sx = data.size();
	for (int x=x0; x<x0+sx; x++) {
		data[x].from_buffer(p);
		p += data[x].buffer_size();
	}
}

template <int dimension, typename type>
void grid<dimension,type>::input(const char* filename)
{
	// file open error check
	std::ifstream input(filename);
	if (!input) {
		std::cerr<<"file input error: "<<filename;
		std::cerr<<" cannot be opened."<<std::endl;
		exit(-1);
	}

	// dimension error check
	int dim;
	input.read(reinterpret_cast<char*>(&dim),sizeof(dim));
	if (dim!=dimension) {
		std::cerr<<"File input error: "<<filename;
		std::cerr<<" has grid of dimension "<<dim<<std::endl;
		exit(-2);
	}

	// read number of fields 
	int nf;
	input.read(reinterpret_cast<char*>(&nf),sizeof(nf));

	// read grid parameters
	for (int i=0; i<dimension; i++)
		input.read(reinterpret_cast<char*>(&nx[i]),sizeof(nx[i]));
	for (int i=0; i<dimension; i++)
		input.read(reinterpret_cast<char*>(&px[i]),sizeof(px[i]));
	for (int i=0; i<dimension; i++)
		input.read(reinterpret_cast<char*>(&dx[i]),sizeof(dx[i]));

	// resize grid data structure
	if (dimension==1) data.resize(nx[0],type(nf));
	if (dimension==2) data.resize(nx[0],type(nf,nx[1]));
	if (dimension==3) data.resize(nx[0],type(nf,nx[1],nx[2]));

	// read grid data
	int pos1 = input.tellg();
	input.seekg(0,std::ios::end);
	int pos2 = input.tellg();
	input.seekg(pos1,std::ios::beg);
	int size = pos2-pos1;
	char* buffer = new char[size];
	input.read(reinterpret_cast<char*>(buffer),size);
	this->from_buffer(buffer);
	delete [] buffer;

	// seed random number generators
	srand(time(NULL));
}

template <int dimension, typename type>
void grid<dimension,type>::output(const char* filename) const
{
	// file open error check
	std::ofstream output(filename);
	if (!output) {
		std::cerr<<"file input error: "<<filename;
		std::cerr<<" cannot be opened."<<std::endl;
		exit(-1);
	}

	// write grid dimension
	int dim = dimension;
	output.write(reinterpret_cast<const char*>(&dim),sizeof(dim));

	// write number of fields
	int nf = fields();
	output.write(reinterpret_cast<const char*>(&nf),sizeof(nf));

	// write grid parameters
	for (int i=0; i<dimension; i++)
		output.write(reinterpret_cast<const char*>(&nx[i]),sizeof(nx[i]));
	for (int i=0; i<dimension; i++)
		output.write(reinterpret_cast<const char*>(&px[i]),sizeof(px[i]));
	for (int i=0; i<dimension; i++)
		output.write(reinterpret_cast<const char*>(&dx[i]),sizeof(dx[i]));

	// write grid data
	int size = this->buffer_size();
	char* buffer = new char[size];
	this->to_buffer(buffer);
	output.write(reinterpret_cast<const char*>(buffer),size);
	delete [] buffer;
}


template < template<typename value_type> class data_type, typename value_type >
class grid1D : public grid<1,data_type<value_type> > {
public:
	// constructors
	grid1D(int x, int nf = 1)
		: grid<1,data_type<value_type> > (nf,x) {}
	grid1D(const char* filename)
		: grid<1,data_type<value_type> > (filename) {}
	grid1D(const char* filename, int id, int np, int ng = 1)
		: grid<1,data_type<value_type> > (filename,id,np,ng) {}

	// utility functions
	const data_type<value_type>& neighbor(int x, int sx) const;
};

template < template<typename value_type> class data_type, typename value_type >
const data_type<value_type>& grid1D<data_type,value_type>::
	neighbor(int x, int sx) const
{
	const grid1D& grid = *this;
	int  nx = grid.nx[0];
	bool px = grid.px[0];

	int xsx = x+sx;
	int ax = (xsx>=nx);
	int bx = (xsx<0);
	int i = xsx+(px)*(nx*(bx-ax))+(!px)*(ax*(nx-1)-(ax|bx)*xsx);

	return grid.data[i];
}


template < template<typename value_type> class data_type, typename value_type >
class grid2D : public grid<2,grid<1,data_type<value_type> > > {
public:
	// constructors
	grid2D(int x, int y, int nf = 1)
		: grid<2,grid<1,data_type<value_type> > > (nf,x,y) {}
	grid2D(const char* filename)
		: grid<2,grid<1,data_type<value_type> > > (filename) {}
	grid2D(const char* filename, int id, int np, int ng = 1)
		: grid<2,grid<1,data_type<value_type> > > (filename,id,np,ng) {}

	// utility functions
	const data_type<value_type>& neighbor(int x, int y, int sx, int sy) const;
};

template < template<typename value_type> class data_type, typename value_type >
const data_type<value_type>& grid2D<data_type,value_type>::
	neighbor(int x, int y, int sx, int sy) const
{
	const grid2D& grid = *this;
	int  nx = grid.nx[0];
	int  ny = grid.nx[1];
	bool px = grid.px[0];
	bool py = grid.px[1];

	int xsx = x+sx;
	int ax = (xsx>=nx);
	int bx = (xsx<0);
	int i = xsx+(px)*(nx*(bx-ax))+(!px)*(ax*(nx-1)-(ax|bx)*xsx);
	
	int ysy = y+sy;
	int ay = (ysy>=ny);
	int by = (ysy<0);
	int j = ysy+(py)*(ny*(by-ay))+(!py)*(ay*(ny-1)-(ay|by)*ysy);

	return grid.data[i][j];
}


template < template<typename value_type> class data_type, typename value_type >
class grid3D : public grid<3,grid<2,grid<1,data_type<value_type> > > > {
public:
	// constructors
	grid3D(int x, int y, int z, int nf = 1)
		: grid<3,grid<2,grid<1,data_type<value_type> > > > (nf,x,y,z) {}
	grid3D(const char* filename)
		: grid<3,grid<2,grid<1,data_type<value_type> > > > (filename) {}
	grid3D(const char* filename, int id, int np, int ng = 1)
		: grid<3,grid<2,grid<1,data_type<value_type> > > > (filename,id,np,ng) {}

	// utility functions
	const data_type<value_type>& neighbor(int x, int y, int z, int sx, int sy, int sz) const;
};

template < template<typename value_type> class data_type, typename value_type >
const data_type<value_type>& grid3D<data_type,value_type>::
	neighbor(int x, int y, int z, int sx, int sy, int sz) const
{
	const grid3D& grid = *this;
	int  nx = grid.nx[0];
	int  ny = grid.nx[1];
	int  nz = grid.nx[2];
	bool px = grid.px[0];
	bool py = grid.px[1];
	bool pz = grid.px[2];

	int xsx = x+sx;
	int ax = (xsx>=nx);
	int bx = (xsx<0);
	int i = xsx+(px)*(nx*(bx-ax))+(!px)*(ax*(nx-1)-(ax|bx)*xsx);
	
	int ysy = y+sy;
	int ay = (ysy>=ny);
	int by = (ysy<0);
	int j = ysy+(py)*(ny*(by-ay))+(!py)*(ay*(ny-1)-(ay|by)*ysy);

	int zsz = z+sz;
	int az = (zsz>=nz);
	int bz = (zsz<0);
	int k = zsz+(pz)*(nz*(bz-az))+(!pz)*(az*(nz-1)-(az|bz)*zsz);

	return grid.data[i][j][k];
}

} // namespace MMSP

#endif 
