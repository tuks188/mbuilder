// MC2vti.cpp
// Convert grids in MMSP MC format to VTK image data format
// Questions/comments to jgruber@andrew.cmu.edu (Jason Gruber)

#include"MCgrid.hpp"
#include<string>
#include<vector>

using namespace std;
using namespace MMSP;

void tokenize(const string& str,
                      std::vector<string>& tokens,
                      const string& delimiters = " ")
{
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);

    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));

        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);

        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}

int main(int argc, char* argv[])
{
  string outputname;
  if (argc==2) {
	// Change the name of the input filename for outout
	string filename = argv[1];
	std::vector<string> tokens;
	string delimeters = "."; // Only a period
	tokenize(filename, tokens, delimeters);
  
	outputname = tokens[0] + ".vti";
  }else if(argc==3){
	outputname = argv[2];
  }else if(argc<2) {
	std::cout<<"Usage: "<<argv[0]<<" file.MC file.vti\n";
	std::cout<<"\tor" <<  argv[0] << "file.MC" << endl;
	exit(-1);
  }
  
  std::ifstream input(argv[1]);
  if (!input) {
	std::cout<<argv[0]<<": File read error. Could not open file"<<argv[1]<<".\n";
	exit(-2);
  }

	int dim;
	input.read(reinterpret_cast<char*>(&dim),sizeof(dim));
	input.close();

	std::ofstream output(outputname.c_str());
	if (!output) {
		std::cout<<argv[0]<<": File write error. Could not open file"<< outputname <<".\n";
		exit(-2);
	}

	std::string byte_order;
	if (0x01 & static_cast<int>(1)) byte_order = "LittleEndian";
	else byte_order = "BigEndian";


	if (dim==1) {
		MCgrid1D grid(argv[1]);
		int   nx = grid.size(0);
		float dx = grid.spacing(0);

		output<<"<?xml version=\"1.0\"?>\n";
		output<<"<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\""<<byte_order<<"\">\n";
		output<<"  <ImageData WholeExtent=\""<<0<<" "<<nx<<" 0 0 0 0\"";
		output<<"   Origin=\"0 0 0\" Spacing=\""<<dx<<" 0 0\">\n";
		output<<"    <Piece Extent=\""<<0<<" "<<nx<<" "<<" 0 0 0 0\">\n";
		output<<"      <CellData Scalars=\"ID\">\n";
		output<<"        <DataArray Name=\"ID\" type=\"Int32\" format=\"ascii\">\n";
		for (int x=0; x<nx; x++)
			output<<grid[x]<<" ";
		output<<"\n";
		output<<"        </DataArray>\n";
		output<<"      </CellData>\n";
		output<<"    </Piece>\n";
		output<<"  </ImageData>\n";
		output<<"</VTKFile>\n";
	}

	if (dim==2) {
		MCgrid2D grid(argv[1]);
		int   nx = grid.size(0);
		int   ny = grid.size(1);
		float dx = grid.spacing(0);
		float dy = grid.spacing(1);

		output<<"<?xml version=\"1.0\"?>\n";
		output<<"<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\""<<byte_order<<"\">\n";
		output<<"  <ImageData WholeExtent=\""<<0<<" "<<nx<<" "<<0<<" "<<ny<<" 0 0\"";
		output<<"   Origin=\"0 0 0\" Spacing=\""<<dx<<" "<<dy<<" 0\">\n";
		output<<"    <Piece Extent=\""<<0<<" "<<nx<<" "<<0<<" "<<ny<<" 0 0\">\n";
		output<<"      <CellData Scalars=\"ID\">\n";
		output<<"        <DataArray Name=\"ID\" type=\"Int32\" format=\"ascii\">\n";
		for (int y=0; y<ny; y++)
			for (int x=0; x<nx; x++)
				output<<grid[x][y]<<" ";
		output<<"\n";
		output<<"        </DataArray>\n";
		output<<"      </CellData>\n";
		output<<"    </Piece>\n";
		output<<"  </ImageData>\n";
		output<<"</VTKFile>\n";
	}

	if (dim==3) {
		MCgrid3D grid(argv[1]);
		int   nx = grid.size(0);
		int   ny = grid.size(1);
		int   nz = grid.size(2);
		float dx = grid.spacing(0);
		float dy = grid.spacing(1);
		float dz = grid.spacing(2);

		output<<"<?xml version=\"1.0\"?>\n";
		output<<"<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\""<<byte_order<<"\">\n";
		output<<"  <ImageData WholeExtent=\""<<0<<" "<<nx<<" "<<0<<" "<<ny<<" "<<0<<" "<<nz<<"\"";
		output<<"   Origin=\"0 0 0\" Spacing=\""<<dx<<" "<<dy<<" "<<dz<<"\">\n";
		output<<"    <Piece Extent=\""<<0<<" "<<nx<<" "<<0<<" "<<ny<<" "<<0<<" "<<nz<<"\">\n";
		output<<"      <CellData Scalars=\"ID\">\n";
		output<<"        <DataArray Name=\"ID\" type=\"Int32\" format=\"ascii\">\n";
		for (int z=0; z<nz; z++) 
			for (int y=0; y<ny; y++)
				for (int x=0; x<nx; x++)
					output<<grid[x][y][z]<<" ";
		output<<"\n";
		output<<"        </DataArray>\n";
		output<<"      </CellData>\n";
		output<<"    </Piece>\n";
		output<<"  </ImageData>\n";
		output<<"</VTKFile>\n";
	}

	output.close();
}

