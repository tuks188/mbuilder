// MCgrid.cpp
// Minimal single processor Monte Carlo code using MMSP
// Questions/comments to jgruber@andrew.cmu.edu (Jason Gruber)

#include"MCgrid.update.hpp"
#include<string>
#include<fstream>
#include<vector>
#include<iostream>
#include<iomanip>
using namespace MMSP;
using namespace std;


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

void  ReadPHFile(string FileName, std::vector<int> &data, int &nx, int &ny, int &nz)
{
  string line; 
  string delimeters(", ;\t"); /* delimeters to split the data */
  std::vector<string> tokens; /* vector to store the split data */
  //std::vector<int> data; /* vector to store the data */

  int error, spin; /* dummy variables */
  //int nx, ny, nz;
  
  ifstream InFile;
  InFile.open(FileName.c_str());
  if(!InFile)
	{
      cout << "Failed to open: " << FileName << endl;
      exit(0);
    }

  getline(InFile, line, '\n');
  tokenize(line,tokens,delimeters);
            
  // Process the header information from the PH file.
  error = 0;
  error += sscanf(tokens[0].c_str(), "%d", &nx);
  error += sscanf(tokens[1].c_str(), "%d", &ny);
  error += sscanf(tokens[2].c_str(), "%d", &nz);
  tokens.clear();

  //  cout << "INFO: PH file grid size: " << nx << "\t" << ny << "\t" << nz << endl;;

  //MCgrid3D* grid = new grid(nx,ny,nz);

  // Get the remaining two lines of the header and ignore
  getline(InFile, line, '\n');
  getline(InFile, line, '\n');

  //The PH file has a unique format of 20 entries on each line. I have
  //now idea who initiated this insanity but I am about to propetuate
  //it.
  //
  //The most simple thing todo is to read the entire dataset into one
  //long vector and then read that vector to assign values to the grid
  
  while(getline(InFile, line, '\n') != NULL)
	{
	  tokens.clear();
	  error = 0;
	  tokenize(line,tokens,delimeters);
// 	  cout << line << endl;
// 	  for(int i=0; i < tokens.size(); i++ )
// 		cout << setw(6) << tokens[i];
// 	  cout << endl;
	  
	  for(int in_spins=0; in_spins < tokens.size(); in_spins++)
		{
		  error += sscanf(tokens[in_spins].c_str(), "%d", &spin);
		  data.push_back(spin);
		}
// 	  if(error != 20)
// 		{
// 		  cout << "ERROR: Invalid number of line entries in PH file" << endl;
// 		}
	}
		  
  tokens.clear();
 
  InFile.close();
  //  return();
}

void WritePHFile(std::string filename, MCgrid3D &grid)
{
  string OutputName;
  
  // Change the name of the input filename for outout
  std::vector<string> tokens;
  string delimeters = "."; // Only a period
  tokenize(filename, tokens, delimeters);
  
  OutputName = tokens[0] +".ph";
  ofstream outfile;
  outfile.open(OutputName.c_str());
  if(!outfile)
	{
	  cout << "Failed to open: " << OutputName << endl;
	  exit(0);
	}

  // Find the unique number of grains
  std::map<int,bool> used;

    for (int x=0; x<grid.size(0); x++)
        for (int y=0; y<grid.size(1); y++)
            for (int z=0; z<grid.size(2); z++)
                used[grid[x][y][z]] = true;

    int grains = 0;
    typedef std::map<int,bool>::iterator iterator;
    for (iterator i = used.begin(); i!=used.end(); i++)
        if ((*i).second == true) grains++;

    //std::cout<<grains<< " " << used.size() << std::endl;


	outfile << "     " << grid.size(0)
			<< "     " << grid.size(1)
			<< "     " << grid.size(2) << endl ;
	outfile << "\'grwXXX\'              52.00  1.000  1.0       " << grains  << endl;
	outfile << " 3.000 0.000 0.000          0        \n"; // << grains << endl;
	
	int count =0;
	for(int k=0; k < grid.size(2); k++)
	  for(int j=0; j < grid.size(1); j++)
		for(int i=0; i < grid.size(0); i++)
		  {
			outfile << setw(6)<< grid[i][j][k] ;
			count++;
		  if(count==20 || (i == (grid.size(0)-1)))
			{
			  outfile << endl;
			  count =0;
			}
		  //			outfile << grid[i][j][k] << endl;
		  }
	outfile << endl;
	outfile.close();
}

void WriteDXFile(std::string filename, MCgrid3D &grid)
{
  string OutputName;
  
  // Change the name of the input filename for outout
  std::vector<string> tokens;
  string delimeters = "."; // Only a period
  tokenize(filename, tokens, delimeters);
  
  OutputName = tokens[0]+".dx";
  ofstream outfile;
  outfile.open(OutputName.c_str());
  if(!outfile)
	{
	  cout << "Failed to open: " << OutputName << endl;
	  exit(0);
	}

   // Find the unique number of grains
  std::map<int,bool> used;

    for (int x=0; x<grid.size(0); x++)
        for (int y=0; y<grid.size(1); y++)
            for (int z=0; z<grid.size(2); z++)
                used[grid[x][y][z]] = true;

    int grains = 0;
    typedef std::map<int,bool>::iterator iterator;
    for (iterator i = used.begin(); i!=used.end(); i++)
        if ((*i).second == true) grains++;


	outfile << "object 1 class gridpositions counts "
			<< grid.size(2)+1 << " "
			<< grid.size(1)+1 << " "
			<< grid.size(0)+1 << endl;	
	outfile << "origin   0    0    0" << endl;
	outfile << "delta    0    0    1" << endl;
	outfile << "delta    0    1    0" << endl;
	outfile << "delta    1    0    0" << endl;
	outfile << endl;
	outfile << "object 2 class gridconnections counts  "
			<< grid.size(2)+1 << " "
			<< grid.size(1)+1 << " "
			<< grid.size(0)+1 << endl;	
	outfile << endl;
	outfile << "object 3 class array type int rank 0 items "<< grid.size(2)*grid.size(1)*grid.size(0)
			<<" data follows" << endl;
	
	int count =0;
	for(int k=0; k < grid.size(2); k++)
	  for(int j=0; j < grid.size(1); j++)
		for(int i=0; i < grid.size(0); i++)
		  {
			outfile << setw(6)<< grid[i][j][k] ;
			count++;
			if(count==20 || (i == (grid.size(0)-1)))
			  {
				outfile << endl;
				count =0;
			  }
		  }

	outfile <<	"attribute \"dep\" string \"connections\"" << endl;
	outfile << endl;
	outfile << "# Create a field." << endl;
	outfile << "object \" " << grains << "  grain microstructure\" class field" <<endl;
	outfile << "component \"positions\" value 1" << endl;
	outfile << "component \"connections\" value 2" << endl;
	outfile << "component \"data\" value 3" << endl;
	outfile << endl;
	outfile << "end" << endl;

	outfile << endl;
	outfile.close();
}

void WriteImplicitFile(std::string filename, MCgrid3D &grid)
{
  string OutputName;
  
  // Change the name of the arguement filename for outout
  std::vector<string> tokens;
  string delimeters = "."; // Only a period
  tokenize(filename, tokens, delimeters);
  
  OutputName = tokens[0] + ".zyx";
  ofstream outfile;
  outfile.open(OutputName.c_str());
  if(!outfile)
	{
	  cout << "Failed to open: " << OutputName << endl;
	  exit(0);
	}

  for(int x=0; x < grid.size(0); x++)
	for(int y=0; y < grid.size(1); y++)
	  for(int z=0; z < grid.size(2); z++)
		outfile << grid[x][y][z] << endl ;

  outfile << endl;
  outfile.close();
}

int main(int argc, char* argv[])
{
  // Local variables
  std::string infile = argv[1];
  std::string outfile = argv[2];
  std::vector<int> data;
  int nx, ny, nz, error;
  bool periodic = true;
  int min_gs;
  int x1, x2, y1, y2, z1 ,z2;
  int nxc, nyc, nzc; // new dimensions for cropped structure

  x1 = x2 = y1 = y2 = z1 = z2 = 0;

  cout << "Enter the crop corners [x1, x2, y1, y2, z1, z2]" << endl;
  cin >> x1 >> x2 >> y1 >> y2 >> z1 >> z2;

  nxc = x2 - x1;
  nyc = y2 - y1;
  nzc = z2 - z1;
  

  // Process the file name and determine what type of input we have.
  std::vector<string> tokens;
  string delimeters = "."; // Only a period
  string extension;
  tokenize(infile, tokens, delimeters);
  
  if (tokens.size() > 2){
	cout << "WARNING: multiple \".\" in file name: " << infile << endl;
	cout << "\t Using: ." << tokens.back() << "As extension" << endl;
	  extension =  tokens.back();
  }else if( tokens.size() == 2 ){
	extension = tokens.back();
  }

  if(extension == "mc" || extension == "MC"){
	MCgrid3D grid(infile.c_str());
	MCgrid3D cgrid(nxc, nyc, nzc);
	for(int z=0; z < nz; z++)
	  for(int y=0; y < ny; y++)
		for(int x=0; x < nx; x++)
		  if(z >= z1 && z < z2 && x >= x1 && x < x2 && y >= y1 && y < y2)
			cgrid[x-x1][y-y1][z-z1] = grid[x][y][z];

	// Output the final structure as binary
	cgrid.output(argv[2]);

	// Output the final structure as ascii in several formats
	WritePHFile(outfile, cgrid);
	WriteImplicitFile(outfile, cgrid);
	WriteDXFile(outfile, cgrid);

	
  } else if(extension == "ph" || extension == "PH"){
	// Read in the initial structure
	ReadPHFile(infile, data, nx, ny, nz);

	// Initialze grid with data from PH file
	MCgrid3D grid(nx,ny,nz);
	for(int z=0; z < nz; z++)
	  for(int y=0; y < ny; y++)
		for(int x=0; x < nx; x++)
		  grid[x][y][z] = data[x+(nx*y)+(ny*nx*z)];

	MCgrid3D cgrid(nxc, nyc, nzc);
	for(int z=0; z < nz; z++)
	  for(int y=0; y < ny; y++)
		for(int x=0; x < nx; x++)
		  if(z >= z1 && z < z2 && x >= x1 && x < x2 && y >= y1 && y < y2)
			cgrid[x-x1][y-y1][z-z1] = grid[x][y][z];	// Set the periodic nature of the boundaries.

	// Output the final structure as binary
	cgrid.output(argv[2]);
  
	// Output the final structure as ascii in several formats
	WritePHFile(outfile, cgrid);
	WriteImplicitFile(outfile, cgrid);
	WriteDXFile(outfile, cgrid);
  }
}
