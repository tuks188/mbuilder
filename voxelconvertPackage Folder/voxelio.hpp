//voxelip.hpp Header file for reading and writing various voxel file
//formats. 
// Questions/comments to stephen.sintay@gmail.com (Stephen Sintay)

//  updated DX read in, ADR Sep 08

#ifndef VOXELIO_HPP
#define VOXELIO_HPP

#include"MCgrid.hpp"
#include"tokenize.hpp"
#include<iomanip>
#include<string.h>

using namespace MMSP;
using namespace std;



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
//        cout << line << endl;
//        for(int i=0; i < tokens.size(); i++ )
//              cout << setw(6) << tokens[i];
//        cout << endl;
          
          for(int in_spins=0; in_spins < tokens.size(); in_spins++)
                {
                  error += sscanf(tokens[in_spins].c_str(), "%d", &spin);
                  data.push_back(spin);
                }
//        if(error != 20)
//              {
//                cout << "ERROR: Invalid number of line entries in PH file" << endl;
//              }
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
                  //                    outfile << grid[i][j][k] << endl;
                  }
        outfile << endl;
        outfile.close();
}


void  ReadDXFile(string FileName, std::vector<int> &data, int &nx, int &ny, int &nz)
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
            
  // Process the header information and look for the string "counts"
  // Then read the data size after that
  int pos1 = 0 ;
  while(pos1 == 0 ) { // continue until we find the keyword
    for(int i=0; i< tokens.size(); i++){
      if( tokens[i] == "counts" ){
	pos1 = i;
      }
    }
    // Read the next line of the header if we did not find the keyword
    // in the line
    if(pos1 == 0){ 
      tokens.clear();
      getline(InFile, line, '\n');
      tokenize(line,tokens,delimeters);
      if (tokens.size() == 20){
	cout << "ERROR: Unable to read data dimensions from the header" << endl;
	InFile.close();
	exit(1);
      }
    }
        
  }

  //  edited ADR to deal with the header line variations in DX files
  int headerType = 0 ;
  cout << "Is this a Z-Y-X header (normal) or a X-Y-Z header (Sukbin stack)? " << endl ;
  cout << "Answer 0 for Z-Y-X header, or 1 for X-Y-Z header:  " ;
  cin >> headerType ;
  if ( headerType < 0 || headerType > 1 ) {
    cout << "Has to be 0 or 1!" << endl ;
    exit(1) ;
  }
  if(pos1 != 0 ){
        error = 0;
	if ( headerType == 0 ) error += sscanf(tokens[pos1+1].c_str(), "%d", &nz) ;
	if ( headerType == 1 ) error += sscanf(tokens[pos1+1].c_str(), "%d", &nx) ;
        error += sscanf(tokens[pos1+2].c_str(), "%d", &ny);
	if ( headerType == 0 ) error += sscanf(tokens[pos1+3].c_str(), "%d", &nx) ;
	if ( headerType == 1 ) error += sscanf(tokens[pos1+3].c_str(), "%d", &nz) ;
        tokens.clear();
        // The dimensions listed in the DX file are always one greater
        // than the actual dimensions
        nx--;
        ny--;
        nz--;
  }

   cout << "INFO: DX data dimensions: " << endl;
   cout << "nz= " << nz << endl;
   cout << "ny= " << ny << endl;
   cout << "nx= " << nx << endl;;

  //The DX file has a unique format of 20 entries on each line. I have
  //no idea who initiated this insanity but I am about to perpetuate
  //it.
  //
  //The most simple thing to do is to read the entire dataset into one
  //long vector and then read that vector to assign values to the grid

  //  ADR:  6 Sep 08; time to make the input much more general!
  //  equivalent to list-direcvted input in Fortran, actually !!

  pos1 = 0 ;
  while(pos1 == 0 ) { // continue until we find the keyword
    for(int i=0; i< tokens.size(); i++){
      if( tokens[i] == "items" ){
	pos1 = i;
      }
    }
    // Read the next line of the header if we did not find the keyword
    // in the line
    if(pos1 == 0){ 
      tokens.clear();
      getline(InFile, line, '\n');
      tokenize(line,tokens,delimeters);
      if (tokens.size() == 20){
	cout << "ERROR: Unable to locate the last header line" << endl;
	InFile.close();
	exit(1);
      }
    }
  }  // when we get here, we are looking at data
  int points = 0 ;
  if(pos1 != 0 ){
        error = 0;
        error += sscanf(tokens[pos1+1].c_str(), "%d", &points);
        tokens.clear();
  }
  cout << "Compare no. points " << points 
       << " with x*y*z: " << nx*ny*nz << endl ;

  bool finished_header, finished_data;
  finished_header = true;
  //  finished_header = false;
  finished_data = false;
  while(getline(InFile, line, '\n') != NULL){
        
    // Get the remaining lines of the header and ignore
    tokens.clear();
    error = 0;
    tokenize(line,tokens,delimeters);
    
    //    if(tokens.size()==20)
    //      finished_header = true;

    if(finished_header && ((tokens[0] == "attribute") || 
			   data.size() == nz*ny*nx))
      finished_data = true ;
        
    if(finished_header && !finished_data){
      for(int in_spins = 0; in_spins < tokens.size() ; in_spins++ ) {
	error += sscanf( tokens[in_spins].c_str() , "%d", &spin ) ;
	data.push_back(spin);
      }
    }
  }

  if(data.size() != (nz*ny*nx)){
        cout << "ERROR: data size does not match header dimensions" << endl;
        cout << "\t" << data.size() << "\t" << nz*nx*ny << endl;
        exit(1);
        InFile.close();
  }
  
  tokens.clear();
  InFile.close();
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

        outfile <<      "attribute \"dep\" string \"connections\"" << endl;
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

void  ReadZYXFile(string FileName, std::vector<int> &data, int &nx, int &ny, int &nz)
{
  // The data is store with one entry per line and the first line is
  // the dimension
  
  string line; 
  string delimeters(", ;\t"); /* delimeters to split the data */
  std::vector<string> tokens; /* vector to store the split data */
  
  int error, spin; /* dummy variables */
    
  ifstream InFile;
  InFile.open(FileName.c_str());
  if(!InFile)
        {
      cout << "Failed to open: " << FileName << endl;
      exit(0);
    }

  getline(InFile, line, '\n');
  tokenize(line,tokens,delimeters);
            
  // Process the header information for the data size "nx ny nz"
  error += sscanf(tokens[0].c_str(), "%d", &nx);
  error += sscanf(tokens[1].c_str(), "%d", &ny);
  error += sscanf(tokens[2].c_str(), "%d", &nz);
  

  //The DX file has a unique format of 20 entries on each line. I have
  //no idea who initiated this insanity but I am about to propetuate
  //it.
  //
  //The most simple thing todo is to read the entire dataset into one
  //long vector and then read that vector to assign values to the grid
  bool finished_header, finished_data;
  finished_header = false;
  finished_data = false;
  while(getline(InFile, line, '\n') != NULL){
        
        // Get the remaining lines of the header and ignore
        tokens.clear();
        error = 0;
        tokenize(line,tokens,delimeters);

        if(tokens.size()==20)
          finished_header = true;

        if(finished_header && ((tokens[0] == "attribute") || data.size() == nz*ny*nx))
          finished_data = true;
        
        if(finished_header && !finished_data){
          for(int in_spins=0; in_spins < tokens.size(); in_spins++){
                error += sscanf(tokens[in_spins].c_str(), "%d", &spin);
                data.push_back(spin);
          }
        }
  }

  if(data.size() != (nz*ny*nx)){
        cout << "ERROR: data size does not match header dimensions" << endl;
        cout << "\t" << data.size() << "\t" << nz*nx*ny << endl;
        exit(1);
        InFile.close();
  }
  
  tokens.clear();
  InFile.close();
}




void WriteZYXFile(std::string filename, MCgrid3D &grid)
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

#endif
