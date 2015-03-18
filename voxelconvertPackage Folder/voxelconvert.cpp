#include"voxelio.hpp"
#include<iostream>
#include<string.h>

int main(int argc, char* argv[])
{
  //Check to see if there are the proper commandline arguments
  if(argc == 1 || argc < 3){
    cout << "USAGE: ./voxelconvert infile.[mc,ph,dx] outfile.[mc,ph,dx]" << endl;
    cout << endl;
    cout << "This program will convert mc grid data between mc, ph, and dx file formats." << endl;
    exit(1);
  }
  string infile = argv[1];
  string outfile = argv[2];
  string inEXT, outEXT;
  int nx, ny, nz;
  std::vector<int> data;
  
  // Process the in file name and determine what type of input we have.
  std::vector<string> tokens;
  string delimeters = "."; // Only a period
  tokenize(infile, tokens, delimeters);
  
  if (tokens.size() > 2){
    cout << "WARNING: multiple \".\" in file name: " << infile << endl;
    cout << "\t Using: ." << tokens.back() << "As extension" << endl;
    inEXT =  tokens.back();
  }else if( tokens.size() == 2 ){
    inEXT = tokens.back();
  }
  
  // Process the out file name and determine what type of output we want.
  tokens.clear();
  tokenize(outfile, tokens, delimeters);
  
  if (tokens.size() > 2){
    cout << "WARNING: multiple \".\" in file name: " << outfile << endl;
    cout << "\t Using: ." << tokens.back() << "As extension" << endl;
    outEXT =  tokens.back();
  }else if( tokens.size() == 2 ){
    outEXT = tokens.back();
  }

  // If the infile is *.mc
  if(inEXT == "mc" || inEXT == "MC"){
    if( outEXT == "mc" || outEXT =="MC"){
      cout << "Infile and outfile are of the same type. Nothing to be done." << endl;
    }else if(outEXT == "ph" || outEXT == "PH"){
      cout << "Read in the initial MC structure" << endl;
      MCgrid3D grid(infile.c_str());
      WritePHFile(outfile, grid);
    }else if(outEXT == "dx" || outEXT == "DX"){
      MCgrid3D grid(infile.c_str());
      WriteDXFile(outfile, grid);
    }

    // If the infile is *.ph
  } else if(inEXT == "ph" || inEXT == "PH"){
    if(outEXT == "ph" || outEXT == "PH"){
      cout << "Infile and outfile are of the same type. Nothing to be done." << endl;
    }else if( outEXT == "mc" || outEXT =="MC"){ 
      cout << "Read in the initial PH structure" << endl;
      ReadPHFile(infile, data, nx, ny, nz) ;
      cout << "Initialize grid with data from PH file" << endl;
      MCgrid3D grid(nx,ny,nz);
      for(int z=0; z < nz; z++)
        for(int y=0; y < ny; y++)
          for(int x=0; x < nx; x++)
            grid[x][y][z] = data[x+(nx*y)+(ny*nx*z)];

      grid.output(argv[2]);
    }else if(outEXT == "dx" || outEXT == "DX"){
      cout << "Read in the initial DX structure" << endl;
      ReadPHFile(infile, data, nx, ny, nz);
      // Initialize grid with data from PH file
      MCgrid3D grid(nx,ny,nz);
      for(int z=0; z < nz; z++)
        for(int y=0; y < ny; y++)
          for(int x=0; x < nx; x++)
            grid[x][y][z] = data[x+(nx*y)+(ny*nx*z)];
          
      WriteDXFile(outfile, grid);
    }
        
  } else if(inEXT == "dx" || inEXT == "DX"){
    if(outEXT == "dx" || outEXT == "DX"){
      cout << "Infile and outfile are of the same type. Nothing to be done." << endl;
    }else if( outEXT == "mc" || outEXT =="MC"){ 
      // Read in the initial structure
      ReadDXFile(infile, data, nx, ny, nz);
      // Initialize grid with data from DX file
      MCgrid3D grid(nx,ny,nz);
      for(int z=0; z < nz; z++)
        for(int y=0; y < ny; y++)
          for(int x=0; x < nx; x++)
            grid[x][y][z] = data[x+(nx*y)+(ny*nx*z)];
      
      grid.output(argv[2]);
    }else if(outEXT == "ph" || outEXT == "PH"){
      // Read in the initial structure
      ReadDXFile(infile, data, nx, ny, nz);
      // Initialize grid with data from DX file
      MCgrid3D grid(nx,ny,nz);
      for(int z=0; z < nz; z++)
        for(int y=0; y < ny; y++)
          for(int x=0; x < nx; x++)
            grid[x][y][z] = data[x+(nx*y)+(ny*nx*z)];
      
      WritePHFile(outfile, grid);
    }
  }
}
