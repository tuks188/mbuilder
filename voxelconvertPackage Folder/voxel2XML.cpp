/*
*****************************************************************************************

Randomize_unique_micro_input.c

This program give randomly a unique spin number for each grain in micro.input.

Sukbin Lee
sukbin@andrew.cmu.edu

Wean Hall, 4340
Materials Science and Engineering Department
Carnegie Mellon University
Pittsburgh, PA 15213

*****************************************************************************************
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <map>
#include <fstream>
#include <iostream>
#include <string.h>
#include "voxelio.hpp"

#define PI 3.14159265
#define num_neigh 26

// constants for procedure ran3
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

// structure voxel stores the spin value and coordinates...
// Each voxel id is the array index of the "point"...
// I'll use the array index from 1 and 0 will be garbage...
struct voxel{
  int spin;
  int newspin;
};

struct volume{
  int spin;
  int size;
  int *member;
};

// Structure neighbor stores a list of neighbor site_ids for each voxel.
// Since we are dealing with 3d, we have 26 neighbors and the array index starts
// from 1...(It's easier than 0 when we deal with periodic boundary conditions...)
struct neighbor{
  int neigh_id[num_neigh+1];
};


// function prototypes...
void get_neighbor_list (struct neighbor *n, int ns, int nsp, int xDim, int yDim, int zDim);
void initialize (struct voxel *p, int ns);
int get_number_grains (struct voxel *p, struct neighbor *n, int ns, int nsp, int xDim, int yDim, int zDim);
void get_grains (struct voxel *p, struct neighbor *n, struct volume *g, 
                 int ns, int nsp, int xDim, int yDim, int zDim);
void get_final_grains (int size_limit, struct voxel *p, struct neighbor *n, struct volume *g, 
                 int ns, int nsp, int xDim, int yDim, int zDim);
void get_output_xml (std::string filname, struct voxel *p, struct neighbor *n, struct volume *g, 
                 int ns, int nsp, int xDim, int yDim, int zDim);
void randomize_grain_spin(struct voxel *p, struct volume *g, int ng);
void get_output_dx (struct voxel *p, int ns, int xDim, int yDim, int zDim);
double ran3 (long *idum);
int q_ran (int q);

// Main function starts...
int main(int argc, char* argv[]){

    //Check to see if there are the proper commandline arguments
  if(argc == 1 || argc < 2){
        cout << "USAGE: ./voxel2XML infile.[mc,ph,dx] " << endl;
        cout << endl;
        cout << "This program will extract an XML file from .mc, .ph, or .dx "<<endl;
        cout << "file formats for use with texturelist." << endl;
        exit(1);
  }
  int NS; // The number of sites(voxels) in the simulation box... 
  int NSP;
  int xnum, ynum, znum;
  int numG;
  int i, j, k;
  int num_input;
  int num_bin;
  int bin_size;
  int size_limit;
  char dummy[255];
  std::vector<int> data;

  std::string infile = argv[1];
  std::string outfile;
  std::string inEXT;
  
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

  outfile = tokens[0] + ".xml";

  // Load the data sukbin data objects
  // If the infile is *.mc
  if(inEXT == "mc" || inEXT == "MC"){
        MCgrid3D grid(infile.c_str());
        xnum = grid.size(0);
        ynum = grid.size(1);
        znum = grid.size(2);

        //copy the grid[][][] into data[]
        for(int z=0; z<znum; z++)
          for(int y=0; y<ynum; y++)
                for(int x=0; x<xnum; x++)
                  data.push_back(grid[x][y][z]);
        
        // If the infile is *.ph
  } else if(inEXT == "ph" || inEXT == "PH"){
        ReadPHFile(infile, data, xnum, ynum, znum);
        
        // If the infile is *.dx
  } else if(inEXT == "dx" || inEXT == "DX"){
        ReadDXFile(infile, data, xnum, ynum, znum);
  } else{
        cout <<"Unknown extension: "<< inEXT << endl;
        exit(1);
  }
  
  printf("\nPlease enter the low limit of size of grains (cutoff size):\n");
  scanf("%d", &size_limit);
  
  NS  = xnum*ynum*znum;
  NSP = xnum*ynum;

  struct neighbor *neigh;
  struct voxel *point;
  struct volume *grain;

  neigh = (struct neighbor *)malloc((NS+1)*sizeof(struct neighbor));
  point = (struct voxel *)malloc((NS+1)*sizeof(struct voxel));

  point[0].spin =  0;

  for(int i=0; i <data.size();i++){ 
        point[i+1].spin = data[i];
        point[i+1].newspin = point[i+1].spin;
  }

   data.clear();
   //   printf("\nReading microstructure...\n");
   //initialize(point, NS);

   //printf("\nFinding neighbors for each site...\n");
  get_neighbor_list(neigh, NS, NSP, xnum, ynum, znum);

  //printf("\nGetting the grains...\n");
  get_final_grains(size_limit, point, neigh, grain, NS, NSP, xnum, ynum, znum);

  //  numG = get_number_grains (point, neigh, NS, NSP, xnum, ynum, znum);
  //grain = (struct volume *)malloc(numG*sizeof(struct volume));
  
  //get_grains(point, neigh, grain, NS, NSP, xnum, ynum, znum);

  //  printf("\nAssining unique, random id number for each grain...\n");
  // randomize_grain_spin(point, grain, numG);
  
  //printf("\nOutput dx files...\n");
  //get_output_dx (point, NS, xnum, ynum, znum);

  get_output_xml(outfile, point, neigh, grain, NS, NSP, xnum, ynum, znum);
  
  printf("Done!\n");

  free(neigh);
  free(point);
  // free(grain);

  return 0;
}

void get_output_xml (std::string filename, struct voxel *p, struct neighbor *n, struct volume *g, 
                 int ns, int nsp, int xDim, int yDim, int zDim)
{
  std::ofstream outfile;
  outfile.open(filename.c_str());
  if(!outfile){
        std::cout << "Failed to open: " << filename << std::endl;
        exit(0);
  }
  
  std::map<int, int> patch_list; // Map from patch ID to patch size
  int csite, numG;
  bool finished;
  finished = false;

  printf("\nCounting the number of grains in the input microstructure...\n");
  numG = get_number_grains (p, n, ns, nsp, xDim, yDim, zDim);

  printf("\tNumber of grains: %6d\n", numG);

  g = (struct volume *)malloc(numG*sizeof(struct volume));
  get_grains(p, n, g, ns, nsp, xDim, yDim, zDim);

  outfile << "<nregions>" << numG << "</nregions>" << std::endl;
  for(int i=0; i < numG; i++){
        outfile<<"<region>" << std::endl;
        outfile<<"\t<id> " << g[i].spin <<" </id>"<<std::endl;
        outfile<<"\t<volume> "<<g[i].size<<" </volume>"<<std::endl;
        outfile<<"\t<approximating_ellipsoid>"<<std::endl;
        outfile<<"\t\t<ecenter> 0 0 0 </ecenter>"<<std::endl;
        outfile<<"\t\t<semi_axis_lengths> 0 0 0 </semi_axis_lenth>"<<std::endl;
        outfile<<"\t\t<axis1> 1.0  0.0 0.0 </axis1>"<<std::endl;
        outfile<<"\t\t<axis2> 0.0  1.0 0.0 </axis2>"<<std::endl;
        outfile<<"\t\t<axis3> 0.0  0.0 1.0 </axis3>"<<std::endl;
        outfile<<"\t\t<orientation> 0 0 0 \t</orientation>"<<std::endl;
        outfile<<"\t<approximating_ellipsoid>"<<std::endl;
                                                                                                                
        // Populated the patch list with patch sizes for grain[i]
        for(int j=0; j < g[i].size; j++){
          csite = g[i].member[j];
          for(int k=1; k <=26; k++){
                // Test to see if the spin of the neighbor is different then
                // the spin of the grain
                if(p[n[csite].neigh_id[k]].spin != g[i].spin){
                  patch_list[p[n[csite].neigh_id[k]].spin]++;
                  }
                }
          }
          
        outfile<<"\t<npatches> "<<patch_list.size()<<" </npatches>"<<std::endl;
        
        typedef std::map<int,int>::iterator iterator;
        for (iterator iter = patch_list.begin(); iter!=patch_list.end(); iter++){
          outfile<<"\t\t<patch>"<<std::endl;
          outfile<<"\t\t\t<neighbor> "<<(*iter).first<<" </neighbor>"<<std::endl;
          outfile<<"\t\t\t<average_outward_normal> 1.0 1.0 1.0 </average_outward_normal>"<<std::endl;
          outfile<<"\t\t\t<area> "<<(*iter).second<<" </area>"<<std::endl;
          outfile<<"\t\t\t<nfacets> "<<(*iter).second<<" </nfacets>"<<std::endl;
          outfile<<"\t\t\t</patch>"<<std::endl;
        }
        outfile<<"</region>"<<std::endl;
          
        patch_list.clear();
  }
  free(g);
}

void get_final_grains (int size_limit, struct voxel *p, struct neighbor *n, struct volume *g, 
                 int ns, int nsp, int xDim, int yDim, int zDim)
{
  std::vector<int> remove_list;
  std::map<int, int> patch_list; // Map from patch ID to patch size
  int new_spin=0;
  int largest_patch = 0;
  int csite, numG;
  bool finished;
  finished = false;

  while(!finished){

        printf("\nCounting the number of grains in the input microstructure...\n");
        numG = get_number_grains (p, n, ns, nsp, xDim, yDim, zDim);

        printf("\tNumber of grains: %6d\n", numG);

        g = (struct volume *)malloc(numG*sizeof(struct volume));
        get_grains(p, n, g, ns, nsp, xDim, yDim, zDim);
        
        //count the number of grains that need to be consumed
        for(int i=0; i < numG; i++)
          if(g[i].size < size_limit)
                remove_list.push_back(i);

        // If the list is empty then continue
        if(remove_list.size() == 0){
          finished = true;
          printf("\nGrain size thresholding complete\n");

          randomize_grain_spin(p, g, numG);
          free(g);

          continue;
        }
          
        printf("\nNumber of small grains: %6d\n",remove_list.size());

        // Absorb the small grains into neighbor with most continuity
        // Loop through the list of small grains
        for(int i=0; i < remove_list.size(); i++){
          
          // Populated the patch list with patch sizes for grain[i]
          for(int j=0; j < g[remove_list[i]].size; j++){
                csite = g[remove_list[i]].member[j];
                for(int k=1; k <=26; k++){
                  // Test to see if the spin of the neighbor is different then
                  // the spin of the grain
                  if(p[n[csite].neigh_id[k]].spin != g[remove_list[i]].spin){
                        patch_list[p[n[csite].neigh_id[k]].spin]++;
                  }
                }
          }
          
          // Find the neighbor with: [the most continuity && a volume larger
          // than the small grain we are testing] This way smaller grains
          // will be consumed only by larger grains
          new_spin= g[remove_list[i]].spin;
          largest_patch = 0;
          typedef std::map<int,int>::iterator iterator;
          for (iterator iter = patch_list.begin(); iter!=patch_list.end(); iter++){
                //              if( ((*iter).second > largest_patch) && (g[(*iter).first].size > g[remove_list[i]].size)){
                if( (*iter).second > largest_patch ){
                  largest_patch = (*iter).second; 
                  new_spin = (*iter).first;
                } 
          }
          
          // Change the id of the sites associated with this grain to the
          // new spin value. 
          g[remove_list[i]].spin = new_spin;
          for(int j=0; j < g[remove_list[i]].size; j++){
                p[g[remove_list[i]].member[j]].spin = new_spin;
                p[g[remove_list[i]].member[j]].newspin = new_spin;
          }
          
          patch_list.clear();
        }
        
        randomize_grain_spin(p, g, numG);
        // Recover the spin...
        for(int i=1; i<=ns; i++){
          p[i].newspin = p[i].spin;
        }
                
        // Free the memory allocated and start again
        free(g);
        remove_list.clear();
  }
}
  

void get_neighbor_list(struct neighbor *n, int ns, int nsp, int xDim, int yDim, int zDim){

  // nsp = number of sites in a plane of xDim by yDim...
  // n[][] = 2 dimensional array storing its site number and neighbors...
  // site_id = id number for each site...starting from 1 to xDim*yDim*zDim....
  //
  // I assumed the square lattice...so the order of neighbors as follows...
  //
  //    4   3   2         13  12  11          22  21  20
  //    5 site  1         14   9  10          23  18  19 
  //    6   7   8         15  16  17          24  25  26
  //    
  //    in-plane          upper plane         lower plane

 int i, j, k;   // loop indices...
 int site_id;   // id number for each site...

  for (k=0; k<=(ns-nsp); k=k+nsp){
    for (j=0; j<=(nsp-xDim); j=j+xDim){
      for (i=1; i<=xDim; i++){
      
        site_id = k + j + i;

        //same plane...
        n[site_id].neigh_id[1]  = k + j + i%xDim + 1;
        n[site_id].neigh_id[2]  = k + (j-xDim+nsp)%nsp + i%xDim + 1 ;
        n[site_id].neigh_id[3]  = k + (j-xDim+nsp)%nsp + i;
        n[site_id].neigh_id[4]  = k + (j-xDim+nsp)%nsp + (i-2+xDim)%xDim + 1;

        n[site_id].neigh_id[5]  = k + j + (i-2+xDim)%xDim + 1;
        n[site_id].neigh_id[6]  = k + (j+xDim)%nsp + (i-2+xDim)%xDim + 1;
        n[site_id].neigh_id[7]  = k + (j+xDim)%nsp + i;
        n[site_id].neigh_id[8]  = k + (j+xDim)%nsp + i%xDim + 1;
        
        //upper plane...
        //move the plane up and use the same scheme...
        n[site_id].neigh_id[9]  = (k-nsp+ns)%ns + j + i;
        n[site_id].neigh_id[10] = (k-nsp+ns)%ns + j + i%xDim + 1;
        n[site_id].neigh_id[11] = (k-nsp+ns)%ns + (j-xDim+nsp)%nsp + i%xDim + 1;

        n[site_id].neigh_id[12] = (k-nsp+ns)%ns + (j-xDim+nsp)%nsp + i;
        n[site_id].neigh_id[13] = (k-nsp+ns)%ns + (j-xDim+nsp)%nsp + (i-2+xDim)%xDim + 1;
        n[site_id].neigh_id[14] = (k-nsp+ns)%ns + j + (i-2+xDim)%xDim + 1;

        n[site_id].neigh_id[15] = (k-nsp+ns)%ns + (j+xDim)%nsp + (i-2+xDim)%xDim + 1;
        n[site_id].neigh_id[16] = (k-nsp+ns)%ns + (j+xDim)%nsp + i;
        n[site_id].neigh_id[17] = (k-nsp+ns)%ns + (j+xDim)%nsp + i%xDim + 1;
        
        //lower plane...
        n[site_id].neigh_id[18] = (k+nsp)%ns + j + i;
        n[site_id].neigh_id[19] = (k+nsp)%ns + j + i%xDim + 1;
        n[site_id].neigh_id[20] = (k+nsp)%ns + (j-xDim+nsp)%nsp + i%xDim + 1;

        n[site_id].neigh_id[21] = (k+nsp)%ns + (j-xDim+nsp)%nsp + i;
        n[site_id].neigh_id[22] = (k+nsp)%ns + (j-xDim+nsp)%nsp + (i-2+xDim)%xDim + 1;
        n[site_id].neigh_id[23] = (k+nsp)%ns + j + (i-2+xDim)%xDim + 1;

        n[site_id].neigh_id[24] = (k+nsp)%ns + (j+xDim)%nsp + (i-2+xDim)%xDim + 1;
        n[site_id].neigh_id[25] = (k+nsp)%ns + (j+xDim)%nsp + i;
        n[site_id].neigh_id[26] = (k+nsp)%ns + (j+xDim)%nsp + i%xDim + 1;  
      }
    }
  }
  /* 
  int i2;
  
  for( i2=1; i2<=ns; i2++){
    printf("site %d\n\n ", i2);
    printf("%d  %d  %d\n",   n[i2].neigh_id[4], n[i2].neigh_id[3], n[i2].neigh_id[2]);
    printf("%d  %d  %d\n",   n[i2].neigh_id[5], i2               , n[i2].neigh_id[1]);
    printf("%d  %d  %d\n\n", n[i2].neigh_id[6], n[i2].neigh_id[7], n[i2].neigh_id[8]);
    
    printf("%d  %d  %d\n",   n[i2].neigh_id[13], n[i2].neigh_id[12],n[i2].neigh_id[11]);
    printf("%d  %d  %d\n",   n[i2].neigh_id[14], n[i2].neigh_id[9], n[i2].neigh_id[10]);
    printf("%d  %d  %d\n\n", n[i2].neigh_id[15], n[i2].neigh_id[16],n[i2].neigh_id[17]);
    
    printf("%d  %d  %d\n",   n[i2].neigh_id[22], n[i2].neigh_id[21], n[i2].neigh_id[20]);
    printf("%d  %d  %d\n",   n[i2].neigh_id[23], n[i2].neigh_id[18], n[i2].neigh_id[19]);
    printf("%d  %d  %d\n\n", n[i2].neigh_id[24], n[i2].neigh_id[25], n[i2].neigh_id[26]);
  }
  */ 
}


void initialize(struct voxel *p, int ns){

  // This function reads micro.input as an input file and stores the spin number in 
  // structure voxel p[i].spin. After that, it generates cartesian coordinates (x, y, z)
  // for each voxel, starting from (0, 0, 0) to (xDim-1, yDim-1, zDim-1), and stores them
  // into p[i].x-y-zcoord.

  FILE *f1;
  int i, j, k, kk, l, x, y, z;

  char name[128];
  char trash[255];

  printf("Enter the name of input file (.dx):\n");
  scanf("%s", name);
 
  if( (f1=fopen(name, "r")) == NULL){
    printf("\nThe input file doesn't exist!\n");
    exit(1);
  }
 
  // Let's get rid of the header lines in micro.input...
  //fscanf(f1, "%d %d %d\n", &x, &y, &z);
  //fgets(trash, 126, f1);
  //fgets(trash, 126, f1);

  // Let's get rid of the 9 header lines in the xxx.dx produced by 3d coarsening code...
  for(kk=0; kk<9; kk++){
    fgets(trash, 255, f1);
  }

  p[0].spin = 0; // Point 0 is a garbage...

  // Let's input the spin numbers for each site from micro.input...
  for(i=1; i<=ns; i++){
    fscanf(f1, "%d",  &p[i].spin);
  }
 
  fclose(f1);

  // Copying the original spins...
  for(j=1; j<=ns; j++){
    p[j].newspin = p[j].spin;
  }

}


int get_number_grains (struct voxel *p, struct neighbor *n, 
                        int ns, int nsp, int xDim, int yDim, int zDim){

  int i, j, k, m, index, ii;
  int nucleus, head, tail, coin;
  int csite, cspin, nsite, nspin, chaser;
  int volume;
  int *burnt_list;
  int num_grain;
  int start;

  int sum;
  sum = 0;

  int tempsite;
  
  burnt_list = (int *)malloc(ns*sizeof(int));

  index = 0;

  tail = 0;
  head = 0;

  for (i=1; i<=ns; i++){
    
    csite = i;
    cspin = p[csite].spin;
    volume = 0;
    
    if(cspin>=0){
      
      nucleus = csite;
      p[nucleus].spin = cspin*(-1) - 10;//burn it...
      burnt_list[tail] = nucleus;
      volume = 1;
      
      start = tail;
      
      do{
        
        chaser = burnt_list[tail];
        
        for(m=1; m<=num_neigh; m++){
          
          nsite = n[chaser].neigh_id[m];
          nspin = p[nsite].spin;
          
          if(nspin == cspin){
            head = head + 1;
            burnt_list[head] = nsite;
            volume = volume + 1;
            p[nsite].spin = nspin*(-1) - 10;
          }
        }
        
        if(tail==head){
          coin = 0;
          tail = tail + 1;
          head = tail;
        }else{
          tail = tail + 1;
          coin = 1;
        }
        
      }while(coin);
      

      index++;
      
    }
      
  }
  
  free(burnt_list);

  // Recover the spin...
  for(ii=1; ii<=ns; ii++){
    p[ii].spin = p[ii].newspin;
  }
  
  return index;

}





void get_grains (struct voxel *p, struct neighbor *n, struct volume *g, 
                                 int ns, int nsp, int xDim, int yDim, int zDim){

  int i, j, k, m, index, ii;
  int nucleus, head, tail, coin;
  int csite, cspin, nsite, nspin, chaser;
  int volume;
  int *burnt_list;
  int start;

  float r3, r;
  
  burnt_list = (int *)malloc(ns*sizeof(int));

  index = 0;

  tail = 0;
  head = 0;

  for (i=1; i<=ns; i++){
    
    csite = i;
    cspin = p[csite].spin;
    volume = 0;
    
    if(cspin>=0){
      
      nucleus = csite;
      p[nucleus].spin = cspin*(-1) - 10;//burn it...
      burnt_list[tail] = nucleus;
      volume = 1;
      
      start = tail;
      
      do{
        
                chaser = burnt_list[tail];
        
                for(m=1; m<=num_neigh; m++){
          
                  nsite = n[chaser].neigh_id[m];
                  nspin = p[nsite].spin;
          
                  if(nspin == cspin){
                        head = head + 1;
                        burnt_list[head] = nsite;
                        volume = volume + 1;
                        p[nsite].spin = nspin*(-1) - 10;
                  }
                }
        
                if(tail==head){
                  coin = 0;
                  tail = tail + 1;
                  head = tail;
                }else{
                  tail = tail + 1;
                  coin = 1;
                }
        
      }while(coin);
      
      g[index].spin = cspin;
      g[index].size = volume;
      g[index].member = (int *)malloc(volume*sizeof(int));

      for(j=0; j<volume; j++){
                g[index].member[j] = burnt_list[start+j];
      }

      index++;
      
    }
      
  }
  
  free(burnt_list);

  // Recover the spin...
  for(ii=1; ii<=ns; ii++){
    p[ii].spin = p[ii].newspin;
  }
  

}



void randomize_grain_spin(struct voxel *p, struct volume *g, int ng){

  int i, j, ii, jj, kk;
  int *Rflag;
  int *new_spin;
  int newS;
  int index;
  int vol;
  int newID;
  int tsite;


  Rflag = (int *)malloc(ng*sizeof(int));
  new_spin = (int *)malloc(ng*sizeof(int));
  for(i=0; i<ng; i++){
    Rflag[i] = 0;
  }

  index = 0;
 
  do{

    newS = q_ran(ng);

    if(Rflag[newS-1]==0){
      new_spin[index] = newS;
      Rflag[newS-1] = 1;
      index++;
    }

  }while(index<ng);


  // Update the spin of the voxels according to the new grain information...
  for(ii=0; ii<ng; ii++){
    
    vol = g[ii].size;
    newID = new_spin[ii];
    
    // update spin of voxels...
    for(jj=0; jj<vol; jj++){
      tsite = g[ii].member[jj];
      p[tsite].newspin = newID;
    }

    g[ii].spin = newID;

  }

  free(Rflag);
  free(new_spin);

}




void get_output_dx(struct voxel *p, int ns, int xDim, int yDim, int zDim){

  // This function creates output microstructure based on q values in voxels...
  // The output file format is dx...(OpenDX input native file format...)

  FILE *f;
 
  char name1[128];
 
  printf("Enter the name of output dx file:\n");
  scanf("%s", name1);
 
  int i, j;

  if ((f = fopen(name1, "w")) == NULL){
    printf("Cannot create the file....\n ");
    exit(0);
  }
  
  fprintf(f, "object 1 class gridpositions counts %d %d %d\n", zDim+1, yDim+1, xDim+1);
  fprintf(f, "origin\t 0\t 0\t 0\n");
  fprintf(f, "delta \t 0\t 0\t 1\n");
  fprintf(f, "delta \t 0\t 1\t 0\n");
  fprintf(f, "delta \t 1\t 0\t 0\n");
  fprintf(f, "\n");

  fprintf(f, "object 2 class gridconnections counts %d %d %d\n", zDim+1, yDim+1, xDim+1);
  fprintf(f, "\n");

  fprintf(f, "object 3 class array type int rank 0 items %d data follows\n", ns); 
  for(i=1; i<ns+1; i++){
    fprintf(f,"%6d", p[i].newspin);
    if((i!=0)&&(i%20==0)){
      fprintf(f, "\n");
    }
  }
  fprintf(f, "\n");
  fprintf(f, "attribute \"dep\" string \"connections\"\n");
  fprintf(f, "\n");

  fprintf(f, "object \"3d MC Coarsening\" class field\n");
  fprintf(f, "component  \"positions\"    value 1\n");
  fprintf(f, "component  \"connections\"  value 2\n");
  fprintf(f, "component  \"data\"         value 3\n");
  fprintf(f, "\n");
  
  fprintf(f, "end\n");

  fclose(f);


  return;
} 




/*
   double ran3(long *idum)
   ============================================================================
   Numerical Recipies random number generator; returns random number in the
   range 0<=n<1.

   idum: seed
   ============================================================================
*/

double ran3(long *idum)
{
        static int inext,inextp;
        static long ma[56];
        static int iff=0;
        long mj,mk;
        int i,ii,k;

        if (*idum < 0 || iff == 0) { 
                iff=1;
                mj=labs(MSEED-labs(*idum)); 
                mj %= MBIG;
                ma[55]=mj;
                mk=1;
                for (i=1;i<=54;i++) {
                        ii=(21*i) % 55;
                        ma[ii]=mk;
                        mk=mj-mk;
                        if (mk < MZ) mk += MBIG;
                        mj=ma[ii];
                }
                for (k=1;k<=4;k++)
                        for (i=1;i<=55;i++) {
                                ma[i] -= ma[1+(i+30) % 55];
                                if (ma[i] < MZ) ma[i] += MBIG;
                        }
                inext=0;
                inextp=31;
                *idum=1;
        }

        if (++inext == 56) inext=1;
        if (++inextp == 56) inextp=1;
        mj=ma[inext]-ma[inextp];
        if (mj < MZ) mj += MBIG; 
        ma[inext]=mj;
        return mj*FAC;
}


int q_ran(int q){
   
   // This function gives a random integer from 1 to q. (1 <= n <= q)

  static int itemp;
  static long seed1;

   //   static long ltime;
   //  ltime = time(NULL);

   seed1 = 1;//*ltime; // just a trick to generate different seed...
   itemp = (int)(q*ran3(&seed1) + 1);
   
   //printf("%6.4f  %5d\n", ran_temp, itemp);
   
   return itemp;
}
