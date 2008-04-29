// Written by Joe Fridy

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

// CXX HEADERS
#include<iostream>
#include<fstream>
#include<cmath>
#include<sstream>
using namespace std;


#define PI 3.14159
#define CB (pow(0.75*((PI/4)-sin(PI/4)),(double)1/3))
#define XMAX 10
#define NCELLS (XMAX*XMAX*XMAX)

#define MIN(a,b)        ((a) < (b) ? (a) : (b))
#define MAX(a,b)        ((a) > (b) ? (a) : (b))
#define MINMIN(a,b,c)   ((a) < (b) ? MIN (a,c) : MIN (b,c))
#define MAXMAX(a,b,c)   ((a) > (b) ? MAX (a,c) : MAX (b,c)) 

typedef double vector[3];
typedef double matrix[3][3];
typedef struct {
	matrix *op;
	int num;
}TSymOp;
typedef struct {
	double a1,a,a2;
}TEuler;
typedef struct{
	double ang;
	vector n;
} TAxAng;

int VectToCell(vector d);
void MinAngMatrix(TSymOp sym, matrix g);
TAxAng GToAA(matrix g);
void Transpose(matrix g, matrix gt);
void MM(matrix g1, matrix g2, matrix g3);
void EToG(TEuler euler, matrix g);
double acos2(double ca);
TSymOp LoadSym(char fname[100]);
void ReadHead(FILE *ipf);
void MM(matrix g1, matrix g2, matrix g3);

int main()
{
  int i,cell,s,snum;
  TAxAng aa;
  FILE *fp,*cf;
  double *odf,c,sum,x,y;
  matrix g,gcs,gca;
  vector d;
  TEuler e1;
  TSymOp sym;
  char fname[100];
  ifstream input;
  string sxx;

  sym = LoadSym("symop.txt");

  odf = (double*)calloc(NCELLS, sizeof(double));

  /*fill odf from data file*/
  printf("odfextract_final.cxx reading in data\n");
  cf=fopen("ctrlfile.txt","r");
  fscanf(cf,"%d",&snum);
  printf("#files=%d \n",snum);
  for(s=0;s<snum;s++){
	  fscanf(cf,"%s",fname);
          printf("FILE=%s \n",fname);
	  fscanf(cf,"%lf%lf%lf",&gcs[0][0],&gcs[0][1],&gcs[0][2]);
	  fscanf(cf,"%lf%lf%lf",&gcs[1][0],&gcs[1][1],&gcs[1][2]);
	  fscanf(cf,"%lf%lf%lf",&gcs[2][0],&gcs[2][1],&gcs[2][2]);
	  
	  input.open(fname);

	  while (input.peek()=='#') {
	    char buffer[500];
	    input.getline(buffer,500);
	  }
	  while(!input.eof()){
	    getline(input,sxx);
	    stringstream ss;
	    ss<<sxx;
	    ss >> e1.a1 >> e1.a >> e1.a2 ;
	    EToG(e1,gca);
	    MM(gcs,gca,g);
	    MinAngMatrix(sym,g);
	    aa = GToAA(g);
	    //printf(" %f   %f   %f  %f  \n",aa.ang,aa.n[0],aa.n[1],aa.n[2]);
	    c = pow(0.75*(aa.ang-sin(aa.ang)),(double)1/3);
	    d[0] = c*aa.n[0];
	    d[1] = c*aa.n[1];
	    d[2] = c*aa.n[2];
	    cell = VectToCell(d);
	    odf[cell] = odf[cell] + 1;
	  }
	  input.close();
  }
  fclose(cf);

  /*load cell volumes and normalize distribution*/
  printf("normalizing odf\n");
  sum = 0.0;
  for(i=0;i<NCELLS;i++){
	  sum = sum + odf[i];
  }

  //char odfoutput[100];
  //ifstream inf1;
  //inf1.open("odfoutname.txt");
  //inf1 >> odfoutput;
  //inf1.close();
  //cout << "New odf file= " << odfoutput << endl;

  fp = fopen("evodf.txt","w");
  for(i=0;i<NCELLS;i++){
	  fprintf(fp,"%d\t%0.8lf\n",i,odf[i]/sum);
  }
  fclose(fp);
  free(odf);
  printf("odf_anneal is done! \n");
}

void ReadHead(FILE *ipf)
{
	char c, s[100];

	while((c=fgetc(ipf))==35){ 
		fscanf(ipf,"%[^\n]",s); 
		c=fgetc(ipf); 
    }
	fseek(ipf,-1,1);
}

int VectToCell(vector d)
{
	int ix,iy,iz,cell;
	
	ix = (int)floor(XMAX*(d[0]+CB)/(2*CB));
	iy = (int)floor(XMAX*(d[1]+CB)/(2*CB));
	iz = (int)floor(XMAX*(d[2]+CB)/(2*CB));
	if(ix==XMAX) ix = XMAX - 1;
	if(iy==XMAX) iy = XMAX - 1;
	if(iz==XMAX) iz = XMAX - 1;
	cell = XMAX*XMAX*iz+XMAX*iy+ix;
	return(cell);
}

void MinAngMatrix(TSymOp sym, matrix g)
{
	int ii,i,j,iispec;
	matrix gf;
	double gdist,tr;

	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			gf[i][j] = g[i][j];
		}
	}
	gdist = -1.0;
	for(ii=0;ii<sym.num;ii++){
		MM(sym.op[ii],gf,g);
		tr = g[0][0] + g[1][1] + g[2][2];
		if(tr>gdist){
			gdist = tr;
			iispec = ii;
		}
	}
	MM(sym.op[iispec],gf,g);
}

TAxAng GToAA(matrix g)
{
	int i,j,k;
	double tr=0.0;
	TAxAng ax;
	double pm[3][3][3] = {{{0,0,0},{0,0,1},{0,-1,0}},{{0,0,-1},{0,0,0},{1,0,0}},{{0,1,0},{-1,0,0},{0,0,0}}};

	for(i=0;i<3;i++){
		tr = tr + g[i][i];
		ax.n[i] = 0.0;
	}
	ax.ang = acos2((tr-1)/2);
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			for(k=0;k<3;k++){
				ax.n[i] = ax.n[i] - pm[i][j][k]*g[j][k];
			}
		}
	}
	if(ax.ang!=0.0){
		for(i=0;i<3;i++){
			ax.n[i] = ax.n[i]/(2*sin(ax.ang));
		}
	}else{
		for(i=0;i<3;i++){
			ax.n[i] = 0.0;
		}
	}
	return(ax);
}

void Transpose(matrix g, matrix gt)
{
	int i,j;

	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			gt[i][j]=g[j][i];
}
     
void MM(matrix g1, matrix g2, matrix g3)
{
  int i,j,k;
  
  for(i=0;i<3;i++)
    for(j=0;j<3;j++){
      g3[i][j] = 0.0;
      for(k=0;k<3;k++)
		g3[i][j] = g3[i][j] + g1[i][k]*g2[k][j];
    }
}
       
void EToG(TEuler euler, matrix g)
{
  g[0][0] = cos(euler.a1)*cos(euler.a2)-sin(euler.a1)*sin(euler.a2)*cos(euler.a);
  g[0][1] = sin(euler.a1)*cos(euler.a2)+cos(euler.a1)*sin(euler.a2)*cos(euler.a);
  g[0][2] = sin(euler.a2)*sin(euler.a);
  g[1][0] = -cos(euler.a1)*sin(euler.a2)-sin(euler.a1)*cos(euler.a2)*cos(euler.a);
  g[1][1] = -sin(euler.a1)*sin(euler.a2)+cos(euler.a1)*cos(euler.a2)*cos(euler.a);
  g[1][2] = cos(euler.a2)*sin(euler.a);
  g[2][0] = sin(euler.a1)*sin(euler.a);
  g[2][1] = -cos(euler.a1)*sin(euler.a);
  g[2][2] = cos(euler.a);
}
      
double acos2(double ca)
{
  if(ca < -1) ca = -1;
  if(ca > 1) ca = 1;
  return acos(ca);
}

TSymOp LoadSym(char fname[100])
{
  FILE *ipf;
  TSymOp sym; 		
  int ii,i,j;

  ipf = fopen(fname,"r");
  fscanf(ipf,"%d",&sym.num);
  sym.op = (matrix*)calloc(sym.num, sizeof(matrix));
  for(ii=0;ii<sym.num;ii++)
    for(i=0;i<3;i++)
      for(j=0;j<3;j++)
		fscanf(ipf,"%lf",&sym.op[ii][i][j]);
  fclose(ipf);
  return(sym);
} 

