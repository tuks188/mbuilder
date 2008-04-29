// Written by Joe Fridy

//Used in conjunction with mdfextract.sh shell script
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


/********************** GLOBALS ****************************/
#define PI 3.14159
#define CB (pow(0.75*((PI/4)-sin(PI/4)),(double)1/3))
#define XMAX 10
#define NCELLS (XMAX*XMAX*XMAX)
#define AMIN (5*(PI/180))

#define MIN(a,b)        ((a) < (b) ? (a) : (b))
#define MAX(a,b)        ((a) > (b) ? (a) : (b))
#define MINMIN(a,b,c)   ((a) < (b) ? MIN (a,c) : MIN (b,c))
#define MAXMAX(a,b,c)   ((a) > (b) ? MAX (a,c) : MAX (b,c)) 
/***********************************************************/

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
typedef struct {
	int x,y;
} TPoint;

void OutputMdf(char fname[100], double *mdf);
void NormalizeMdf(double *mdf);
void AddToMdf(TSymOp sym, double *mdf, matrix g1, matrix g2);
void FillMdf(TSymOp sym, TPoint oim_size, TEuler **data, double *mdf);
double Misort(TSymOp sym, matrix g1, matrix g2);
void GetData(char oimf_name[100], TPoint oim_size, TEuler **data);
TPoint GetOimXY(char oimf_name[100]);
void ReadHead(FILE *ipf);
int VectToCell(vector d);
void MinAngMatrix(TSymOp sym, matrix g);
TAxAng GToAA(matrix g);
void Transpose(matrix g, matrix gt);
void MM(matrix g1, matrix g2, matrix g3);
void EToG(TEuler euler, matrix g);
double acos2(double ca);
TSymOp LoadSym(char fname[100]);

int main()
{
  double *mdf;
  int i,snum,s;
  TSymOp sym;
  TPoint oim_size;
  TEuler **data;
  char infname[100],outfname[100];
  matrix gcs;
  FILE *cf, *outf;

  strcpy(outfname,"evmdf.txt");
  printf("reading in data\n");
  sym = LoadSym("symop.txt");
  cf=fopen("ctrlfile.txt","r");
  fscanf(cf,"%d",&snum);

  for(s=0;s<snum;s++){
	  fscanf(cf,"%s",infname);
	  fscanf(cf,"%lf%lf%lf",&gcs[0][0],&gcs[0][1],&gcs[0][2]);
	  fscanf(cf,"%lf%lf%lf",&gcs[1][0],&gcs[1][1],&gcs[1][2]);
	  fscanf(cf,"%lf%lf%lf",&gcs[2][0],&gcs[2][1],&gcs[2][2]);
	  
	  mdf = (double*)calloc(NCELLS, sizeof(double));
	  oim_size = GetOimXY(infname);
	  data = (TEuler**)calloc(oim_size.x,(sizeof(TEuler*)));
	  for (i=0;i<oim_size.x;i++){
		data[i]=(TEuler*)calloc(oim_size.y,(sizeof(TEuler)));
	  }
	  GetData(infname,oim_size,data);
	  FillMdf(sym,oim_size,data,mdf);
  }
  NormalizeMdf(mdf);
  //outf=fopen("mdfoutname.txt","r");
  //outf=fopen("evmdf.txt","r");
  //fscanf(outf,"%s",outfname);
  OutputMdf(outfname,mdf);
  free(mdf);

  printf("MDF_anneal is done! \n");
}

void OutputMdf(char fname[100], double *mdf)
{
	FILE *fp;
	int i;

	fp = fopen(fname,"w");
	for(i=0;i<NCELLS;i++){
		fprintf(fp,"%d\t%lf\n",i,mdf[i]);
	}
	fclose(fp);
}

void NormalizeMdf(double *mdf)
{
	int i;
	double sum;
	
	sum = 0.0;
	for(i=0;i<NCELLS;i++){
		sum = sum + mdf[i];
	}
	for(i=0;i<NCELLS;i++){
		mdf[i] = mdf[i]/sum;
	}
}

void AddToMdf(TSymOp sym, double *mdf, matrix g1, matrix g2)
{
	matrix g,gt;
	int cell;
	vector d;
	double c;
	TAxAng aa;

	Transpose(g2,gt);
	MM(g1,gt,g);
	MinAngMatrix(sym,g);
	aa = GToAA(g);
	c = pow(0.75*(aa.ang-sin(aa.ang)),(double)1/3);
	aa.n[0] = fabs(aa.n[0]);
	aa.n[1] = fabs(aa.n[1]);
	aa.n[2] = fabs(aa.n[2]);
	d[0] = MAXMAX(aa.n[0],aa.n[1],aa.n[2]);
	d[2] = MINMIN(aa.n[0],aa.n[1],aa.n[2]);
	d[1] = aa.n[0]+aa.n[1]+aa.n[2]-d[0]-d[2];
	d[0] = c*d[0];
	d[1] = c*d[1];
	d[2] = c*d[2];
	cell = VectToCell(d);
	mdf[cell] = mdf[cell] + 1;
}

void FillMdf(TSymOp sym, TPoint oim_size, TEuler **data, double *mdf)
{
	int i,j,a,b;
	matrix g1,g2;
	double ang;
	
	for(j=0;j<oim_size.y;j++){
		for(i=0;i<oim_size.x;i++){
			EToG(data[i][j],g1);
			for(b=0;b<=1;b++)
				for(a=0;a<=1;a++)
					if((i+a<oim_size.x)&&(j+b<oim_size.y))
						if(a!=b){
							EToG(data[i+a][j+b],g2);
							ang = Misort(sym,g1,g2);
							if(ang>=AMIN) AddToMdf(sym,mdf,g1,g2);
						}
		}
	}
}

double Misort(TSymOp sym, matrix g1, matrix g2)
{
	int ii,i,j,k;
	double omeg,gdist,ca,angle;

	gdist = -1;
	for(ii=0;ii<sym.num;ii++){
		omeg = 0;
		for(i=0;i<3;i++)
			for(j=0;j<3;j++)
				for(k=0;k<3;k++)
					omeg = omeg+sym.op[ii][i][j]*g1[i][k]*g2[j][k];
		if(omeg>gdist) gdist = omeg;
	}
	ca = (gdist-1)/2;
    angle = acos2(ca);
    return(angle);
}

void GetData(char oimf_name[100], TPoint oim_size, TEuler **data)
{
	FILE *ipf;
	int i,j;
	double x,y,iq,ci,t1,t2;
	// C++ Vars
	double phi1,phi,phi2;
	ifstream input;
	string s;

	input.open(oimf_name);

        while (input.peek()=='#') {
                char buffer[500];
                input.getline(buffer,500);
        }

	for(j=0;j<oim_size.y;j++){
	  for(i=0;i<oim_size.x;i++){		  
	    getline(input,s);
	    stringstream ss;
	    ss<<s;
	    ss >> data[i][j].a1 >> data[i][j].a >> data[i][j].a2 ;
	    //ss >> phi1 >> phi >> phi2;
	    //data[i][j].a1 = phi1;
	    //data[i][j].a = phi;
	    //data[i][j].a2 = phi2;
	  }
	}
		  
	input.close();
}


TPoint GetOimXY(char oimf_name[100])
{
	TPoint oim_size;
	int i,j;
	double x,y,iq,ci,t1,t2,x1,y1,a1,a,a2;
	double ang[5];
	double xcoor,ycoor;
	ifstream input;
	string s;
	
	i = j = 1;
	x1 = y1 = 0;

	input.open(oimf_name);

        while (input.peek()=='#') {
                char buffer[500];
                input.getline(buffer,500);
        }

        while (!input.eof()) {
                getline(input,s);
		stringstream ss;
                ss<<s;
                ss>>ang[0]>>ang[1]>>ang[2]>>ang[3]>>ang[4];
		if (ang[3]>x1) {
			x1 = ang[3];
			i++;
		}
		if (ang[4]>y1) {
			y1 = ang[4];
			j++;
		}
        }

	input.close();
	oim_size.x = i;
	oim_size.y = j;
	cout << oim_size.x << " "<< oim_size.y << endl;
	return(oim_size);
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
	
	ix = (int)floor(XMAX*(d[0])/(CB));
	iy = (int)floor(XMAX*(d[1])/(CB));
	iz = (int)floor(XMAX*(d[2])/(CB));
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


