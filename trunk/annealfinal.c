/* Program to assign orientations to a Fridy style microstructural geometry file
such that the ODF and MDF matched 
================================================================================
INPUT FILES:
	cellIdealization.xml: "Fridy style" geometry file
	evodf.txt: ODF in homochoric parameterization 
	evmdf.txt: MDF in homochoric parameterization
	symop.txt: symmetry operators
OUTPUT FILES
	orts.txt: assigned grain orientations (Rodrigues vectors)
	rodf.txt: final odf (Homochoric)
	rmdf.txt: final mdf (Homochoric)
	anneal.log: log file
================================================================================
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

/*
   ============================================================================
   constants
   ============================================================================
*/

#define TFACTR 0.8 // temperature reduction
#define NSTEPS 50 // annealing steps
#define MINERROR 0.010 // minimum error which the program exits regardless of
	            // number of steps that have been performed.

#define XMAX 10 // divisions of both ODF and MDF bounding boxes in each dimension
#define NCELLS (XMAX*XMAX*XMAX) // total cells for both ODF and MDF
#define CB (pow(0.75*((PI/4)-sin(PI/4)),(double)1/3)) // extent of bounding boxes
#define EPS 10E-7 
#define PI 3.14159

// constants for procedure ran3
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

/*
   ============================================================================
   macros
   ============================================================================
*/

static double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
#define MIN(a,b)        ((a) < (b) ? (a) : (b))
#define MAX(a,b)        ((a) > (b) ? (a) : (b))
#define MINMIN(a,b,c)   ((a) < (b) ? MIN (a,c) : MIN (b,c))
#define MAXMAX(a,b,c)   ((a) > (b) ? MAX (a,c) : MAX (b,c)) 

/*
   ============================================================================
   data structures
   ============================================================================
*/


typedef double vector[3];
typedef double matrix[3][3];
typedef struct {
	double a1,a,a2;
}TEuler;
typedef struct{
	int neighbor;
	double area;
	vector n;
}TPatch;
typedef struct{
	double vol;
	int npatch;
	vector ort;
	TPatch *patch;
}TRegion;
typedef struct {
	matrix *op;
	int num;
}TSymOp;
typedef struct{
	double ang;
	vector n;
} TAxAng;

/*
   ============================================================================
   function definitions
   ============================================================================
*/

double *LoadDF(char fname[100]);
TRegion *LoadGeometryData(char fname[100], int *nregions);
double *InitODF(TSymOp sym, TRegion *data, int nregions, double *tvol);
double *InitMDF(TSymOp sym, TRegion *data, int nregions, double *tarea);
void RandE(TSymOp sym, vector d);
int metrop(double de, double t);
double ran3(long *idum);
void DFSave(double *odf, char fname[100]);
void MinAngMatrix(TSymOp sym, matrix g);
TAxAng GToAA(matrix g);
void Transpose(matrix g, matrix gt);
void MM(matrix g1, matrix g2, matrix g3);
void EToG(TEuler euler, matrix g);
double acos2(double ca);
TSymOp LoadSym(char fname[100]);
int OVectToCell(vector d);
int MVectToCell(vector d);
void VectToG(vector n, matrix g);
void MisVect(TSymOp sym, matrix g1, matrix g2, vector d);
void VConvert(vector r, vector d);

/*
   ============================================================================
   main program
   ============================================================================
*/

int main()
{
	TRegion *data;
	double *odf,*mdf,*codf,*cmdf,*otemp,*mtemp;
	int ans,nover,nlimit,nsucc,nregions;
	int i,n,nn,id,id1,id2;
	long idum;
	double error,de,t,tvol,tarea;
	FILE *ofp,*ifp;
	time_t now;
	TSymOp sym;
	vector df,tmp1,tmp2,dp;
	int io,fo,kk,icell,fcell;
	double atemp;
	matrix gi,gf,gn;
	vector dd;

	// open log file
	ifp = fopen("anneal.log","w");
	now = time(NULL);
	fprintf(ifp,"%s\n",ctime(&now));
	fflush(ifp);

	// load data
	sym = LoadSym("symop.txt");

	fprintf(ifp,"reading in grain geometry data ...\n");
	data = LoadGeometryData("cellIdealization.xml",&nregions);
	fprintf(ifp,"%d grains in microstructure\n",nregions);

	fprintf(ifp,"reading in ODF data ...\n");
	odf = LoadDF("evodf.txt");

	fprintf(ifp,"reading in MDF data ...\n");
	mdf = LoadDF("evmdf.txt");
	fflush(ifp);
	
	// initialize
	fprintf(ifp,"getting initial ODF ...\n");
	codf = InitODF(sym,data,nregions,&tvol);

	fprintf(ifp,"getting initial MDF ...\n");
	cmdf = InitMDF(sym,data,nregions,&tarea);
	fflush(ifp);


	/*initial error*/
	error = 0.0;
	for(i=0;i<NCELLS;i++){
		error = error + SQR(codf[i]-odf[i]) + SQR(cmdf[i]-mdf[i]);
	}
	fprintf(ifp,"error = %0.6lf\n\n",error);
	
	/*start the anneal*/
	fprintf(ifp,"starting simulated annealing\n");
	now = time(NULL);
	fprintf(ifp,"%s\n",ctime(&now));
	fflush(ifp);

	nover = 100*nregions; //maximum number of iterations per temp step
	nlimit = 10*nregions; //maximum number of successful iterations per temp step
	t = 0.00001; //initial annealing temp
	idum = -1;
	for(n=1;n<=NSTEPS;n++){
 		nsucc = 0;
	
		for(nn=1;nn<=nover;nn++){

			otemp = (double*)calloc(NCELLS, sizeof(double));
			mtemp = (double*)calloc(NCELLS, sizeof(double));

			if(ran3(&idum)<=0.5){ /*change*/

				/*assign random grain and new orientation*/
				id = (int)floor(nregions*ran3(&idum));
				RandE(sym,df);

				/*calculate de from changing orts*/
				VConvert(data[id].ort,tmp1);
				VConvert(df,tmp2);
				io = OVectToCell(tmp1);
				fo = OVectToCell(tmp2);
				otemp[io] = otemp[io] - data[id].vol;
				otemp[fo] = otemp[fo] + data[id].vol;
		
				/*calculate de from changing misorts*/
				VectToG(data[id].ort,gi);
				VectToG(df,gf);
				for(kk=0;kk<data[id].npatch;kk++){
					atemp = data[id].patch[kk].area;
					VectToG(data[data[id].patch[kk].neighbor].ort,gn);

					MisVect(sym,gi,gn,dd);
					VConvert(dd,dp);
					icell = MVectToCell(dp);
					MisVect(sym,gf,gn,dd);
					VConvert(dd,dp);
					fcell = MVectToCell(dp);

					mtemp[icell] = mtemp[icell] - atemp;
					mtemp[fcell] = mtemp[fcell] + atemp;
				}
	
				de = 0.0;
				for(i=0;i<NCELLS;i++){
					if(otemp[i]!=0.0){
						de = de - SQR(codf[i]-odf[i]) + SQR(codf[i]+otemp[i]-odf[i]);
					}
					if(mtemp[i]!=0.0){
						de = de - SQR(cmdf[i]-mdf[i]) + SQR(cmdf[i]+mtemp[i]-mdf[i]);
					}
				}

				ans = metrop(de,t);
				/* if true - make changes */
				if(ans){
					++nsucc;
					error = error + de;

					for(i=0;i<NCELLS;i++){
						if(otemp[i]!=0.0){
							codf[i] = codf[i] + otemp[i];
						}
						if(mtemp[i]!=0.0){
							cmdf[i] = cmdf[i] + mtemp[i];
						}
					}
					data[id].ort[0] = df[0];
					data[id].ort[1] = df[1];
					data[id].ort[2] = df[2];
			
				}

			} else { /*swap*/

				id1 = (int)floor(nregions*ran3(&idum));
				id2 = (int)floor(nregions*ran3(&idum));

				VConvert(data[id1].ort,tmp1);
				VConvert(data[id2].ort,tmp2);
				io = OVectToCell(tmp1);
				fo = OVectToCell(tmp2);
				otemp[io] = otemp[io] - data[id1].vol;
				otemp[fo] = otemp[fo] + data[id1].vol;
				otemp[fo] = otemp[fo] - data[id2].vol;
				otemp[io] = otemp[io] + data[id2].vol;

				/*calculate de from changing misorts*/
				VectToG(data[id1].ort,gi);
				VectToG(data[id2].ort,gf);
				for(kk=0;kk<data[id1].npatch;kk++){
					if(id2 != data[id1].patch[kk].neighbor){
						atemp = data[id1].patch[kk].area;
						VectToG(data[data[id1].patch[kk].neighbor].ort,gn);

						MisVect(sym,gi,gn,dd);
						VConvert(dd,dp);
						icell = MVectToCell(dp);
						MisVect(sym,gf,gn,dd);
						VConvert(dd,dp);
						fcell = MVectToCell(dp);

						mtemp[icell] = mtemp[icell] - atemp;
						mtemp[fcell] = mtemp[fcell] + atemp;
					}
				}
				VectToG(data[id2].ort,gi);
				VectToG(data[id1].ort,gf);
				for(kk=0;kk<data[id2].npatch;kk++){
					if(id1 != data[id2].patch[kk].neighbor){
						atemp = data[id2].patch[kk].area;
						VectToG(data[data[id2].patch[kk].neighbor].ort,gn);

						MisVect(sym,gi,gn,dd);
						VConvert(dd,dp);
						icell = MVectToCell(dp);
						MisVect(sym,gf,gn,dd);
						VConvert(dd,dp);
						fcell = MVectToCell(dp);

						mtemp[icell] = mtemp[icell] - atemp;
						mtemp[fcell] = mtemp[fcell] + atemp;
					}
				}

				de = 0.0;
				for(i=0;i<NCELLS;i++){
					if(otemp[i]!=0.0){
						de = de - SQR(codf[i]-odf[i]) + SQR(codf[i]+otemp[i]-odf[i]);
					}
					if(mtemp[i]!=0.0){
						de = de - SQR(cmdf[i]-mdf[i]) + SQR(cmdf[i]+mtemp[i]-mdf[i]);
					}
				}
			
				ans = metrop(de,t);
				/* if true - make changes */
				if(ans){

					++nsucc;
					error = error + de;

					for(i=0;i<NCELLS;i++){
						if(otemp[i]!=0.0){
							codf[i] = codf[i] + otemp[i];
						}
						if(mtemp[i]!=0.0){
							cmdf[i] = cmdf[i] + mtemp[i];
						}
					}

					df[0] = data[id2].ort[0];
					df[1] = data[id2].ort[1];
					df[2] = data[id2].ort[2];

					data[id2].ort[0] = data[id1].ort[0];
					data[id2].ort[1] = data[id1].ort[1];
					data[id2].ort[2] = data[id1].ort[2];

					data[id1].ort[0] = df[0];
					data[id1].ort[1] = df[1];
					data[id1].ort[2] = df[2];
				}
			}

			free(otemp);
			free(mtemp);
			if(nsucc >= nlimit) break;
		}

		//step statistics to log file
		fprintf(ifp,"\n %s %d %s %d\n","Step ",n," of ",NSTEPS);
		fprintf(ifp,"%s %10.6lf %s %12.8lf \n","T = ",t,"error = ",error);
		fprintf(ifp,"Successful moves: %6d\n",nsucc);
		now = time(NULL);
		fprintf(ifp,"%s\n",ctime(&now));
		fflush(ifp);

		//write out data after each step (crash protection)
		ofp = fopen("orts.txt","w");
		for(i=0;i<nregions;i++){
			fprintf(ofp,"%d\t%0.8lf\t%0.8lf\t%0.8lf\n",i,data[i].ort[0],data[i].ort[1],data[i].ort[2]);
		}
		fclose(ofp);
		fprintf(ifp,"outputting results to file ...\n");
		ofp = fopen("rodf.txt","w");
		for(i=0;i<NCELLS;i++){
			fprintf(ofp,"%d\t%0.8lf\n",i,codf[i]);
		}
		fclose(ofp);
		ofp = fopen("rmdf.txt","w");
		for(i=0;i<NCELLS;i++){
			fprintf(ofp,"%d\t%0.8lf\n",i,cmdf[i]);
		}
		fclose(ofp);

		t *= TFACTR;
		if(nsucc==0) return(1);
                //If combined MDF & ODF error is less than MINERROR
		//Exit this loop.
                if(error<MINERROR) break;
	}
	
	free(odf);
	free(mdf);
	free(data);
	free(codf);
	free(cmdf);

	fprintf(ifp,"done!\n");
	now = time(NULL);
	fprintf(ifp,"%s\n",ctime(&now));
	fclose(ifp);

	return(0);
}

/*
   TRegion *LoadGeometryData(char fname[100], int *nregions)	
   ============================================================================
   Loads a "Fridy style" geometry file and returns a TRegion structure that 
   contains the necessary geometric information

   fname: Fridy style geometry file
   nregions: number of grains
   ============================================================================
*/
TRegion *LoadGeometryData(char fname[100], int *nregions)
{
	FILE *ifp;
	char s[100],s1[100],s2[100];
    int i,j,count=0;
	TRegion *data;
	
    ifp = fopen(fname,"r");
    if(!ifp){printf("Error opening %s\n",fname); exit(1);}
    fscanf(ifp,"%*s %d %*s",nregions);

    data = (TRegion*)calloc(*nregions, sizeof(TRegion));
    for(i=0;i<*nregions;i++){
		fscanf(ifp,"%*s");
	    fscanf(ifp,"%*s %*s %*s");
	    fscanf(ifp,"%*s %s %*s",s);
	    data[i].vol = atof(s);
	
		fscanf(ifp,"%*s %*s %*s %*s %*s");

		fscanf(ifp,"%*s");
		fscanf(ifp,"%*s %*s %*s %*s %*s");
		fscanf(ifp,"%*s %*s %*s %*s %*s");
		fscanf(ifp,"%*s %*s %*s %*s %*s");
		fscanf(ifp,"%*s %*s %*s %*s %*s");
		fscanf(ifp,"%*s %*s %*s %*s %*s");
		fscanf(ifp,"%*s %*s %*s %*s %*s");
		fscanf(ifp,"%*s");

		fscanf(ifp,"%*s %d %*s",&data[i].npatch);
		count = count + data[i].npatch;
		data[i].patch = (TPatch*)calloc(data[i].npatch, sizeof(TPatch));
		for(j=0;j<data[i].npatch;j++){
		  fscanf(ifp,"%*s");
		  fscanf(ifp,"%*s %d %*s",&data[i].patch[j].neighbor);
		  fscanf(ifp,"%*s %s %s %s %*s",s,s1,s2);
		  data[i].patch[j].n[0] = atof(s);
		  data[i].patch[j].n[1] = atof(s1);
		  data[i].patch[j].n[2] = atof(s2);
          fscanf(ifp,"%*s %s %*s",s);
		  data[i].patch[j].area = atof(s);
		  fscanf(ifp,"%*s %*s %*s");
		  fscanf(ifp,"%*s");
		}
		fscanf(ifp,"%*s");
	}
	fclose(ifp);
	return(data);
}

/*
   double *LoadDF(char fname[100])	
   ============================================================================
   Loads a "Saylor style" orientation or misorientation distribution file and 
   returns the distribution.

   fname: Saylor style distribution file
   ============================================================================
*/
double *LoadDF(char fname[100])
{
	FILE *ifp;
	int i,temp;
	double *odf;

	odf = (double*)calloc(NCELLS, sizeof(double));
	ifp = fopen(fname,"r");
        if(!ifp){printf("Error opening %s\n",fname); exit(1);}
	for(i=0;i<NCELLS;i++){
		fscanf(ifp,"%d%lf",&temp,&odf[i]);
	}
	fclose(ifp);
	return(odf);
}

/*
   double *InitODF(TSymOp sym, TRegion *data, int nregions, double *tvol)
   ============================================================================
   Initializes the starting ODF by assigning random orientations to each grain 
   in the structure; returns an ODF.

   sym: symmetry operators
   *data: geometry structure
   nregions: number of grains
   *tvol: total volume
   ============================================================================
*/
double *InitODF(TSymOp sym, TRegion *data, int nregions, double *tvol)
{
	double *codf;
	int i,cell;
	vector d,dp;

	codf = (double*)calloc(NCELLS, sizeof(double));
	*tvol = 0.0;
	for(i=0;i<nregions;i++){
		RandE(sym,d);
		data[i].ort[0] = d[0];
		data[i].ort[1] = d[1];
		data[i].ort[2] = d[2];
		VConvert(d,dp);
		cell = OVectToCell(dp);
		*tvol = *tvol + data[i].vol;
		codf[cell] = codf[cell] + data[i].vol;
	}
	for(i=0;i<NCELLS;i++){
		codf[i] = codf[i]/(*tvol);
	}
	for(i=0;i<nregions;i++){
		data[i].vol = data[i].vol/(*tvol);
	}
	return(codf);
}

/*
   double *InitMDF(TSymOp sym, TRegion *data, int nregions, double *tarea)
   ============================================================================
   Calculates the starting MDF and returns it.

   sym: symmetry operators
   *data: geometry structure
   nregions: number of grains
   *tarea: total area
   ============================================================================
*/
double *InitMDF(TSymOp sym, TRegion *data, int nregions, double *tarea)
{
	double *cmdf;
	int i,j,cell;
	matrix g1,g2;
	vector d,dp;

	cmdf = (double*)calloc(NCELLS, sizeof(double));
	*tarea = 0.0;
	for(i=0;i<nregions;i++){
		VectToG(data[i].ort,g1);
		for(j=0;j<data[i].npatch;j++){
			VectToG(data[data[i].patch[j].neighbor].ort,g2);
			MisVect(sym,g1,g2,d);
			VConvert(d,dp);
			cell = MVectToCell(dp);	
			cmdf[cell] = cmdf[cell] + data[i].patch[j].area;
			*tarea = *tarea + data[i].patch[j].area;
		}
	}
	for(i=0;i<NCELLS;i++){
		cmdf[i] = cmdf[i]/(*tarea);
	}
	for(i=0;i<nregions;i++){
		for(j=0;j<data[i].npatch;j++){
			data[i].patch[j].area = 2*data[i].patch[j].area/(*tarea);
		}
	}
	return(cmdf);
}

/*
   void RandE(TSymOp sym, vector d)
   ============================================================================
   Generates a random orientation as a Rodrigues vector in the fundamental zone.

   sym: symmetry operators
   d: Rodrigues vector orientation
   ============================================================================
*/
void RandE(TSymOp sym, vector d)
{
	double c;
	TAxAng aa;
	static long gljdum = 1;
	TEuler e;
	matrix g;

	e.a1 = 2*PI*ran3(&gljdum);
	e.a = acos2(2*ran3(&gljdum)-1);
	e.a2 = 2*PI*ran3(&gljdum);
	EToG(e,g);
	MinAngMatrix(sym,g);
	aa = GToAA(g);
	c = tan(aa.ang/2);
	d[0] = c*aa.n[0];
	d[1] = c*aa.n[1];
	d[2] = c*aa.n[2];
}

/*
   int metrop(double de, double t)
   ============================================================================
   "Oracle" that determines whether an orientation change or swap is successful 

   de: energy change
   t: temperature
   ============================================================================
*/
int metrop(double de, double t)
{
	static long gljdum = 1;

	return de < 0.0 || ran3(&gljdum) < exp(-de/t);
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

/*
   void DFSave(double *odf, char fname[100])
   ============================================================================
   Saves a distribution (ODF or MDF) to file

   *odf: distribution
   fname: file name
   ============================================================================
*/
void DFSave(double *odf, char fname[100])
{
	FILE *ifp;
	int i;

	ifp = fopen(fname,"w");
	for(i=0;i<NCELLS;i++){
		fprintf(ifp,"%0.4lf\n",odf[i]);
	}
	fclose(ifp);
}

/*
   TAxAng GToAA(matrix g)
   ============================================================================
   Converts orientation matrix to axis/angle pair.

   g: orientaiton matrix
   ============================================================================
*/
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

/*
   void Transpose(matrix g, matrix gt)
   ============================================================================
   Takes the transpose of matrix g.

   g: matrix
   gt: transposed matrix
   ============================================================================
*/
void Transpose(matrix g, matrix gt)
{
	int i,j;

	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			gt[i][j]=g[j][i];
}
 
/*
   void MM(matrix g1, matrix g2, matrix g3)
   ============================================================================
   Multiplies two 3x3 matrices.
   g1 x g2 = g3
   ============================================================================
*/    
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
 
/*
   void EToG(TEuler euler, matrix g)
   ============================================================================
   Converts Euler angles to orientation matrix
   
   euler: Euler angles
   g: orientation matrix
   ============================================================================
*/          
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
 
/*
   double acos2(double ca)
   ============================================================================
   arccos function with protection against precision errors
   ============================================================================
*/           
double acos2(double ca)
{
  if(ca < -1) ca = -1;
  if(ca > 1) ca = 1;
  return acos(ca);
}

/*
   TSymOp LoadSym(char fname[100])
   ============================================================================
   Loads symmetry operators from file; returns symmetry operator structure
   
   fname: file name
   ============================================================================
*/         
TSymOp LoadSym(char fname[100])
{
  FILE *ipf;
  TSymOp sym; 		
  int ii,i,j;

  ipf = fopen(fname,"r");
  if(!ipf){ printf("Error opening %s\n",fname); exit(1);}
  fscanf(ipf,"%d",&sym.num);
  sym.op = (matrix*)calloc(sym.num, sizeof(matrix));
  for(ii=0;ii<sym.num;ii++)
    for(i=0;i<3;i++)
      for(j=0;j<3;j++)
		fscanf(ipf,"%lf",&sym.op[ii][i][j]);
  fclose(ipf);
  return(sym);
} 

/*
   void MinAngMatrix(TSymOp sym, matrix g)
   ============================================================================
   Finds the symmetrically equivalent misorientation matrix that yields the
   minimum misorientation angle.
   
   sym: symmetry operators
   g: misorientation matrix
   ============================================================================
*/     
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

/*
   int OVectToCell(vector d)
   ============================================================================
   Determines the orientation cell index corresponding to HOMOCHORIC vector d; 
   returns cell index.
   
   d: homochoric orientatin vector
   ============================================================================
*/    
int OVectToCell(vector d)
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

/*
   int MVectToCell(vector d)
   ============================================================================
   Determines the misorientation cell index corresponding to HOMOCHORIC vector d; 
   returns cell index.
   
   d: homochoric misorientatin vector
   ============================================================================
*/    
int MVectToCell(vector d)
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

/*
   void VectToG(vector n, matrix g)
   ============================================================================
   Converts Rodrigues vector into orientation matrix
   
   n: Rodrigues vector
   g: orientation matrix
   ============================================================================
*/   
void VectToG(vector n, matrix g)
{
	double nmag,phi;
	vector nd;

	nmag = sqrt(SQR(n[0])+SQR(n[1])+SQR(n[2]));
	nd[0]= n[0]/nmag;
	nd[1]= n[1]/nmag;
	nd[2]= n[2]/nmag;
	phi = 2*atan(nmag);
	g[0][0] = cos(phi) + (1-cos(phi))*SQR(nd[0]);
    g[0][1] = -nd[2]*sin(phi)+(1-cos(phi))*nd[0]*nd[1];
    g[0][2] = +nd[1]*sin(phi)+(1-cos(phi))*nd[0]*nd[2];
    g[1][0] = +nd[2]*sin(phi)+(1-cos(phi))*nd[1]*nd[0];
    g[1][1] = cos(phi) + (1-cos(phi))*SQR(nd[1]);
    g[1][2] = -nd[0]*sin(phi)+(1-cos(phi))*nd[1]*nd[2];
    g[2][0] = -nd[1]*sin(phi)+(1-cos(phi))*nd[2]*nd[0];
    g[2][1] = +nd[0]*sin(phi)+(1-cos(phi))*nd[2]*nd[1];
    g[2][2] = cos(phi) + (1-cos(phi))*SQR(nd[2]);
}

/*
   void MisVect(TSymOp sym, matrix g1, matrix g2, vector d)
   ============================================================================
   Given orientation matrices g1 and g2, finds the misorientation in terms of 
   a Rodrigues vector in the fundamental zone.  This routine will only work for
   cubic symmetry.
   
   d: Rodrigues vector
   g: orientation matrix
   ============================================================================
*/
void MisVect(TSymOp sym, matrix g1, matrix g2, vector d)
{
	matrix g,gt;
	TAxAng aa;
	double c;

	Transpose(g2,gt);
	MM(g1,gt,g);
	MinAngMatrix(sym,g);
	aa = GToAA(g);
	c = tan(aa.ang/2);
	aa.n[0] = fabs(aa.n[0]);
	aa.n[1] = fabs(aa.n[1]);
	aa.n[2] = fabs(aa.n[2]);
	d[0] = MAXMAX(aa.n[0],aa.n[1],aa.n[2]);
	d[2] = MINMIN(aa.n[0],aa.n[1],aa.n[2]);
	d[1] = aa.n[0]+aa.n[1]+aa.n[2]-d[0]-d[2];
	d[0] = c*d[0];
	d[1] = c*d[1];
	d[2] = c*d[2];
}

/*
   void VConvert(vector r, vector d)
   ============================================================================
   Converts the Rodrigues vector representation to Homochoric vector.
   
   r: Rodrigues vector
   d: Homochoric vector
   ============================================================================
*/
void VConvert(vector r, vector d)
{
	double nmag,c,phi;

	nmag = sqrt(SQR(r[0])+SQR(r[1])+SQR(r[2]));
	phi = 2*atan(nmag);
	c = pow(0.75*(phi-sin(phi)),(double)1/3);
	d[0]= c*r[0]/nmag;
	d[1]= c*r[1]/nmag;
	d[2]= c*r[2]/nmag;
}


	
