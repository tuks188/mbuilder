#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define PI acos(-1)
#define EPS 10E-7 

static double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

typedef double vector[3];
typedef double matrix[3][3];
typedef struct {
	double a1,a,a2;
}TEuler;


void VectToG(vector n, matrix g);
TEuler GToE(matrix g);
double acos2(double ca);

int main()
{ 
	FILE *fp,*fp2;
	int id;
	vector d;
	matrix g;
	TEuler f;

	fp = fopen("orts.txt","r");
	fp2 = fopen("EAorts.txt","w");
	
	while(fscanf(fp,"%d%lf%lf%lf",&id,&d[0],&d[1],&d[2])!=EOF){
		VectToG(d,g);
		f = GToE(g);
		fprintf(fp2,"%d\t%0.4lf\t%0.4lf\t%0.4lf\n",id,f.a1,f.a,f.a2);
	}

	fclose(fp);
	fclose(fp2);

	return(0);
}

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

TEuler GToE(matrix g)
{
  double x,sf;
  TEuler f;
 
  x = g[2][2];
  if(x*x<1) 
    sf = sqrt(1-x*x);
  else 
    sf = 0.0;
  if(sf>EPS){
    f.a1 = -g[2][1]/sf;
    f.a1 = acos2(f.a1);
    if(g[2][0]<0) f.a1 = 2*PI-f.a1;
    f.a = acos2(x);
    f.a2 = g[1][2]/sf;
    f.a2 = acos2(f.a2);
    if(g[0][2]<0) f.a2 = 2*PI-f.a2;
  }
  else{
    f.a1 = g[1][1];
    f.a1 = acos2(f.a1);
    if(g[0][1]<0) f.a1 = 2*PI-f.a1;
    f.a = 0;
    f.a2 = 0;
  }
  return(f);
}

double acos2(double ca)
{
  if(ca < -1) ca = -1;
  if(ca > 1) ca = 1;
  return acos(ca);
}


