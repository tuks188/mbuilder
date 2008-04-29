// Modified 2/23/2005 by Steve Sintay to change the output
// to that of a *.wts file

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define PI acos(-1)
#define EPS 10E-7 
#define RAD_TO_DEG    57.29577951308

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
	FILE *fp,*fp2,*fp3;
	int id;
        int idmax=15000;
	vector d;
	matrix g;
	TEuler f;
        double vol;
        double odf[idmax];
        int sysstat;
        double totalvol=0.0;

	fp = fopen("orts.txt","r");
	fp2 = fopen("orts.wts","w");

        sysstat=system("awk ' { if( $1 ~ /<id>/ ) x=$2 } { if( $1 ~ /<nregions>/) xmax=$2 } { if( $1 ~ /<volume>/ ) arr[x]=$2 } END { for (i=0; i<=xmax ; i++) print i,arr[i] } ' cellIdealization.xml > tempfreq.txt ");
        if(sysstat!=0){
         printf("Problem executing UNIX utility AWK \n");
         exit(12);
        }

	fp3 = fopen("tempfreq.txt","r"); //  temp file(ID, Vol Fraction)

        if(!fp3 && !fp2 && !fp){
        printf("problem opening fp3 \n");
        }

// Read the evodf.txt file which contains VF of each orientation 
        while(fscanf(fp3,"%d%lf",&id,&vol)!=EOF){
          odf[id]=vol;
          if(vol>1.00000){totalvol+=vol;}
          //printf("%d   %0.8lf \n" ,id,odf[id]);
          if(id>idmax){
            printf("# of grain id's exceeds idmax \nIncrease the value of idmax \n");
            exit(1);
          } 
        }
 
   fprintf(fp2,"orts.wts file from MB fitted texture wts=volume\n");
   fprintf(fp2,"  Evm    F11    F12    F13    F21    F22    F23    F31    F32    F33\n");
   fprintf(fp2,"  0.000  1.000  0.000  0.000  0.000  1.000  0.000  0.000  0.000  1.000\n");
   fprintf(fp2,"Bunge:Psi  Theta   phi    weight");
   fprintf(fp2,"   (up to 6 state parameters, f8.2)   XYZ= 1 2 3\n");

   while(fscanf(fp,"%d%lf%lf%lf",&id,&d[0],&d[1],&d[2])!=EOF){
		VectToG(d,g);
		f = GToE(g);
                f.a1 *= RAD_TO_DEG;
                f.a *= RAD_TO_DEG;
                f.a2 *= RAD_TO_DEG;
                //printf("%d %0.8lf \n",id,odf[id]);
		//fprintf(fp2,"%d\t%0.4lf\t%0.4lf\t%0.4lf\n",id,f.a1,f.a,f.a2);
   fprintf(fp2,"%7.4lf  %7.4lf  %7.4lf  %0.8lf\n",f.a1,f.a,f.a2,odf[id]/totalvol);
	}

	fclose(fp);
	fclose(fp2);
	fclose(fp3);
        sysstat=system("rm tempfreq.txt");

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


