#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXACTIVE 10
#define MAXTRANSITIONS 100
#define ENGULFPACKET 64


/* structure definitions... */

struct xyzBuffer {
  double *abscissas; /* the point values in either x,y, or z... */
  int *indices; /* the indices associated with the sorted arrays of 
		      abscissas */
  int *reverseIndices; /* reverseIndices[i] is the index of i in indices... */
};

struct ellipsoid {
  double axisLengths[3]; /* the three axis lengths... */
  double axes[3][3]; /* the three axes, as unit vectors.
			axes[0][*] corresponds to axisLength[0],
			etc... */
};

struct occupant{
  int id; /* the id of the point contained... */
  double eFunc; /* the ellipsoid function of point id in the current 
		   ellipsoid... */
  double offset[3]; /* the offset required to bring point id into the current
		       ellipsoid... */
};

struct occupantList {
  struct occupant *theOccupant;
  struct occupantList *next;
};

struct hash{
  int basis; /* the hash function is index % basis... */
  int **hash; /* hash[i][0] is the index that hashes to the i'th location,
		 hash[i][1] is the location in the target array... */
};


struct foamCell {
  int id; /* the index leading to this particular foamCell... */
  int actives; /* the number of enclosing cells that are active... */
  int activeEngulfers[MAXACTIVE]; /* the indices of the enclosing cells... */
  double point[3]; /* the center of the ellipsoid... */
  struct ellipsoid e; /* the desired ellipsoidal cell... */
  int nContained; /* the number of points my ellipsoid contains... */
  struct occupant **containedIn; /* the points that my ellipse contains,
				    along with their ellipsoidFunctions... */
  struct hash *containedInHash; /* a hash to allow index containedIn lookup 
				   easily... */
  int nEngulfers; /* the number of foamCells whose ellipsoids contain this
		     point... */
  int *engulfers; /* the indices of the foamCells containing  this
		     foamCell... */
  struct hash  *engulferHash; /* a hash to allow engulfer index lookup 
				 easily... */
  double orientation[3]; /* the desired orientation, in bunge notation... */
  
};


struct bbox{
  double ranges[3][2]; /* the range in x,y,z...  
			  xmin = ranges[0][0], xmax = ranges[0][1], 
			  ymin = ranges[1][0], ymax = ranges[1][1], 
			  zmin = ranges[2][0], zmax = ranges[2][1]   */
};


/* a few convenient globals... */
long ranseed = 0l;
int transitionCount = 0;
double consumeAward = 1.0;
double overlapEncouragement = 1.0;
double zeroPenalty = 0.95;
double systemState = 0.0;
int maxTransitions;
double sumTransitions=0.0;
int currentOperation;
struct bbox *theBbox;
double paddingFraction=0.5;
int metaInformation=0;
int *immortals;
FILE *logFile;


/* function prototypes... */
double ellipsoidFunction(double *point, struct foamCell *theCell);
double addCost(int index, struct foamCell **theCells,int *activeList);
double subtractCost(int index, struct foamCell **theCells,int *activeList);
double swapCost(int markedForDeath,int reincarnated,
		struct foamCell **theCells,int *activeList);
void commitAdd(int index,struct foamCell **theCells,int *activeList);
void commitSubtract(int index,struct foamCell **theCells,int *activeList);
void commitSwap(int markedForDeath,int reincarnated,
		struct foamCell **theCells, int *activeList);
void updateActives(struct foamCell **theCells,int index,int *activeList);
double findEFunc(struct foamCell *me, int index, int *found);
void add(int nCells,struct foamCell **theCells,int *activeList);
void swap(int nCells,struct foamCell **theCells,int *activeList);
void jog(int nCells,struct foamCell **theCells,int *activeList);
void delete(int nCells,struct foamCell **theCells,int *activeList);
int oracle(double cost);
double annealingSchedule(double cost);
void genrandom(long *ranseed,int count,double *ranvalues);
void logState(double cost);
void reorderAxes(struct foamCell *theCell);
void sortBuffer(int count,struct xyzBuffer *buffer);
void findContained(int count,struct bbox *theBbox,
		   struct xyzBuffer *xBuffer,
		   struct xyzBuffer *yBuffer,
		   struct xyzBuffer *zBuffer,
		   int index,
		   struct foamCell **theCells);
int bisectFind(int count, double value, struct xyzBuffer *theBuffer);
int intervalCount(int nsegs,int *intervals);
int inInterval(int probe,int nsegs,int *intervals);
void findThisInterval(int nCells,
		      struct bbox *theBbox,
		      int theCoord,
		      struct xyzBuffer *theBuffer,
		      double vmin,
		      double vmax,
		      int *nsegs,int *breakPoints,
		      double *offsets);
void updateEngulfers(int id,struct foamCell **theCells);
void buildHashes(int id, struct foamCell **theCells);
void readPoints(char *fnam,int *nCells,struct foamCell ***theCells);
void repackEngulfers(struct foamCell *me);
void twiddle(int nCells,struct foamCell **theCells,int *activeList);
double overlap(double EFunc);
void outputCell(FILE *outfile,int noOffset,struct foamCell *theCell,
		struct foamCell **theCells,struct bbox *theBbox);
void genNewCells(struct foamCell **newCells,int *next,
		 int index,struct foamCell **oldCells,int *replacements);
void output(int nCells,struct foamCell **theCells,int *activeList,
	    struct bbox *theBbox);
void logState(double cost) {

  static int initialized=0;
  static FILE *logFile;
  static int lastCount;

  if (!initialized) {
    initialized = 1;
    logFile = fopen("annealing.log","w");
    lastCount = 0;
  }
  if ((int)(transitionCount/1000) > (int)(lastCount/1000)) {
    lastCount = transitionCount;
    fprintf(logFile,"%d %d %f %f %f %f\n",transitionCount,currentOperation,
	    sumTransitions,sumTransitions/transitionCount,cost,systemState);
    fflush(logFile);
  }
}




double overlap(double EFunc) {

  /* given the ellipsoid function of a point in two cells, calculate
     the penalty or benefit of the overlap.  the overlap function
     is -k*EFunc^2 + q, where k and q are determined by 
     zeroPenalty and overlapEncouragement... */

  static int initialized=0; /* have i calculated the parabolic overlap 
			       function... */
  static double k; 
  static double q;

  if (!initialized) {
    
    /* solve for k and q.  at EFunc == zeroPenalty the overlap function
       is 0.  at EFunc == 1 the overlap function is -overlapEncouragement. */

    initialized = 1;
    q = -overlapEncouragement/(1.0 - (1.0/(zeroPenalty*zeroPenalty)));
    k = q/(zeroPenalty*zeroPenalty);
  }
  return(-k*EFunc*EFunc + q);
}
			       

    

double ellipsoidFunction(double *point, struct foamCell *theCell) {

  /* returns the ellipsoid function (x/a)^2 + (y/b)^2 + (z/c)^2 of 
     point in the ellipsoidal cell theCell.  values of < 1 mean
     point is inside theCell, > 1 mean point is outside theCell, and
     exactly 1 means point lies on the surface of the ellipsoid... */

  int i,j; /* trusty indices... */
  double retval; /* (x/a)^2 + (y/b)^2 + (z/c)^2 */
  double offset[3]; /* point - center of ellipsoid... */
  double local[3]; /* offset in the local ellipsoid coordinate system... */

  for(i=0;i<3;i++) {
    offset[i] = point[i]-theCell->point[i];
  }
  retval = 0.0;
  for(i=0;i<3;i++) {
    local[i] = 0.0;
    for(j=0;j<3;j++) {
      local[i] += offset[j]*(theCell->e.axes[i][j]);
    }
    retval +=  (local[i]/(theCell->e.axisLengths[i]))*
      (local[i]/(theCell->e.axisLengths[i]));
  }
  return(retval);
}
  
double addCost(int index, struct foamCell **theCells,int *activeList) {

  /* returns the cost to make the index'th ellipsoid active... */

  int i,j; /* trusty indices... */
  double retval; /* the return value... */
  struct foamCell *me; /* shorthand for the index'th cell... */
  int neighbor; /* the contained point... */
  int theRival; /* another cell that contains the current point... */
  int found; /* have i found the ellipsoid function of the point in 
		question...*/
  
  me = theCells[index];
  retval = 0.0;

  for (i=0;i<me->nContained;i++) {
    neighbor = me->containedIn[i]->id;
    if (theCells[neighbor]->actives == 0) {
      retval -= consumeAward;
    }
    else {
      retval += overlap(me->containedIn[i]->eFunc);

      /* and include the collateral damage on the other side... */

      for (j=0;j<theCells[neighbor]->actives;j++) {
	theRival = theCells[neighbor]->activeEngulfers[j];
	retval += overlap(findEFunc(theCells[theRival],neighbor,&found));
	if (!found) {
	  fprintf(stderr,"Inconsistent structure detected\n");
	  fprintf(stderr,"foamCell:%d, neighbor:%d\n",theRival,neighbor);
	}
      }

    }
  }

  return(retval);

}

void commitAdd(int index,struct foamCell **theCells,int *activeList) {

  /* adds index to the activeList, and makes theCells reflect the change... */

  int i; /* trusty indices... */
  struct foamCell *me;
  int neighbor; /* a point contained in me... */

  me = theCells[index];

  activeList[index] = 1;
  for(i=0;i<me->nContained;i++) {
    neighbor= me->containedIn[i]->id;
    updateActives(theCells,neighbor,activeList);
  }  
}

double subtractCost(int index, struct foamCell **theCells,
		    int *activeList) {

  /* returns the cost to make the index'th ellipsoid inActive... */

  int i,j; /* trusty indices... */
  double retval; /* the return value... */
  struct foamCell *me; /* shorthand for the index'th cell... */
  int neighbor; /* the contained point... */
  int theRival; /* another active cell that contains neighbor... */
  int found; /* have i found the ellipsoid function for neighbor? */

  if (!activeList[index]) {
    fprintf(stderr,"subtractCost called on inactive cell...\n");
    fprintf(stderr,"index:%d\n",index);
  }
  
  me = theCells[index];
  retval = 0.0;

  for (i=0;i<me->nContained;i++) {
    neighbor = me->containedIn[i]->id;
    if (theCells[neighbor]->actives == 1) {
      if (theCells[neighbor]->activeEngulfers[0] == index) {
	retval += consumeAward;
      }
      else {
	fprintf(stderr,"%d thinks it is active and contains %d\n",
		index,neighbor);
	fprintf(stderr,"but %d thinks it is contained by %d\n",
		neighbor,theCells[neighbor]->activeEngulfers[0]);
      }
    }
    else {
      if (theCells[neighbor]->actives > 1) {
	retval -= overlap(me->containedIn[i]->eFunc);

	/* and include the collateral damage on the other side... */

	for (j=0;j<theCells[neighbor]->actives;j++) {
	  theRival = theCells[neighbor]->activeEngulfers[j];
	  if (theRival != index) {
	    retval -= overlap(findEFunc(theCells[theRival],neighbor,&found));
	    if (!found) {
	      fprintf(stderr,"Inconsistent structure detected\n");
	      fprintf(stderr,"foamCell:%d, index:%d\n",theRival,index);
	    }
	  }
	}
      }
      else {
	fprintf(stderr,"Inconsistent structure detected\n");
	fprintf(stderr,"neighbor:%d doesn't know he's in active cell:%d\n",
		neighbor,index);
      }
    }
  }

  return(retval);

}


void commitSubtract(int index,struct foamCell **theCells,int *activeList) {

  /* removes index from the activeList, and makes theCells reflect the 
     change... */

  int i; /* trusty index... */
  struct foamCell *me;
  int neighbor; /* a point contained in me... */

  me = theCells[index];

  activeList[index] = 0;
  for(i=0;i<me->nContained;i++) {
    neighbor= me->containedIn[i]->id;
    updateActives(theCells,neighbor,activeList);
  }  
}


void updateActives(struct foamCell **theCells,int index,int *activeList) {
  
  int i; /* trusty index... */
  int next; /* the next index in activeEngulfers to fill... */
  struct foamCell *me;
  
  me = theCells[index];

  for(i=0;i<MAXACTIVE;i++) {
    me->activeEngulfers[i] = -1;
  }

  next = 0;
  me->actives = 0;
  for(i=0;i<me->nEngulfers;i++) {
    if (activeList[me->engulfers[i]]) {
      me->actives++;
      me->activeEngulfers[next++] = me->engulfers[i];
    }
  }
}


double findEFunc(struct foamCell *me, int index, int *found) {

  int probe; /* the place to start looking in the containedInHash... */
  int i; /* trusty index... */
  int maxVal; /* the size of the containedInHash... */

  maxVal = me->containedInHash->basis;
  probe = index % (me->containedInHash->basis);
  for(i=probe;i<maxVal;i++) {
    if (me->containedInHash->hash[i][0] == index) {
      *found =me->containedInHash->hash[i][1] + 1;
      return(me->containedIn[me->containedInHash->hash[i][1]]->eFunc);
    }
    else {
      if (me->containedInHash->hash[i][0] < 0) {
	*found = 0;
	return(0.0);
      }
    }
  }

  /* could get here from a collision at the last element... */
  for(i=0;i<probe;i++) {
    if (me->containedInHash->hash[i][0] == index) {
      *found = me->containedInHash->hash[i][1] + 1;
      return(me->containedIn[me->containedInHash->hash[i][1]]->eFunc);
    }
    else {
      if (me->containedInHash->hash[i][0] < 0) {
	*found = 0;
	return(0.0);
      }
    }
  }
  *found = 0;
  return(0.0);
}


double swapCost(int markedForDeath,int reincarnated,
		struct foamCell **theCells,int *activeList) {


  /* the cost of killing marked for death and replacing him with 
     reincarnated.  be careful to correct for the effect of markedForDeath
     in the returned cost... */

  int i,j; /* trusty indices... */
  double retval; /* the return value... */
  struct foamCell *me; /* shorthand for the current cell... */
  int neighbor; /* the contained point... */
  int theRival; /* another active cell that contains neighbor... */
  int found; /* have i found the ellipsoid function for neighbor? */
  
  if ((activeList[markedForDeath]) && (!activeList[reincarnated])) {
    
  
    me = theCells[reincarnated];
    retval = 0.0;

    for (i=0;i<me->nContained;i++) {
      neighbor = me->containedIn[i]->id;
      if (theCells[neighbor]->actives == 0) {
	retval -= consumeAward;
      }
      else {
	if ((theCells[neighbor]->actives == 1) &&
	    (theCells[neighbor]->activeEngulfers[0] == markedForDeath)) {
	  retval -= consumeAward;
	}
	else {
	  retval += overlap(me->containedIn[i]->eFunc);
	}

	/* and include the collateral damage on the other side... */

	for (j=0;j<theCells[neighbor]->actives;j++) {
	  theRival = theCells[neighbor]->activeEngulfers[j];
	  if (theRival != markedForDeath) {
	    retval += 
	      overlap(findEFunc(theCells[theRival],neighbor,&found));
	    if (!found) {
	      fprintf(stderr,"Inconsistent structure detected\n");
	      fprintf(stderr,"foamCell:%d, neighbor:%d\n",theRival,neighbor);
	    }
	  }
	}
      }
    }
    me = theCells[markedForDeath];

    for (i=0;i<me->nContained;i++) {
      neighbor = me->containedIn[i]->id;
      if (theCells[neighbor]->actives == 1) {
	if (theCells[neighbor]->activeEngulfers[0] == markedForDeath) {
	  retval += consumeAward;
	}
	else {
	  fprintf(stderr,"%d thinks it is active and contains %d\n",
		  markedForDeath,neighbor);
	  fprintf(stderr,"but %d thinks it is contained by %d\n",
		  neighbor,theCells[neighbor]->activeEngulfers[0]);
	}
      }
      else {
	if (theCells[neighbor]->actives > 1) {
	  retval -= overlap(me->containedIn[i]->eFunc);

	  /* and include the collateral damage on the other side... */

	  for (j=0;j<theCells[neighbor]->actives;j++) {
	    theRival = theCells[neighbor]->activeEngulfers[j];
	    if (theRival != markedForDeath) {
	      retval -= 
		overlap(findEFunc(theCells[theRival],neighbor,&found));
	      if (!found) {
		fprintf(stderr,"Inconsistent structure detected\n");
		fprintf(stderr,"foamCell:%d, neighbor:%d\n",theRival,neighbor);
	      }
	    }
	  }
	}
	else {
	  fprintf(stderr,"Inconsistent structure detected\n");
	  fprintf(stderr,"neighbor:%d doesn't know he's in active cell:%d\n",
		  neighbor,markedForDeath);
	}
      }
    }

    return(retval);

  }
  else {
    fprintf(stderr,"Inconsistent call to swapCost:\n");
    fprintf(stderr,"markedForDeath:%d, activeList[markedForDeath]:%d\n",
	    markedForDeath,activeList[markedForDeath]);
    fprintf(stderr,"reincarnated:%d, activeList[reincarnated]:%d\n",
	    reincarnated,activeList[reincarnated]);
    return(HUGE);
  }

}



void add(int nCells,struct foamCell **theCells,
	 int *activeList) {

  /* selects an unenclosed cell at random, evaluates the cost of adding it,
     adds it if the oracle says yes, and updates the systemState... */

  int next; /* the next foamCell to check... */
  int theIndex; /* the place to start looking for a cell to turn on... */
  double cost; /* the change in systemState from adding this cell... */
  double deviate; /* a random deviate... */
  int found; /* have i found a candidate to add... */
  int priorIndex; /* the predecessor to the starting position... */
  
  found = 0;
  genrandom(&ranseed,1,&deviate);
  
  theIndex = (int)(nCells*deviate);
  if (theIndex != 0) {
    priorIndex = theIndex-1;
  }
  else {
    priorIndex = nCells - 1;
  }
  
  next = theIndex;
  while ((!found) && (next != priorIndex)) {
    if (theCells[next]->actives == 0) {
      found = 1;
      cost = addCost(next,theCells,activeList);
      if (oracle(cost)) {
	commitAdd(next,theCells,activeList);
	systemState += cost;
	logState(cost);
	return;
      }
    }
    next++;
    if (next == nCells) {
      next = 0;
    }
  }
}


void swap(int nCells,struct foamCell **theCells,int *activeList) {

  /* selects an active and an inactive cell at random, evaluates the cost of 
     swapping their state, commits if the oracle says yes, 
     and updates the systemState... */

  int next; /* the next foamCell to check... */
  int theIndex; /* the place to start looking for a cell to turn on... */
  double cost; /* the change in systemState from adding this cell... */
  double deviate; /* a random deviate... */
  int found; /* have i found a candidate to add... */
  int markedForDeath; /* the cell to turn off... */
  int reincarnated; /* the cell to turn on... */
  int priorIndex; /* the predecessor to the starting position... */
  
  found = 0;
  genrandom(&ranseed,1,&deviate);
  
  theIndex = (int)(nCells*deviate);
  if (theIndex != 0) {
    priorIndex = theIndex-1;
  }
  else {
    priorIndex = nCells - 1;
  }
  
  next = theIndex;
  while ((!found) && (next != priorIndex)) {
    if ((activeList[next]) && (!immortals[next])) {
      found = 1;
      markedForDeath = next;
    }
    next++;
    if (next == nCells) {
      next = 0;
    }
  }
  if (found) {
    found = 0;
    genrandom(&ranseed,1,&deviate);
  
    theIndex = (int)(nCells*deviate);
    if (theIndex != 0) {
      priorIndex = theIndex-1;
    }
    else {
      priorIndex = nCells - 1;
    }
    
    next = theIndex;
    while ((!found) && (next != priorIndex)) {
      if (!activeList[next]) {
	found = 1;
	reincarnated = next;
      }
      next++;
      if (next == nCells) {
	next = 0;
      }
    }
    if (found) {
      cost = swapCost(markedForDeath,reincarnated,theCells,activeList);
      if (oracle(cost)) {
	commitSwap(markedForDeath,reincarnated,theCells,activeList);
	systemState += cost;
	logState(cost);
      }
    }
  }
}



void jog(int nCells,struct foamCell **theCells,int *activeList) {

  /* selects an active cell at random, evaluates the cost of 
     swapping the cell with one contained in it, commits if the oracle says 
     yes, and updates the systemState... */

  int next; /* the next foamCell to check... */
  int theIndex; /* the place to start looking for a cell to turn on... */
  double cost; /* the change in systemState from adding this cell... */
  double deviate; /* a random deviate... */
  int found; /* have i found a candidate to add... */
  int markedForDeath; /* the cell to turn off... */
  int reincarnated; /* the cell to turn on... */
  int priorIndex; /* the predecessor to the starting position... */
  
  found = 0;
  genrandom(&ranseed,1,&deviate);
  
  theIndex = (int)(nCells*deviate);
  if (theIndex != 0) {
    priorIndex = theIndex-1;
  }
  else {
    priorIndex = nCells - 1;
  }
  
  next = theIndex;
  while ((!found) && (next != priorIndex)) {
    if ((activeList[next]) && (!immortals[next])) {
      found = 1;
      markedForDeath = next;
    }
    next++;
    if (next == nCells) {
      next = 0;
    }
  }
  if (found) {

    found = 0;
    genrandom(&ranseed,1,&deviate);
  
    theIndex = (int)((theCells[markedForDeath]->nContained)*deviate);
    if (theIndex != 0) {
      priorIndex = theIndex-1;
    }
    else {
      priorIndex = (theCells[markedForDeath]->nContained) - 1;
    }
    
    next = theIndex;
    while ((!found) && (next != priorIndex)) {
      reincarnated = theCells[markedForDeath]->containedIn[next]->id;
      if ((theCells[reincarnated]->actives == 1) && 
	  (theCells[reincarnated]->activeEngulfers[0] == markedForDeath) &&
	  (reincarnated != markedForDeath)) {
	found = 1;
      }
      next++;
      if (next == theCells[markedForDeath]->nContained) {
	next = 0;
      }
    }
    if (found) {
      cost = swapCost(markedForDeath,reincarnated,theCells,activeList);
      if (oracle(cost)) {
	commitSwap(markedForDeath,reincarnated,theCells,activeList);
	systemState += cost;
	logState(cost);
      }
    }
  }
}

void delete(int nCells,struct foamCell **theCells,int *activeList) {

  /* chooses an active cell at random, calculates the cost of deleting it,
     commits the delete if the oracle says yes, and updates the 
     systemState... */

  int next; /* the next foamCell to check... */
  int theIndex; /* the place to start looking for a cell to turn on... */
  double cost; /* the change in systemState from adding this cell... */
  double deviate; /* a random deviate... */
  int found; /* have i found a candidate to add... */
  int priorIndex; /* the predecessor to the starting position... */
  
  found = 0;
  genrandom(&ranseed,1,&deviate);
  
  theIndex = (int)(nCells*deviate);
  if (theIndex != 0) {
    priorIndex = theIndex-1;
  }
  else {
    priorIndex = nCells - 1;
  }
  
  next = theIndex;
  while ((!found) && (next != priorIndex)) {
    if ((activeList[next]) && (!immortals[next])) {
      found = 1;
      cost = subtractCost(next,theCells,activeList);
      if (oracle(cost)) {
	commitSubtract(next,theCells,activeList);
	systemState += cost;
	logState(cost);
	return;
      }
    }
    next++;
    if (next == nCells) {
      next = 0;
    }
  }
}


      
void commitSwap(int markedForDeath,int reincarnated,
		struct foamCell **theCells, int *activeList) {

  /* updates theCells to reflect the deletion of markedForDeath from 
     activeList, and the addition of reincarnated... */

  commitSubtract(markedForDeath,theCells,activeList);
  commitAdd(reincarnated,theCells,activeList);

}



int oracle(double cost) {

  double prob; /* the acceptance probability... */
  double deviate; /* a probe... */


  sumTransitions += fabs(cost);
  transitionCount++;
  if (cost < 0.0) {
    return(1);
  }
  else {
    prob = exp(-cost/annealingSchedule(cost));
  }
  
  genrandom(&ranseed,1,&deviate);

  if (deviate < prob) {
    return(1);
  }
  else {
    return(0);
  }
}


double annealingSchedule(double cost) {
  
  /* calculates an annealing schedule, which is a specification of the
     how likely the program is to accept a positive transition (from a better
     to a worse state).  the annealing schedule is specified in terms of 
     a fraction of the average transition magnitude, as a function of 
     transitionCount.  for instance, if there is no annealing.schedule file
     the default annealing schedule is 
     
     0 0.4
     100000 0.0001

     which means that to begin with, positive transitions up to 0.4 times the 
     average transition magnitude have a 50% chance of being accepted, but by 
     100000 transitions, positive transitions of only 0.0001*the average 
     transition magnitude have a 50% chance of being accepted... */

  static int nbreaks; /* the number of break points in the annealing 
			 schedule... */
  static int xs[MAXTRANSITIONS]; /* the x values of the annealing schedule...*/
  static double ys[MAXTRANSITIONS]; /* the y values of the annealing 
				       schedule... */
  static int initialized=0; /* have i initialized xs and ys... */
  double frac; /* the fraction of the average transition at which the
		  probability of acceptance is 50%... */
  double avg; /* the average positive transition... */
  int i; /* trusty index... */
  int theBreak; /* the interval of the annealing schedule in which 
		   we currently fall... */
  FILE *annealingFile; /* potentially a file with an annealing schedule */
  
  if (!initialized) {
    initialized = 1;
    annealingFile = fopen("annealing.schedule","r");
    if (annealingFile == NULL) {
      nbreaks = 2;
      xs[0] = 0;
      xs[1] = 100000;
      ys[0] = 0.4;
      ys[1] = 0.0001;
    }
    else {
      fscanf(annealingFile,"%d",&nbreaks);
      for(i=0;i<nbreaks;i++){
	fscanf(annealingFile,"%d %lf",&(xs[i]),&(ys[i]));
      }
    }
    maxTransitions = xs[nbreaks-1];
  }
  if (transitionCount > 0) {
    avg = sumTransitions/transitionCount;
  }
  else {
    avg = sumTransitions;
  }
  theBreak = 0;
  if (transitionCount <= xs[0]) {
    frac = ys[0];
  }
  else {
    if (transitionCount >= xs[nbreaks-1]) {
      frac = ys[nbreaks-1];
    }
    else {
      for(i=0;i<nbreaks-1;i++) {
	if ((xs[i] <= transitionCount)  && (xs[i+1] > transitionCount)) {
	  frac = ys[i] + ((double)(transitionCount - xs[i])/
			  (double)(xs[i+1]-xs[i]))*(ys[i+1]-ys[i]);
	  break;
	}
      }
    }
  }
  
  /* 1/2 = exp(-frac*avg/retval) => retval = (frac*avg)/0.69 */
  
  return(frac*avg/0.69);
}


#define STOREDSTATE 97

void genrandom(long *ranseed,int count,double *ranvalues)

{

  /* this function generates count random numbers between 0 and 1, shuffled
     to avoid serial correlation.  the generator keeps a cache of random 
     numbers, which it accesses and replenishes in a random order.  the cache
     is initialized at first call.  if *ranseed is zero, reinitialize the 
     random number cache, and the ranseed, which is initialized from the 
     current time... */

  static int initialized=0; /* has the history vector been initialized... */
  static double history[STOREDSTATE]; /* a bunch of random doubles to choose
					 from... */
  int i; /* an index */
  double probe; /* random deviate used for calculating selected value... */
  int selection; /* the index into history of the selected item... */
  
  time_t *tloc; /* the current time value, if needed */


  if ((!initialized) || (ranseed == 0)) {
    initialized = 1;

    if (*ranseed == 0) {

      tloc = (time_t *)malloc(sizeof(time_t));
      *ranseed = (long)time(tloc);
    }

    /* initialize the random number generator... */

    srand48(*ranseed); 
    /* srand(*ranseed); */
    
    /* fill the history vector... */

    for (i=0;i<STOREDSTATE;i++) 
       history[i] = drand48(); 
    /* history[i] = ((double)rand())/((double)RAND_MAX+1); */

  }

  for (i=0;i<count;i++) {
    probe = drand48();
    /* probe = ((double)rand())/((double)RAND_MAX+1); */
    selection = (int)(STOREDSTATE*drand48()); 
    /* selection = (int)(STOREDSTATE*((double)rand())/((double)RAND_MAX+1));*/
    ranvalues[i] = history[selection];
    history[selection] = drand48(); 
    /* history[selection] = ((double)rand())/((double)RAND_MAX+1);*/
  }

}


void readPoints(char *fnam,int *nCells,struct foamCell ***theCells) {
  
  /* read the points, along with their associated ellipsoids and orientations.
     allocate foamCell, fill in the ellipsoid and the orientation.

     the format of the input file is:

     foreach point:

        point (3 doubles),
        ellipsoid axis lengths (3 doubles),
	ellipsoid axes (9 doubles, unit vector for 
	                three axes),
	orientation (3 doubles)

  */

  int count; /* the current count of points... */
  double record[18]; /* one point record... */
  int i,j; /* trusty indices... */
  int fullrecord; /* did we have a partial record at the end of the file? */
  FILE *infile; /* the input file pointer... */
  struct xyzBuffer *xBuffer; /* the sorted x coordinates of the centers... */
  struct xyzBuffer *yBuffer; /* the sorted y coordinates of the centers... */
  struct xyzBuffer *zBuffer; /* the sorted z coordinates of the centers... */

  infile = fopen(fnam,"r");
  if (infile == NULL) {
    fprintf(stderr,"Input file:%s cannot be read.  Aborting...\n",
	    fnam);
    exit(1);
  }
  count = 0;
  while(fscanf(infile,"%lf",&(record[0])) == 1) {
    fullrecord = 1;
    for(i=1;i<18;i++) {
      if (fscanf(infile,"%lf",&(record[i])) != 1) {
	fullrecord = 0;
      }
    }
    if (fullrecord) {
      count++;
    }
  }
  fclose(infile);

  infile = fopen(fnam,"r");
  *nCells = count;
  *theCells = (struct foamCell **)malloc(count*sizeof(struct foamCell *));
  if (*theCells == NULL) {
    fprintf(stderr,"Malloc failure...\n");
  }
  
  for(i=0;i<count;i++) {
    for(j=0;j<18;j++) {
      fscanf(infile,"%lf",&(record[j]));
    }

    (*theCells)[i] = (struct foamCell *)malloc(sizeof(struct foamCell));
    if ((*theCells)[i] == NULL) {
      fprintf(stderr,"Malloc failure...\n");
    }

    for(j=0;j<3;j++) {
      (*theCells)[i]->point[j] = record[j];
      (*theCells)[i]->e.axisLengths[j] = record[j+3];
      (*theCells)[i]->e.axes[0][j] = record[j+6];
      (*theCells)[i]->e.axes[1][j] = record[j+9];
      (*theCells)[i]->e.axes[2][j] = record[j+12];
      (*theCells)[i]->orientation[j] = record[j+15];
      (*theCells)[i]->id = i;

      /* initialize stuff to be populated later... */
      (*theCells)[i]->nContained = 0;
      (*theCells)[i]->actives = 0;
      (*theCells)[i]->nEngulfers = 0;
      (*theCells)[i]->engulfers = NULL;
    }
    reorderAxes((*theCells)[i]);
  }

  /* build xBuffer, yBuffer, and zBuffer to hold sorted x, y, and z 
     buffers... */

  xBuffer = (struct xyzBuffer *)malloc(sizeof(struct xyzBuffer));
  xBuffer->abscissas = (double *)malloc(count*sizeof(double));
  xBuffer->indices = (int *)malloc(count*sizeof(int));
  xBuffer->reverseIndices = (int *)malloc(count*sizeof(int));
  yBuffer = (struct xyzBuffer *)malloc(sizeof(struct xyzBuffer));
  yBuffer->abscissas = (double *)malloc(count*sizeof(double));
  yBuffer->indices = (int *)malloc(count*sizeof(int));
  yBuffer->reverseIndices = (int *)malloc(count*sizeof(int));
  zBuffer = (struct xyzBuffer *)malloc(sizeof(struct xyzBuffer));
  zBuffer->abscissas = (double *)malloc(count*sizeof(double));
  zBuffer->indices = (int *)malloc(count*sizeof(int));
  zBuffer->reverseIndices = (int *)malloc(count*sizeof(int));

  for(i=0;i<count;i++) {
    xBuffer->indices[i] = i;
    xBuffer->abscissas[i] = (*theCells)[i]->point[0];
    yBuffer->indices[i] = i;
    yBuffer->abscissas[i] = (*theCells)[i]->point[1];
    zBuffer->indices[i] = i;
    zBuffer->abscissas[i] = (*theCells)[i]->point[2];
  }

  /* ok, now sort the buffers... */

  sortBuffer(count,xBuffer);
  sortBuffer(count,yBuffer);
  sortBuffer(count,zBuffer);

  /* and build the reverseIndices... */
  for(i=0;i<count;i++) {
    xBuffer->reverseIndices[xBuffer->indices[i]] = i;
    yBuffer->reverseIndices[yBuffer->indices[i]] = i;
    zBuffer->reverseIndices[zBuffer->indices[i]] = i;
  }

  /* now find the bounding box... */

  if (!metaInformation) {
    theBbox = (struct bbox *)malloc(sizeof(struct bbox));
    
    theBbox->ranges[0][0] = xBuffer->abscissas[0];
    theBbox->ranges[0][1] = xBuffer->abscissas[count-1];
    theBbox->ranges[1][0] = yBuffer->abscissas[0];
    theBbox->ranges[1][1] = yBuffer->abscissas[count-1];
    theBbox->ranges[2][0] = zBuffer->abscissas[0];
    theBbox->ranges[2][1] = zBuffer->abscissas[count-1];
  }
  else {
    if (xBuffer->abscissas[0] < theBbox->ranges[0][0]) {
      fprintf(stderr,"Requested xmin of %f overridden to %f...\n",
	      theBbox->ranges[0][0],xBuffer->abscissas[0]);
      theBbox->ranges[0][0] = xBuffer->abscissas[0];
    }
    if (xBuffer->abscissas[count-1] > theBbox->ranges[0][1]) {
      fprintf(stderr,"Requested xmax of %f overridden to %f...\n",
	      theBbox->ranges[0][1],xBuffer->abscissas[count-1]);
      theBbox->ranges[0][1] = xBuffer->abscissas[count-1];
    }

    if (yBuffer->abscissas[0] < theBbox->ranges[1][0]) {
      fprintf(stderr,"Requested ymin of %f overridden to %f...\n",
	      theBbox->ranges[1][0],yBuffer->abscissas[0]);
      theBbox->ranges[1][0] = yBuffer->abscissas[0];
    }
    if (yBuffer->abscissas[count-1] > theBbox->ranges[1][1]) {
      fprintf(stderr,"Requested ymax of %f overridden to %f...\n",
	      theBbox->ranges[1][1],yBuffer->abscissas[count-1]);
      theBbox->ranges[1][1] = yBuffer->abscissas[count-1];
    }

    if (zBuffer->abscissas[0] < theBbox->ranges[2][0]) {
      fprintf(stderr,"Requested zmin of %f overridden to %f...\n",
	      theBbox->ranges[2][0],zBuffer->abscissas[0]);
      theBbox->ranges[2][0] = zBuffer->abscissas[0];
    }
    if (zBuffer->abscissas[count-1] > theBbox->ranges[2][1]) {
      fprintf(stderr,"Requested zmax of %f overridden to %f...\n",
	      theBbox->ranges[2][1],zBuffer->abscissas[count-1]);
      theBbox->ranges[2][1] = zBuffer->abscissas[count-1];
    }
  }


  

  for(i=0;i<count;i++) {

    /* now find out who's inside of whom... */
    findContained(count,theBbox,xBuffer,yBuffer,zBuffer,i,*theCells);

    /* and make sure that everybody knows who they're inside... */
    updateEngulfers(i,*theCells);

    if (i == 0) {
      fprintf(stderr, "Building data structure #:\n");
      fprintf(stderr, "Number of Ellipsoids: %d\n",count);
    }
    if (i<10) {
      fprintf(stderr, "%d\n",i);
    }
    else {
      if (i<100) {
	if ((i % 10) == 0) {
	  fprintf(stderr,"%d\n",i);
	}
      }
      else {
	if (i < 1000) {
	  if ((i % 100) == 0) {
	    fprintf(stderr,"%d\n",i);
	  }
	}
	else {
	  if ((i % 1000) == 0) {
	    fprintf(stderr,"%d\n",i);
	  }
	}
      }
    }
  }

  for(i=0;i<count;i++) {
    buildHashes(i,*theCells);
  }

}
	

void buildHashes(int id, struct foamCell **theCells) {
  
  /* builds hashes of containedIn and engulfers.  just use linear insertion
     for collisions... */

  struct foamCell *me; /* shorthand for theCells[id]... */
  struct hash *thisHash; /* the current hash... */
  int i; /* trusty index... */
  int hashIndex; /* a hash index for the current id... */
  int firstIndex; /* the first hash index... */
  int found; /* have i found a home for this item yet? */

  me = theCells[id];

  thisHash = (struct hash *)malloc(sizeof(struct hash));
  thisHash->basis = 2*me->nContained - 1;
  thisHash->hash = (int **)malloc((thisHash->basis)*sizeof(int *));
  for(i=0;i<thisHash->basis;i++) {
    thisHash->hash[i] = (int *)malloc(2*sizeof(int));
    thisHash->hash[i][0] = -1;
    thisHash->hash[i][1] = -1;
  }
  for(i=0;i<me->nContained;i++) {
    hashIndex = me->containedIn[i]->id % thisHash->basis;
    firstIndex = hashIndex;
    if (thisHash->hash[hashIndex][0] < 0) {
      thisHash->hash[hashIndex][0] = me->containedIn[i]->id;
      thisHash->hash[hashIndex][1] = i;
    }
    else {
      found = 0;
      hashIndex++;
      while((!found) && (hashIndex < thisHash->basis)) {
	if (thisHash->hash[hashIndex][0] < 0) {
	  found = 1;
	  thisHash->hash[hashIndex][0] = me->containedIn[i]->id;
	  thisHash->hash[hashIndex][1] = i;
	}
	else {
	  hashIndex++;
	}
      }
      if (!found) {
	hashIndex = 0;
	while((!found) && (hashIndex < firstIndex)) {
	  if (thisHash->hash[hashIndex][0] < 0) {
	    found = 1;
	    thisHash->hash[hashIndex][0] = me->containedIn[i]->id;
	    thisHash->hash[hashIndex][1] = i;
	  }
	  else {
	    hashIndex++;
	  }
	}
	if (!found) {
	  fprintf(stderr,
		  "Violation of the pigeonhole principle in buildHash...\n");
	  fprintf(stderr,"Aborting... \n");
	  exit(4);
	}
      }
    }
  }
  me->containedInHash = thisHash;

  thisHash = (struct hash *)malloc(sizeof(struct hash));
  thisHash->basis = 4*me->nEngulfers - 1;
  thisHash->hash = (int **)malloc((thisHash->basis)*sizeof(int *));
  for(i=0;i<thisHash->basis;i++) {
    thisHash->hash[i] = (int *)malloc(2*sizeof(int));
    thisHash->hash[i][0] = -1;
    thisHash->hash[i][1] = -1;
  }
  for(i=0;i<me->nEngulfers;i++) {
    hashIndex = me->engulfers[i] % thisHash->basis;
    firstIndex = hashIndex;
    if (thisHash->hash[hashIndex][0] < 0) {
      thisHash->hash[hashIndex][0] = me->engulfers[i];
      thisHash->hash[hashIndex][1] = i;
    }
    else {
      found = 0;
      hashIndex++;
      while((!found) && (hashIndex < thisHash->basis)) {
	if (thisHash->hash[hashIndex][0] < 0) {
	  found = 1;
	  thisHash->hash[hashIndex][0] = me->engulfers[i];
	  thisHash->hash[hashIndex][1] = i;
	}
	else {
	  hashIndex++;
	}
      }
      if (!found) {
	hashIndex = 0;
	while((!found) && (hashIndex < firstIndex)) {
	  if (thisHash->hash[hashIndex][0] < 0) {
	    found = 1;
	    thisHash->hash[hashIndex][0] = me->engulfers[i];
	    thisHash->hash[hashIndex][1] = i;
	  }
	  else {
	    hashIndex++;
	  }
	}
	if (!found) {
	  fprintf(stderr,
		  "Violation of the pigeonhole principle in buildHash...\n");
	  fprintf(stderr,"Aborting... \n");
	  exit(4);
	}
      }
    }
  }
  me->engulferHash = thisHash;
}







void updateEngulfers(int id,struct foamCell **theCells) {

  /* makes sure that all the cells inside cell id know that they're inside
     cell id... */

  
  struct foamCell *me; /* shorthand for theCells[id]... */
  int engulfee; /* the id of a point contained in me... */
  int i; /* trusty index... */

  me = theCells[id];

  for(i=0;i<me->nContained;i++) {
    engulfee = me->containedIn[i]->id;
    if ((theCells[engulfee]->nEngulfers % ENGULFPACKET) == 0) {
      
      /* increase the length of the engulfers list by ENGULFPACKET... */
      repackEngulfers(theCells[engulfee]);
    }
    theCells[engulfee]->engulfers[theCells[engulfee]->nEngulfers++] = id;
  }
}

void repackEngulfers(struct foamCell *me) {

  /* malloc a new list of engulfers, with ENGULFPACKET more storage
     locations, copy all the old engulfers into place, and 
     free the old engulfers memory.  Assumes that the current array
     is full... */

  int currentSize; /* the current size of the engulfers array... */
  int nextSize; /* the new size of the engulfers array... */
  int *newArray; /* the new array of engulfers... */
  int i; /* trusty index... */

  currentSize = (me->nEngulfers);
  nextSize = currentSize+ENGULFPACKET;
  newArray = (int *)malloc(nextSize*sizeof(int));

  for(i=0;i<me->nEngulfers;i++) {
    newArray[i] = me->engulfers[i];
  }
  for(i=me->nEngulfers;i<nextSize;i++) {
    newArray[i] = -1;
  }

  free(me->engulfers);
  me->engulfers = newArray;

}
  

void findContained(int count,struct bbox *theBbox,
		   struct xyzBuffer *xBuffer,
		   struct xyzBuffer *yBuffer,
		   struct xyzBuffer *zBuffer,
		   int index,
		   struct foamCell **theCells) {

  /* updates the containedIn data structure of theCells[index].
     takes into account the periodicity of space, using the bounding box
     theBbox... */

  struct foamCell *me; /* shorthand for theCells[index]... */
  int nxr; /* number of intervals in x (can be 1 or 2) */
  int nyr; /* ditto for y... */
  int nzr; /* ditto for z... */
  int ixr[4]; /* intervals in xBuffer->indices that are close enough to be
		 considered for inclusion in me... */
  int iyr[4]; /* ditto for yBuffer->indices... */
  int izr[4]; /* ditto for zBuffer->indices... */
  double xOffset[4]; /* offset for each of the indices in ixr... */
  double yOffset[4]; /* offset for each of the indices in iyr... */
  double zOffset[4]; /* offset for each of the indices in izr... */

  double theEFunc; /* the ellipsoid function of the current candidate... */
  struct occupantList *theOccupants; /* a linked list from which i will build
					an array... */
  struct occupantList *head; /* something to hold onto, while inserting
				from the front... */
  int nOccupants; /* the current number of occupants in the linked list... */
  struct occupant *theOccupant; /* the current occupant... */
  int smallestInterval; /* do i want to traverse the interval in xBuffer,
			   yBuffer, or zBuffer... */
  int i,j; /* trusty indices... */
  int nchunk; /* the block of indices could conceivably be split into two
		 chunks.  which chunk am i in... */
  int next; /* the next element in containedIn to be filled... */
  int theCandidate; /* a point possibly close enough to be tested against
		       the ellipsoidFunction... */
  double theOffset[3]; /* how much should i shift this point to make the 
			  copy i want... */
  double thePoint[3]; /* the point to be tested against ellipsoidFunction.
			 it could be shifted by a boxwidth in any direction. */
  int xChunk,yChunk,zChunk; /* which interval in x, y, and z am i in... */

  theOccupants = NULL;
  nOccupants = 0;
  head = NULL;

  me = theCells[index];
  findThisInterval(count,theBbox,0,xBuffer,
		   me->point[0] - me->e.axisLengths[0],
		   me->point[0] + me->e.axisLengths[0],
		   &nxr,ixr,xOffset);

  findThisInterval(count,theBbox,1,yBuffer,
		   me->point[1] - me->e.axisLengths[0],
		   me->point[1] + me->e.axisLengths[0],
		   &nyr,iyr,yOffset);

  findThisInterval(count,theBbox,2,zBuffer,
		   me->point[2] - me->e.axisLengths[0],
		   me->point[2] + me->e.axisLengths[0],
		   &nzr,izr,zOffset);



  /* scan the smallest interval... */

  if (intervalCount(nxr,ixr) < intervalCount(nyr,iyr)) {
    if (intervalCount(nxr,ixr) < intervalCount(nzr,izr)) {
      smallestInterval = 0;
    }
    else {
      smallestInterval = 2;
    }
  }
  else {
    if (intervalCount(nyr,iyr) < intervalCount(nzr,izr)) {
      smallestInterval = 1;
    }
    else {
      smallestInterval = 2;
    }
  }

  if (smallestInterval == 0) {
    /* traverse the x interval, and build the collection of points that
       fall in all three intervals, keeping track of the appropriate 
       offsets... */

    for(nchunk=0;nchunk<nxr;nchunk++) {
      xChunk = nchunk+1;
    
      for(i=ixr[2*nchunk];i<ixr[2*nchunk+1];i++) {
	theCandidate = xBuffer->indices[i];
	if ((yChunk = inInterval(yBuffer->reverseIndices[theCandidate],
				nyr,iyr))) {
	  if ((zChunk = inInterval(zBuffer->reverseIndices[theCandidate],
				  nzr,izr))) {
	  
	    /* ok, this point *could* be inside the ellipsoid.
	       build the point and calculate its ellipsoid function... */
	    
	    theOffset[0] = xOffset[(xChunk-1)*2];
	    theOffset[1] = yOffset[(yChunk-1)*2];
	    theOffset[2] = zOffset[(zChunk-1)*2];

	    for(j=0;j<3;j++) {
	      thePoint[j] = theCells[theCandidate]->point[j] + theOffset[j];
	    }
	    theEFunc = ellipsoidFunction(thePoint,me);
	    
	    if (theEFunc <= 1.0) {
	      /* whoo-hoo, i'm inside... */
	      
	      theOccupant = (struct occupant *)malloc(sizeof(struct occupant));
	      theOccupant->id = theCandidate;
	      theOccupant->eFunc = theEFunc;
	      theOccupant->offset[0] = theOffset[0];
	      theOccupant->offset[1] = theOffset[1];
	      theOccupant->offset[2] = theOffset[2];
	      
	      head = 
		(struct occupantList *)malloc(sizeof(struct occupantList));
	      head->theOccupant = theOccupant;
	      head->next = theOccupants;
	      theOccupants = head;
	      nOccupants++;
	    }
	  }
	}
      }
    }
  }


  if (smallestInterval == 1) {

    /* traverse the y interval, and build the collection of points that
       fall in all three intervals, keeping track of the appropriate 
       offsets... */

    for(nchunk=0;nchunk<nyr;nchunk++) {
      yChunk = nchunk+1;
    
      for(i=iyr[2*nchunk];i<iyr[2*nchunk+1];i++) {
	theCandidate = yBuffer->indices[i];
	if ((xChunk = inInterval(xBuffer->reverseIndices[theCandidate],
				nxr,ixr))) {
	  if ((zChunk = inInterval(zBuffer->reverseIndices[theCandidate],
				  nzr,izr))) {
	  
	    /* ok, this point *could* be inside the ellipsoid.
	       build the point and calculate its ellipsoid function... */
	    
	    theOffset[0] = xOffset[(xChunk-1)*2];
	    theOffset[1] = yOffset[(yChunk-1)*2];
	    theOffset[2] = zOffset[(zChunk-1)*2];

	    for(j=0;j<3;j++) {
	      thePoint[j] = theCells[theCandidate]->point[j] + theOffset[j];
	    }
	    theEFunc = ellipsoidFunction(thePoint,me);
	    
	    if (theEFunc <= 1.0) {
	      /* whoo-hoo, i'm inside... */
	      
	      theOccupant = (struct occupant *)malloc(sizeof(struct occupant));
	      theOccupant->id = theCandidate;
	      theOccupant->eFunc = theEFunc;
	      theOccupant->offset[0] = theOffset[0];
	      theOccupant->offset[1] = theOffset[1];
	      theOccupant->offset[2] = theOffset[2];
	      
	      head = 
		(struct occupantList *)malloc(sizeof(struct occupantList));
	      head->theOccupant = theOccupant;
	      head->next = theOccupants;
	      theOccupants = head;
	      nOccupants++;
	    }
	  }
	}
      }
    }
  }


  if (smallestInterval == 2) {

    /* traverse the z interval, and build the collection of points that
       fall in all three intervals, keeping track of the appropriate 
       offsets... */

    for(nchunk=0;nchunk<nzr;nchunk++) {
      zChunk = nchunk+1;
    
      for(i=izr[2*nchunk];i<izr[2*nchunk+1];i++) {
	theCandidate = zBuffer->indices[i];
	if ((xChunk = inInterval(xBuffer->reverseIndices[theCandidate],
				nxr,ixr))) {
	  if ((yChunk = inInterval(yBuffer->reverseIndices[theCandidate],
				   nyr,iyr))) {
	    
	    /* ok, this point *could* be inside the ellipsoid.
	       build the point and calculate its ellipsoid function... */
	    
	    
	    theOffset[0] = xOffset[(xChunk-1)*2];
	    theOffset[1] = yOffset[(yChunk-1)*2];
	    theOffset[2] = zOffset[(zChunk-1)*2];

	    for(j=0;j<3;j++) {
	      thePoint[j] = theCells[theCandidate]->point[j] + theOffset[j];
	    }
	    theEFunc = ellipsoidFunction(thePoint,me);
	    
	    if (theEFunc <= 1.0) {
	      /* whoo-hoo, i'm inside... */
	      
	      theOccupant = (struct occupant *)malloc(sizeof(struct occupant));
	      if (theOccupant == NULL) {
		fprintf(stderr,"Malloc failure...\n");
	      }
	      theOccupant->id = theCandidate;
	      theOccupant->eFunc = theEFunc;
	      theOccupant->offset[0] = theOffset[0];
	      theOccupant->offset[1] = theOffset[1];
	      theOccupant->offset[2] = theOffset[2];
	      
	      head = 
		(struct occupantList *)malloc(sizeof(struct occupantList));
	      if (head == NULL) {
		fprintf(stderr,"Malloc failure...\n");
	      }
	      head->theOccupant = theOccupant;
	      head->next = theOccupants;
	      theOccupants = head;
	      nOccupants++;
	    }
	  }
	}
      }
    }
  }

  /* ok, now lets free the linked list and build the array of occupants... */

  me->containedIn = 
    (struct occupant **)malloc(nOccupants*sizeof(struct occupant *));
  if (me->containedIn == NULL) {
    fprintf(stderr,"Malloc failure... \n");
  }

  next = 0;
  while (theOccupants != NULL) {
    me->containedIn[next++] = theOccupants->theOccupant;
    head = theOccupants;
    theOccupants = theOccupants->next;
    free(head);
  }

  me->nContained = nOccupants;


}

int bisectFind(int count, double value, struct xyzBuffer *theBuffer) {

  /* find the abscissa in theBuffer immediately less than value... */

  int low,mid,high; /* the bisection bracket... */

  if (value <= theBuffer->abscissas[0]) {
    return(0);
  }
  else {
    low = 0;
  }
  if (value >= theBuffer->abscissas[count-1]) {
    return(count-1);
  }
  else {
    high = count-1;
  }
  while ((high-low) > 1) {
    mid = (high+low)/2;
    if (theBuffer->abscissas[mid] > value) {
      high = mid;
    }
    else {
      low = mid;
    }
  }
  return(low);
}


int intervalCount(int nsegs,int *intervals) {
  
  int retval;

  retval = intervals[1] - intervals[0];
  if (nsegs > 1) {
    retval += intervals[3] - intervals[2];
  }
  return(retval);
}

int inInterval(int probe,int nsegs,int *intervals) {

  if ((probe >= intervals[0]) && (probe < intervals[1])) {
    return(1);
  }
  if (nsegs > 1) {
    if ((probe >= intervals[2]) && (probe < intervals[3])) {
      return(2);
    }
  }
  return(0);
}


void findThisInterval(int nCells,
		      struct bbox *theBbox,
		      int theCoord,
		      struct xyzBuffer *theBuffer,
		      double vmin,
		      double vmax,
		      int *nsegs,int *breakPoints,
		      double *offsets) {

  /* given a min and max value, this routine finds the interval(s) in
     theBuffer that fall between the min and max value.  this task is 
     complicated by the fact that space is toroidal, so if vmin or vmax 
     hang over either end, they wrap around at the other end.  offsets 
     keeps track of the offset necessary for this wrapping operation... */


  int lowBracket; /* the index of the first element of theBuffer greater
		     than vmin... */
  int highBracket; /* the index of the first element of theBuffer greater 
		      than vmax... */

  *nsegs = 1;
  if (vmin < theBbox->ranges[theCoord][0]) {
    if (vmax > theBbox->ranges[theCoord][1]) {
      *nsegs = 1;
      breakPoints[0] = 0;
      breakPoints[1] = nCells;
      offsets[0] = 0.0;
      offsets[1] = 0.0;
      return;
    }
    *nsegs = 2;
    offsets[0] = theBbox->ranges[theCoord][1] - theBbox->ranges[theCoord][0];
    offsets[1] = offsets[0];
    lowBracket = bisectFind(nCells,vmin+offsets[0],theBuffer);
    breakPoints[0] = lowBracket;
    breakPoints[1] = nCells;
    highBracket = bisectFind(nCells,vmax,theBuffer);
    offsets[2] = 0.0;
    offsets[3] = 0.0;
    breakPoints[2] = 0;
    breakPoints[3] = highBracket+1;
    return;
  }
  else {
    lowBracket = bisectFind(nCells,vmin,theBuffer);
    if (vmax > theBbox->ranges[theCoord][1]) {
      *nsegs = 2;
      offsets[0] = 0.0;
      offsets[1] = 0.0;
      breakPoints[0] = lowBracket;
      breakPoints[1] = nCells;
      offsets[2] = theBbox->ranges[theCoord][0] - theBbox->ranges[theCoord][1];
      offsets[3] = offsets[2];
      highBracket = bisectFind(nCells,vmax + offsets[2],theBuffer);
      breakPoints[2] = 0;
      breakPoints[3] = highBracket+1;
      return;
    }
  }
  *nsegs = 1;
  highBracket = bisectFind(nCells,vmax,theBuffer);
  offsets[0] = 0.0;
  offsets[1] = 0.0;
  breakPoints[0] = lowBracket;
  breakPoints[1] = highBracket+1;
  return;
}
  

void sortBuffer(int count,struct xyzBuffer *buffer) {
  
  /* modified from numeric recipes heapsort routine sort2... */

  int l,j,ir,i;
  double rra;
  int rrb;

  l = (count >> 1)+1;
  ir = count;
  for (;;) {
    if (l > 1) {
      rra = buffer->abscissas[--l-1];
      rrb = buffer->indices[l-1];
    }
    else {
      rra = buffer->abscissas[ir-1];
      rrb = buffer->indices[ir-1];
      buffer->abscissas[ir-1] = buffer->abscissas[0];
      buffer->indices[ir-1] = buffer->indices[0];
      if (--ir == 1) {
	buffer->abscissas[0] = rra;
	buffer->indices[0] = rrb;
	return;
      }
    }
    i = l;
    j = l << 1;
    while (j <= ir) {
      if (j < ir && buffer->abscissas[j-1] < buffer->abscissas[j]) ++j;
      if (rra < buffer->abscissas[j-1]) {
	buffer->abscissas[i-1] = buffer->abscissas[j-1];
	buffer->indices[i-1] = buffer->indices[j-1];
	j += (i=j);
      }
      else j=ir+1;
    }
    buffer->abscissas[i-1] = rra;
    buffer->indices[i-1] = rrb;
  }
}


void reorderAxes(struct foamCell *theCell) {

  /* reorders the axes and axes lengths, so that the longest axis comes
     first... */

  int i,j; /* trusty indices... */
  int order[3]; /* the appropriate order for the axisLengths... */
  struct ellipsoid myE; /* the rearranged ellipsoid... */
  struct ellipsoid *e; /* shorthand for theCell->e */

  e = &(theCell->e);

  if (e->axisLengths[0] >= e->axisLengths[1]) {
    if (e->axisLengths[0] >= e->axisLengths[2]) {
      order[0] = 0;
      if (e->axisLengths[1] >= e->axisLengths[2]) {
	order[1] = 1;
	order[2] = 2;
      }
      else {
	order[1] = 2;
	order[2] = 1;
      }
    }
    else {
      order[0] = 2;
      order[1] = 0;
      order[2] = 1;
    }
  }
  else {
    if (e->axisLengths[1] >= e->axisLengths[2]) {
      order[0] = 1;
      if (e->axisLengths[0] >= e->axisLengths[2]) {
	order[1] = 0;
	order[2] = 2;
      }
      else {
	order[1] = 2;
	order[2] = 0;
      }
    }
    else {
      order[0] = 2;
      order[1] = 1;
      order[2] = 0;
    }
  }

  /* now build the temporary copy in myE... */

  for(i=0;i<3;i++) {
    myE.axisLengths[i] = e->axisLengths[order[i]];
    for(j=0;j<3;j++) {
      myE.axes[i][j] = e->axes[order[i]][j];
    }
  }

  /* and copy over the reordered info... */

  for(i=0;i<3;i++) {
    if (order[i] != i) {
      e->axisLengths[i] = myE.axisLengths[i];
      for(j=0;j<3;j++) {
	e->axes[i][j] = myE.axes[i][j];
      }
    }
  }
}


void twiddle(int nCells,struct foamCell **theCells,int *activeList) {
  
  /* chooses an operator randomly, and performs it.  modifies the
     global variable currentOperator, to allow the logging routine to
     keep track of what operators are being used... */

  static int initialized=0; /* have i initialized the operator probs? */
  static double opMarginal[4]; /* the marginal probility distribution for 
				  selecting the operator... */
  static int abort=0; /* does the file "enough.already" exist? */
  double deviate; /* a uniform random deviate... */
  FILE *testfile; /* used to sense the existence of "enough.already"... */

  if (!initialized) {
    initialized = 1;

    /* opMarginal[0] corresponds to an add operation,
       opMarginal[1] corresponds to a subtract operation,
       opMarginal[2] corresponds to a swap operation,
       opMarginal[3] corresponds to a jog operation... */

    opMarginal[0] = 0.25;
    opMarginal[1] = 0.50;
    opMarginal[2] = 0.75;
    opMarginal[3] = 1.0;
  }

  testfile = fopen("enough.already","r");
  if (testfile != NULL) {
    abort = 1;
    fclose(testfile);
  }

  if (abort) {
    return;
  }
  
  genrandom(&ranseed,1,&deviate);
  
  if (deviate < opMarginal[0]) {
    add(nCells,theCells,activeList);
    currentOperation = 0;
  }
  if ((deviate > opMarginal[0]) && (deviate <= opMarginal[1])) {
    delete(nCells,theCells,activeList);
    currentOperation = 1;
  }
  if ((deviate > opMarginal[1]) && (deviate <= opMarginal[2])) {
    swap(nCells,theCells,activeList);
    currentOperation = 2;
  }
  if (deviate > opMarginal[2]) {
    jog(nCells,theCells,activeList);
    currentOperation = 3;
  }
}
  

void output(int nCells,struct foamCell **theCells,int *activeList,
	    struct bbox *theBbox){

  /* output eliminates any points that aren't in any active cell,
     eliminates any points in multiple active cells and replaces them
     with nearby points in individual cells.  output then writes out 
     the cells (with their data updated to the information from the master
     cell), and the list of active cells, along with which points they contain,
     and then 27 different copies of the points for tesselation (to allow 
     for periodicity)...*/

  struct foamCell **newCells; /* new cells must be created to resolve 
				 points in multiple cells... */

  double suicide,murder; /* for comparing the cost of two alternatives for
			    eliminating active cells whose centers are contained
			    in other active cells... */
  int next; /* the next point to fill... */
  int npts; /* the total number of points... */
  int nNew; /* the number of new cells (for resolving points in multiple cells)... */
  int activeCount; /* the number of active cells... */
  int *replacements; /* the indices of the replacement cells in newCells... */
  int *newNames; /* the renumbered indices... */
  double **thePoints; /* the list of points... */
  int newInd; /* an index into newCells... */
  int candidate; /* an old cell id... */

  int i,j,k,el; /* trusty indices... */
  double dv[3]; /* half the bounding box... */
  int replicantCount; /* the total number of replicated (and original cells)
			 in paddedPts... */

  FILE *allCells; /* the elliptical.cells file, listing all the modified cell 
		     information... */
  FILE *cellist; /* the active.cells file, listing all the active cells and which
		    points they contain... */
  FILE *metaData; /* the number of points and the bounding box... */
  FILE *paddedPoints; /* all the real points and enough padding to eliminate
			 boundary effects... */
  FILE *pointKey; /* a key for determining which points are equivalent... */


  /* repair any situations wherein an active cell contains 
     the center of another active cell... */

  for(i=0;i<nCells;i++) {
    if (theCells[i]->actives > 1) {
      if (activeList[i]) {

	/* ok, the center of an active cell is in some other cell.
	   somebody has to die... */
	
	suicide = 0.0;
	murder = 0.0;
	for(j=0;j<theCells[i]->actives;j++) {
	  if (theCells[i]->activeEngulfers[j] != i) {
	    murder += subtractCost(theCells[i]->activeEngulfers[j],
				   theCells,activeList);
	  }
	  else {
	    suicide += subtractCost(i,theCells,activeList);
	  }
	}
	if (suicide <= murder) {
	  commitSubtract(i,theCells,activeList);
	  systemState += suicide;
	}
	else {
	  for(j=0;j<theCells[i]->actives;j++) {
	    if (theCells[i]->activeEngulfers[j] != i) {
	      commitSubtract(theCells[i]->activeEngulfers[j],
			     theCells,activeList);
	    }
	  }
	  systemState += murder;
	}
      }
    }
  }

  /* count all the active cells, and count any points that appear in more than one 
     cell... */
  

  nNew = 0;
  activeCount = 0;
  for(i=0;i<nCells;i++) {
    if (activeList[i]) {
      activeCount++;
    }
    if (theCells[i]->actives > 1) {
      nNew += theCells[i]->actives;
    }
  }

  newCells = (struct foamCell **)malloc(nNew*sizeof(struct foamCell *));
  replacements = (int *)malloc(nCells*sizeof(int));
  for(i=0;i<nCells;i++) {
    replacements[i] = -1;
  }

  nNew = 0;

  /* now generate multiple cells for all the points contained in more than one
     active cell... */

  for(i=0;i<nCells;i++) {
    if (theCells[i]->actives > 1) {
      genNewCells(newCells,&nNew,i,theCells,replacements);
    }
  }

  /* ok, now lets figure out the new ordering of points, and output the cells... */

  allCells = fopen("elliptical.cells","w");

  newNames = (int *)malloc(nCells*sizeof(int));
  for(i=0;i<nCells;i++) {
    newNames[i] = -1;
  }

  thePoints = (double **)malloc((nCells + nNew)*sizeof(double *));
  next = 0;
  for(i=0;i<nCells;i++) {
    if (theCells[i]->actives == 1) {
      newNames[i] = next;
      thePoints[next] = (double *)malloc(3*sizeof(double));
      thePoints[next][0] = theCells[i]->point[0];
      thePoints[next][1] = theCells[i]->point[1];
      thePoints[next][2] = theCells[i]->point[2];
      outputCell(allCells,0,theCells[i],theCells,theBbox);
      next++;
    }
  }
  for(i=0;i<nNew;i++) {
    newCells[i]->id = next;
    thePoints[next] = (double *)malloc(3*sizeof(double));
    thePoints[next][0] = newCells[i]->point[0];
    thePoints[next][1] = newCells[i]->point[1];
    thePoints[next][2] = newCells[i]->point[2];
    outputCell(allCells,1,newCells[i],theCells,theBbox);
    next++;
  }

  npts = next;

  cellist = fopen("active.list","w");

  /* ok, now lets output the active cells, and which other points they contain... */

  fprintf(stderr,"Active Cell Count: %d\n",activeCount);
  fprintf(cellist,"%d\n",activeCount);
  for(i=0;i<nCells;i++) {
    if (activeList[i]) {
      fprintf(cellist,"%d\n",newNames[i]);
      fprintf(cellist,"%d\n",theCells[i]->nContained);
      for(j=0;j<theCells[i]->nContained;j++) {
	candidate = theCells[i]->containedIn[j]->id;
	if (replacements[candidate] >= 0) {
	  newInd = replacements[candidate];
	  while(newCells[newInd]->activeEngulfers[0] != i) {
	    newInd++;
	  }
	  fprintf(cellist,"%d ",newCells[newInd]->id);
	}
	else {
	  if (newNames[candidate] >= 0) {
	    fprintf(cellist,"%d ",newNames[candidate]);
	  }
	  else {
	    fprintf(stderr,"cell %d has no newName...\n",candidate);
	    fprintf(cellist,"%d ",newNames[candidate]);
	  }
	}
      }
      fprintf(cellist,"\n");
    }
  }

  fclose(cellist);

  /* now lets output a padded set of points to ensure periodicity... */

  for(i=0;i<3;i++) {
    dv[i] = (theBbox->ranges[i][1] - theBbox->ranges[i][0])/2.0;
  }

  /* output all of the real points first... */

  metaData = fopen("metaData","w");
  fprintf(metaData,"%d\n",npts);
  fprintf(metaData,"%f %f %f %f %f %f\n",theBbox->ranges[0][0],
	  theBbox->ranges[0][1],theBbox->ranges[1][0],
	  theBbox->ranges[1][1],theBbox->ranges[2][0],
	  theBbox->ranges[2][1]);
  fclose(metaData);
  paddedPoints = fopen("paddedPoints","w");
  pointKey = fopen("pointKey","w");

  /* leave a blank buffer to write the dimension and the count into... */
  fprintf(paddedPoints,"                                     \n");
  replicantCount = 0;
  for(el=0;el<npts;el++) {
    fprintf(paddedPoints,"%f %f %f\n",thePoints[el][0],
	    thePoints[el][1],thePoints[el][2]);
    fprintf(pointKey,"%d 0 0 0\n",el);
    replicantCount++;
  }

  for(i=-1;i<2;i++) {
    for(j=-1;j<2;j++) {
      for(k=-1;k<2;k++) {
    	if ((fabs(i)+fabs(j)+fabs(k)) > 0) {
	  for(el=0;el<npts;el++) {
	    if ((thePoints[el][0] + 2.0*i*dv[0] > 
		 (theBbox->ranges[0][0] - paddingFraction*2.0*dv[0])) &&
		(thePoints[el][0] + 2.0*i*dv[0] <
		 (theBbox->ranges[0][1] + paddingFraction*2.0*dv[0])) &&
		(thePoints[el][1] + 2.0*j*dv[1] > 
		 (theBbox->ranges[1][0] - paddingFraction*2.0*dv[1])) &&
		(thePoints[el][1] + 2.0*j*dv[1] <
		 (theBbox->ranges[1][1] + paddingFraction*2.0*dv[1])) &&
		(thePoints[el][2] + 2.0*k*dv[2] > 
		 (theBbox->ranges[2][0] - paddingFraction*2.0*dv[2])) &&
		(thePoints[el][2] + 2.0*k*dv[2] <
		 (theBbox->ranges[2][1] + paddingFraction*2.0*dv[2]))) {  
	      fprintf(paddedPoints,"%f ",thePoints[el][0] + 2.0*i*dv[0]);      
	      fprintf(paddedPoints,"%f ",thePoints[el][1] + 2.0*j*dv[1]);
	      fprintf(paddedPoints,"%f ",thePoints[el][2] + 2.0*k*dv[2]);
	      fprintf(paddedPoints,"\n");
	      fprintf(pointKey,"%d %d %d %d\n",el,i,j,k);
	      replicantCount++;
	    }
	  }
	}
      }
    }
  }	
  fseek(paddedPoints,0,SEEK_SET);
  fprintf(paddedPoints,"%d %d",3,replicantCount);
  fclose(paddedPoints); 
  fclose(pointKey);

}


void genNewCells(struct foamCell **newCells,int *next,
		 int index,struct foamCell **oldCells,int *replacements) {

  /* given oldCells[index], which is contained in more than one active cell,
     generate a collection of points near oldCells[index]->point that are each in
     just one cell, and create cells around them... */

  int nActive; /* the number of active cells oldCells[index] was in... */
  int insider; /* the cell we're generating a pt inside... */
  int i,j; /* trusty indices... */
  double high,low,mid; /* the bounds of a bisection interval... */
  double pt[3]; /* the current bisection point... */
  int insiderIndex; /* used for looking up the correct offset... */
  int stillInside; /* used for determining which half of the bisection i'm in... */
  double theFunc; /* the ellipsoid function in one of the enclosing ellipsoids... */
  struct foamCell *me; /* shorthand for newCells[*next]... */
  static double abortTol=0.001; /* stop bisecting here... */

  replacements[index] = *next;
  nActive = oldCells[index]->actives;
  for(i=0;i<nActive;i++) {
    insider = oldCells[index]->activeEngulfers[i];
    findEFunc(oldCells[insider],index,&insiderIndex);
    if (insiderIndex) {
      insiderIndex--;
    }
    else {
      fprintf(stderr,"Inconsistent state in genNewCells.  cell %d thinks it's inside\n",
	      index);
      fprintf(stderr,"cell %d, but the feeling is not reciprocated...\n",insider);
    }
    high = 1.0;
    low = 0.0;

    while ((high - low) > abortTol) { 
      mid = (high + low)/2.0;

      for(j=0;j<3;j++) {
	pt[j] = oldCells[insider]->point[j] + 
	  mid*(oldCells[index]->point[j] + 
	       oldCells[insider]->containedIn[insiderIndex]->offset[j] - 
	       oldCells[insider]->point[j]);
      }

      stillInside = 0;
      for(j=0;j<nActive;j++) {
	if (j != i) {
	  theFunc = ellipsoidFunction(pt,oldCells[oldCells[index]->activeEngulfers[j]]);
	  if ((theFunc >= 0.0) && (theFunc < 1.0)) {
	    stillInside = 1;
	  }
	}
      }
      if (stillInside) {
	high = mid;
      }
      else {
	low = mid;
      }
    }
    for(j=0;j<3;j++) {
      pt[j] = oldCells[insider]->point[j] + 
	low*(oldCells[index]->point[j] +
	     oldCells[insider]->containedIn[insiderIndex]->offset[j] - 
	     oldCells[insider]->point[j]);
    }
    newCells[*next] = (struct foamCell *)malloc(sizeof(struct foamCell));
    me = newCells[*next];
    for(j=0;j<3;j++) {
      me->point[j] = pt[j];
      me->orientation[j] = oldCells[insider]->orientation[j];
      me->e.axisLengths[j] = oldCells[insider]->e.axisLengths[j];
      me->e.axes[0][j] = oldCells[insider]->e.axes[0][j];
      me->e.axes[1][j] = oldCells[insider]->e.axes[1][j];
      me->e.axes[2][j] = oldCells[insider]->e.axes[2][j];
    }
    me->actives = 1;
    me->activeEngulfers[0] = insider;
    me->nContained = 0;
    (*next)++;
  }
}
    

void outputCell(FILE *outfile,int noOffset,struct foamCell *theCell,
		struct foamCell **theCells,struct bbox *theBbox) {

  /* modifies all cell i's information to match the cell that contains it,
     translates cell i, if necessary, and writes the data to outfile... */

  
  int index; /* the index of cell i in its engulfer's containedIn array... */
  int engulfer; /* the index of my engulfer... */
  int j,k; /* trusty indices... */

  if (theCell->actives == 1) {
    engulfer = theCell->activeEngulfers[0];
    if (!noOffset) {
      findEFunc(theCells[engulfer],theCell->id,&index);
      if (index) {
	
	/* findEFunc actually passes index+1 back... */
	index--;
      }
      else {
	fprintf(stderr,"Confusion in outputCell.  theCells[%d] says %d contains it\n",
		theCell->id,engulfer);
	fprintf(stderr,"but cell %d doesn't know a thing about it...\n",engulfer);
	noOffset = 1;
      }
    }
    if (!noOffset) {
      fprintf(outfile,"%f %f %f\n",
	      theCell->point[0] + theCells[engulfer]->containedIn[index]->offset[0],
	      theCell->point[1] + theCells[engulfer]->containedIn[index]->offset[1],
	      theCell->point[2] + theCells[engulfer]->containedIn[index]->offset[2]);
    }
    else {
      fprintf(outfile,"%f %f %f\n",theCell->point[0],theCell->point[1],
	      theCell->point[2]);
    }
    fprintf(outfile,"%f %f %f\n",theCells[engulfer]->e.axisLengths[0],
	    theCells[engulfer]->e.axisLengths[1],theCells[engulfer]->e.axisLengths[2]);
    for(j=0;j<3;j++) {
      for(k=0;k<3;k++) {
	fprintf(outfile,"%f ",theCells[engulfer]->e.axes[j][k]);
      }
      fprintf(outfile,"\n");
    }
    fprintf(outfile,"%f %f %f\n",theCells[engulfer]->orientation[0],
	    theCells[engulfer]->orientation[1],theCells[engulfer]->orientation[2]);
  }
  else {
    fprintf(stderr,"Confusion in outputCell.  Cell %d thinks it has %d engulfers...\n",
	    theCell->id,theCell->actives);
  }
}




void writeOutput(int nCells,struct foamCell **theCells,int *activeList) {

  FILE *outfile;
  int nActives;
  int i,j; 

  /* count the number of active cells... */
  
  nActives = 0;
  for(i=0;i<nCells;i++) {
    if (activeList[i]) {
      nActives++;
    }
  }

  outfile = fopen("activeFile","w");
  fprintf(outfile,"%d\n",nActives);
  for(i=0;i<nCells;i++) {
    if (activeList[i]) {
      fprintf(outfile,"%d\n",i);
      fprintf(outfile,"%d\n",theCells[i]->nContained);
      for(j=0;j<theCells[i]->nContained;j++) {
	fprintf(outfile,"%d ",theCells[i]->containedIn[j]->id);
      }
      fprintf(outfile,"\n");
    }
  }
}

int main(int argc, char *argv[]) {

  /* read a set of points, with associated cell sizes and orientations.
     then generate a consistent set of foamCells (which know which points are
     inside which cells).  Then calculate the average number of points in each
     cell.  this defines an approximation for the total number of active 
     cells (npts/avgPerCell).  call add that number of times, and then start
     a random rotation of add, subtract, swap and jog.  stop when you've 
     finished the annealing cycle, or when a file called "enough.already" 
     appears in the current directory.  upon completion, eliminate any points
     in multiple cells, replacing them with collections of nearby points in 
     enclosing cells.  output a new set of points with associated cell sizes
     and orientations, where each point takes the desired orientation of its
     cell.  also output a list of active cells, with the points they contain...
  */

  FILE *infile; /* the input file describing morphology and texture... */
  FILE *ctrlFile; /* a file containing descriptive and control information...*/
  int abort; /* no valid input file... */
  int fullhelp; /* give full help message? */
  struct foamCell **theCells; /* the foam cells... */
  int nCells; /* the total number of cells... */
  int i; /* trusty index... */
  double sumContained; /* the sum of all the nContained fields... */
  double avgContained; /* the avg of all the nContained fields... */
  double sumEngulfers; /* the sum of all the nEngulfers fields... */
  double avgEngulfers; /* the avg of all the nEngulfers fields... */
  int *activeList; /* the list of who's active and who's not... */
  int startingSet; /* an approximation of the number of desired 
		      active cells... */
  double xmin,xmax; /* if the "cell.ctrl" exists, these variables are used
		       to read the bounding box... */
  double ymin,ymax; /* if the "cell.ctrl" exists, these variables are used
		       to read the bounding box... */
  double zmin,zmax; /* if the "cell.ctrl" exists, these variables are used
		       to read the bounding box... */
  int nImmortals; /* if "cell.ctrl exists, this variable holds the number of
		     mandatory cells... */
  int nextImmortal; /* the next mandatory cell... */
  double nextParam; /* the next control parameter... */
  double cost; /* used to insert mandatory cells... */

  if (argc > 1) {
    infile = fopen(argv[argc-1],"r");
    if (infile != NULL) {
      abort = 0;
    }
    else {
      abort = 1;
    }
  }
  else {
    abort = 1;
  }
  if ((argc>=2) && (strcmp(argv[1],"-h") == 0)) {
    fullhelp = 1;
  }
  else {
    fullhelp = 0;
  }

  if (abort && fullhelp) {

    fprintf(stderr,"Usage::ellipticalFoam inFile\n");
    fprintf(stderr,"       where\n");
    fprintf(stderr,
	    "            inFile contains a full description of the desired\n");
    fprintf(stderr,
	    "            cell geometry, as a function of randomly sampled\n");
    fprintf(stderr,
	    "            spatial coordinates.\n");
    fprintf(stderr,
	    "            Foreach of npts sample points:\n");
    fprintf(stderr,
	    "               point (3 doubles)\n");
    fprintf(stderr,
	    "               axisLengths (3 doubles)\n");
    fprintf(stderr,
	    "               axes (three unit vectors, 9 doubles in all)\n");
    fprintf(stderr,
	    "               crystallographic orientation (3 doubles)\n\n");
    
    fprintf(stderr,
	    "          inFile defines morphology and texture as a \n");
    fprintf(stderr,
	    "          (sampled) function of space.  The ellipsoidal cells\n");
    fprintf(stderr,
	    "          defined around the points in inFile overlap (and \n");
    fprintf(stderr,
	    "          could conceivably have gaps too).  ellipticalFoam\n");
    fprintf(stderr,
	    "          chooses a subset of the cells in inFile, and labels\n");
    fprintf(stderr,
	    "          them active.  Each active cell will be represented \n");
    fprintf(stderr,
	    "          as the union of the voronoi regions of all the \n");
    fprintf(stderr,
	    "          points inside the cell.  The trick is to choose the\n");
    fprintf(stderr,
	    "          active cells such that space is filled (as much as \n");
    fprintf(stderr,
	    "          possible), and such that the cells overlap as \n");
    fprintf(stderr,   
	    "          little as possible.  This is accomplished by a \n");
    fprintf(stderr,
	    "          process of simulated annealing.  After the \n");
    fprintf(stderr,
	    "          annealing schedule is completed (or the user \n");
    fprintf(stderr,
	    "          decides to terminate the run), ellipticalFoam will\n");
    fprintf(stderr,
	    "          produce two output files. \n\n");
    fprintf(stderr,
	    "          The first (\"elliptical.cells\") is a file of the \n");
    fprintf(stderr,
	    "          same format as infile, but the data differs, in \n");
    fprintf(stderr,
	    "          that all points within an active cell have the same\n");
    fprintf(stderr,
	    "          crystallographic orientation and cell morphology \n");
    fprintf(stderr,
	    "          as the master point whose cell they're in.\n");
    fprintf(stderr,
	    "          Points contained in more than one cell are \n");
    fprintf(stderr,
	    "          eliminated and replaced with collections of nearby \n");
    fprintf(stderr,
	    "          points in individual cells.  Points not contained \n");
    fprintf(stderr,
	    "          in any active cell are also eliminated \n");
    fprintf(stderr,
	    "          (and not replaced with anything).\n\n");
    fprintf(stderr,
	    "          The second output file (\"active.list\") is a list \n");
    fprintf(stderr,
	    "          of active cells, followed by a list of all the \n");
    fprintf(stderr,
	    "          points each cell contains.  The indices used to \n");
    fprintf(stderr,
	    "          identify points in \"active.list\" are indices into\n");
    fprintf(stderr,
	    "          \"elliptical.cells\".  More explicitly, the file \n");
    fprintf(stderr,
	    "          format of \"active.list\" is:\n");
    fprintf(stderr,
	    "               nCells (1 int - number of active cells)\n");
    fprintf(stderr,
	    "               Foreach of nCells cells:\n");
    fprintf(stderr,
	    "                   cellID (1 int - index of active cell in\n");
    fprintf(stderr,
	    "                           \"elliptical.cells\")\n");
    fprintf(stderr,
	    "                   nPts (1 int - number of points contained \n");
    fprintf(stderr,
	    "                         cell cellID.  cell cellID will be \n");
    fprintf(stderr,  
	    "                         represented as the union of the \n");
    fprintf(stderr,
	    "                         voronoi regions of these points)\n");
    fprintf(stderr,
	    "                   Foreach of nPts points:\n");
    fprintf(stderr,
	    "                       pointID (1 int - index of this point\n");
    fprintf(stderr,
	    "                                in \"elliptical.cells\".  \n");
    fprintf(stderr,
	    "                                N.B. pointID's and cellID's \n");
    fprintf(stderr,
	    "                                are the same thing - indices\n");
    fprintf(stderr,
	    "                                into \"elliptical.cells\")\n\n");
    fprintf(stderr,
	    "          ellipticalFoam also uses four incidental files \n");
    fprintf(stderr,
	    "          to control and monitor the simulated annealing \n");
    fprintf(stderr,
	    "          process.  The first, \"annealing.schedule\" defines\n");
    fprintf(stderr,
	    "          how many transitions to run, and how fast to cool \n");
    fprintf(stderr,
	    "          the structure.  The annealing schedule is \n");
    fprintf(stderr, 
	    "          a piecewise linear curve, whose y axis defines\n");
    fprintf(stderr,
	    "          that positive transition, (expressed as a fraction\n");
    fprintf(stderr,
	    "          of the average positive transition), which has a \n");
    fprintf(stderr, 
	    "          one in two chance of being accepted, and whose \n");
    fprintf(stderr,
	    "          x axis is the transitionCount.\n");
    fprintf(stderr,
	    "          Explicitly the file format for \n");
    fprintf(stderr,
	    "          \"annealing.schedule\" follows: \n");
    fprintf(stderr,
	    "              nBreaks (1 int - the number of knot points\n");
    fprintf(stderr,
	    "                       in the annealing schedule)\n");
    fprintf(stderr,
	    "              Foreach of nBreaks knot points:\n");
    fprintf(stderr,
	    "                  transitionCount (1 int)\n");
    fprintf(stderr,
	    "                  halfLifeFraction (1 double)\n\n");
    fprintf(stderr,
	    "          For instance, the default annealing schedule,\n");
    fprintf(stderr,
	    "          used if \"annealing.schedule\" does not exist \n");
    fprintf(stderr,
	    "          would be defined by the following file\n\n");
    fprintf(stderr,"          2\n");
    fprintf(stderr,"          0 0.4\n");
    fprintf(stderr,"          100000 0.0001\n\n");
    fprintf(stderr,
	    "          which specifies that at the beginning of the run\n");
    fprintf(stderr,
	    "          positive transitions which were 0.4 times the \n");
    fprintf(stderr,
	    "          average positive transition would be accepted half \n");
    fprintf(stderr,
	    "          the time, but by the end of the run (at 100000 \n");
    fprintf(stderr,
	    "          transitions), positive transitions that were \n");
    fprintf(stderr,
	    "          0.0001 times the average positive transition would \n");
    fprintf(stderr,
	    "          be accepted half the time.\n\n");
    fprintf(stderr,
	    "          The second incidental file is the file \n");
    fprintf(stderr,
	    "          \"enough.already\", which if present terminates the\n");
    fprintf(stderr,
	    "          run and generates the output files.\n\n");
    fprintf(stderr,
	    "          The third incidental file is the file \n");
    fprintf(stderr,
	    "          \"annealing.log\", which records the \n");
    fprintf(stderr,
	    "          transitionCount, cost, and total system energy\n");
    fprintf(stderr,"          for each successful transition.\n\n");
    fprintf(stderr,
	    "          The fourth incidental file is \"cell.ctrl\".\n");
    fprintf(stderr,
	    "          If present, \"cell.ctrl\" first specifies \n");
    fprintf(stderr,
	    "                 the bounding box to replicate about.\n");
    fprintf(stderr,
	    "                 (if unspecified, ellipticalFoam uses the \n");
    fprintf(stderr,
	    "                  bounding box of the data).\n");
    fprintf(stderr,
	    "                 The format of the bounding box specification\n");
    fprintf(stderr,
	    "                 follows:\n");
    fprintf(stderr,
	    "                      xmin,xmax (2 doubles)\n");
    fprintf(stderr,
	    "                      ymin,ymax (2 doubles)\n");
    fprintf(stderr,
	    "                      zmin,zmax (2 doubles)\n");
    fprintf(stderr,
	    "          After the bounding box, \"cell.ctrl\" allows the \n");
    fprintf(stderr,
	    "                 user to specify a list of input ellipsoids  \n");
    fprintf(stderr,
	    "                 must be active.\n");
    fprintf(stderr,
	    "                      nMandatory (1 int)\n");
    fprintf(stderr,
	    "                      foreach of nMandatory mandatory cells:\n");
    fprintf(stderr,
	    "                           index (1 int)\n");
    fprintf(stderr,
	    "          After the list of mandatory cells, \"cell.ctrl\" \n");
    fprintf(stderr,
	    "                 allows the user to change the penalties that\n");
    fprintf(stderr,
	    "                 control the simulated annealing process.  \n");
    fprintf(stderr,
	    "                 These are, in order:\n");
    fprintf(stderr,
	    "                      consumeAward (1 double)\n");
    fprintf(stderr,
	    "                      overlapEncouragement (1 double)\n");
    fprintf(stderr,
	    "                      zeroPenalty (1 double)\n");

    exit(1);

  }
  if (abort && !fullhelp) {
    fprintf(stderr,"Usage:: ellipticalFoam infile\n");
    fprintf(stderr,"        reads infile, produces \"elliptical.cells\"\n");
    fprintf(stderr,"        and \"active.list\".\n\n");
    fprintf(stderr,"        For complete help, type:\n\n");
    fprintf(stderr,"ellipticalFoam -h\n");
    exit(1);
  }

  ctrlFile = fopen("cell.ctrl","r");
  if (ctrlFile == NULL) {
    metaInformation = 0;
  }
  else {
    metaInformation = 0;
    if (fscanf(ctrlFile,"%lf %lf",&xmin,&xmax) == 2) {
      if (fscanf(ctrlFile,"%lf %lf",&ymin,&ymax) == 2) {
	if (fscanf(ctrlFile, "%lf %lf",&zmin,&zmax) == 2) {
	  metaInformation = 1;
	  theBbox = (struct bbox *)malloc(sizeof(struct bbox));
    
	  theBbox->ranges[0][0] = xmin;
	  theBbox->ranges[0][1] = xmax;
	  theBbox->ranges[1][0] = ymin;
	  theBbox->ranges[1][1] = ymax;
	  theBbox->ranges[2][0] = zmin;
	  theBbox->ranges[2][1] = zmax;

	}
      }
    }
  }
  
  readPoints(argv[argc-1],&nCells,&theCells);

  immortals = (int *)malloc(nCells*sizeof(int));
  for(i=0;i<nCells;i++) {
    immortals[i] = 0;
  }

  if (metaInformation) {
    if (fscanf(ctrlFile,"%d",&nImmortals) == 1) {
      for(i=0;i<nImmortals;i++) {
	fscanf(ctrlFile,"%d",&nextImmortal);
	if ((nextImmortal >= 0) && (nextImmortal < nCells)) {
	  immortals[nextImmortal] = 1;
	}
      }
    }
    if (fscanf(ctrlFile,"%lf",&nextParam) == 1) {
      consumeAward = nextParam;
    }
    if (fscanf(ctrlFile,"%lf",&nextParam) == 1) {
      overlapEncouragement = nextParam;
    }
    if (fscanf(ctrlFile,"%lf",&nextParam) == 1) {
      zeroPenalty = nextParam;
    }
    fclose(ctrlFile);
  }

  /* allocate the active list, and initialize everyone to off... */

  activeList = (int *)malloc(nCells*sizeof(int));
  for(i=0;i<nCells;i++) {
    activeList[i] = 0;
  }

  /* calculate the average number of points in a cell... */

  sumContained = 0.0;
  sumEngulfers = 0.0;
  for(i=0;i<nCells;i++) {
    sumContained += theCells[i]->nContained;
    sumEngulfers += theCells[i]->nEngulfers;
  }
  avgContained = sumContained/nCells;
  avgEngulfers = sumEngulfers/nCells;

  startingSet = nCells/avgContained;

  for(i=0;i<nCells;i++) {
    if (immortals[i]) {
      cost = addCost(i,theCells,activeList);
      commitAdd(i,theCells,activeList);
      systemState += cost;
      logState(cost);
    }
  }

  for(i=0;i<startingSet;i++) {
    currentOperation = 0;
    add(nCells,theCells,activeList);
  }

  for(i=0;i<maxTransitions;i++) {
    twiddle(nCells,theCells,activeList);
  }

  output(nCells,theCells,activeList,theBbox);
  exit(0);
}
