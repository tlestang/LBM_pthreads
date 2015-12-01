#include <pthread.h>
#include <iostream>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <ctype.h>
#include <stdio.h>

#ifndef __global__
#define __global__
#include "global.h"
#endif

using namespace std;

#define NUM_THREADS 2  // according to the number of processor cores

void *do_it(void *);

typedef struct shared{
  int xx, yy, kk, xblcksize, yblcksize, tid;
  double *in, *out, *rho, *ux, *uy, *uz;
  double omega, beta;
} shared;

void streamingAndCollision_POSIX(double *fin, double *fout, double *rho, double *ux, double *uy, double beta, double tau)//, int Dx, int Dy)
{
  int xblcksize=Dx, yblcksize=Dy/NUM_THREADS;
  int threadIdx = 0;
  int rc;
  void *status;
  pthread_t thread[NUM_THREADS];
  shared data[NUM_THREADS];
  /*Set thread joinable attribute.*/
  /*Not sure it is required*/
  pthread_attr_t attr;
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  
  for(int xx=0;xx<Dx;xx+=xblcksize+1)
    {

      threadIdx = 0;
      for(int yy=0;yy<Dy;yy+=yblcksize+1)
	{
	  data[threadIdx].xx = xx;
	  data[threadIdx].yy = yy;
	  data[threadIdx].xblcksize = xblcksize;
	  data[threadIdx].yblcksize = yblcksize;
	  data[threadIdx].in = fin;
	  data[threadIdx].out = fout;
	  data[threadIdx].rho = rho;
	  data[threadIdx].ux = ux;
	  data[threadIdx].uy = uy;
	  data[threadIdx].omega = 1.0/tau;
	  data[threadIdx].beta = beta;
	  data[threadIdx].tid = threadIdx;

	  /*Launches threads*/
	  rc = pthread_create(&thread[threadIdx], &attr, &do_it, (void *) &data[threadIdx]);
	  if(rc)
	    {
	      perror("pthread_create"), exit(-1);
	    }
	  threadIdx++;
	}
      /*Now threadIdx equals the total number of threads*/
      /*Wait for all threads to complete*/
      for(int t=0; t<threadIdx;t++)
	{
	  rc = pthread_join(thread[t], &status);
	  if (rc){
	    std::cout << "ERROR: return code from pthread_join() is" << rc << std::endl;
	    exit(-1);
	  }
	}
    }
}

void *do_it(void *arg0)
{
  struct shared *arg = (struct shared *)arg0;
  int xStart, yStart, nx, ny, xblcksize, yblcksize;
  double *fin, *fout, *rho, *ux, *uy;
  double beta, omega;
  double ftemp, feq, rhoDum, uxDum, uyDum, eu, eueu, u2;
  double omega1, coeff_forcing, force_driving;
  
  xStart = arg->xx;
  yStart = arg->yy;

  fin = arg->in;
  fout = arg->out;
  rho = arg->rho;
  ux = arg->ux;
  uy = arg->uy;

  omega = arg->omega;
  omega1 = 1.0-omega; coeff_forcing = 1.0- 0.5*omega;
  beta = arg->beta;
  xblcksize = arg->xblcksize;
  yblcksize = arg->yblcksize;

  for(int x=xStart;x<min(Dx, xStart+xblcksize+1);x++)
    {
      for(int y=yStart;y<min(Dy,yStart+yblcksize+1);y++)
	{
	  rhoDum = 0.0;
	  uxDum = 0.0; uyDum = 0.0;
	  for(int k=0;k<9;k++)
	    {
	      ftemp = fin[IDX(x,y,k)];
	      rhoDum += ftemp;
	      uxDum += ftemp*c[k][0];
	      uyDum += ftemp*c[k][1];
	    }
	  uxDum = uxDum/rhoDum + /*Influence of the force*/0.5*beta/rhoDum; uyDum /= rhoDum;
	  u2 = -1.5*(uxDum*uxDum + uyDum*uyDum);
	  for(int k=0;k<9;k++)
	    {
	      eu = c[k][0]*uxDum + c[k][1]*uyDum;
	      eueu = 4.5*eu*eu;
	      feq = w[k]*rhoDum*(1.0+3.0*eu+eueu+u2);
	      force_driving = w[k]*coeff_forcing*beta*(3.*(c[k][0]-uxDum)+9.*c[k][0]*eu);
	      fin[IDX(x,y,k)] = fin[IDX(x,y,k)]*omega1+feq*omega+force_driving;
	      /*Streaming*/
	      nx = (x + c[k][0]+Dx)%Dx; ny = (y + c[k][1]+Dy)%Dy;
	      fout[IDX(nx,ny,k)] = fin[IDX(x,y,k)];
	    }
	  /*Assign fields values to lattice arrays*/
	  ux[idx(x,y)] = uxDum; uy[idx(x,y)] = uyDum;
	  rho[idx(x,y)] = rhoDum;
	}
    }
  /*Kill thread*/
  pthread_exit(NULL);
}
