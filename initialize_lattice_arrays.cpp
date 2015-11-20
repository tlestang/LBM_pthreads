#ifndef __global__
#define __global__
#include "global.h"
#endif

void initializePopulations(double *fin, int Dx, int Dy)
{

  for(int x=0;x<Dx;x++)
    {
      for(int y=0;y<Dy;y++)
	{
	  for(int k=0;k<9;k++)
	    {
	      fin[IDX(x,y,k)] = w[k];
	    }
	}
    }
}

void initializeFields(double *rho, double *ux, double *uy, int Dx, int Dy)
{
  double rho0 = 1.0; double zeroVelocity[2] = {0.0, 0.0};
  for(int x=0;x<Dx;x++)
    {
      for(int y=0;y<Dy;y++)
	{
	  rho[idx(x,y)] = rho0;
	  ux[idx(x,y)] = zeroVelocity[0];
	  uy[idx(x,y)] = zeroVelocity[1];
      	}
    }
}
