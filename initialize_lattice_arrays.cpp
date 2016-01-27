#ifndef __global__
#define __global__
#include "global.h"
#endif
#include <stdlib.h>
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

void initializeFields(double *fin, double *rho, double *ux, double *uy, int Dx, int Dy)
{
  double rho_, u, v, f;
  for(int x=0;x<Dx;x++)
    {
      for(int y=0;y<Dy;y++)
	{
	  u = v = rho_ = 0.0;
	  for(int k=0;k<9;k++)
	    {
	      f = fin[IDX(x,y,k)];
	      rho_ += f;
	      u += f*c[k][0];
	      v += f*c[k][1];
	    }
	  rho[idx(x,y)] = rho_;
	  ux[idx(x,y)] = u/rho_;
	  uy[idx(x,y)] = v/rho_;
      	}
    }
}
