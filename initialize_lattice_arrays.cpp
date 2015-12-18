#ifndef __global__
#define __global__
#include "global.h"
#endif
#include <stdlib.h>

void initializePopulationsRandom(double *fin, int Dx, int Dy)
{
  double rho = 1.0;
  double ux, uy, u2, eu, eueu;
  for(int x=0;x<Dx;x++)
    {
      for(int y=0;y<Dy;y++)
	{
	  ux = 1.0e-2 + (5.0e-3)*(2*drand48()-1);
	  uy = 1.0e-2 + (5.0e-3)*(2*drand48()-1);
	  u2 = -1.5*(ux*ux + uy*uy);
	  for(int k=0;k<9;k++)
	    {
	      eu = c[k][0]*ux + c[k][1]*uy;
	      eueu = 4.5*eu*eu;
	      fin[IDX(x,y,k)] = w[k]*rho*(1.0+3.0*eu+eueu+u2);
	    }
	}
    }
}
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
