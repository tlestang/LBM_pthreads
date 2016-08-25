#include <iostream>

#ifndef __global__
#define __global__
#include "global.h"
#endif

double computeForceOnWalls(double *f, double omega)
{
  double ot = 1./3.;
  double rho0 = 1.0;
  int x0, y0;
  double eu, eueu, u2;
  double fEast, fWest, fNorth, fSouth;
  double rho_, ux, uy, Pi_xx, Pi_xy;
  double fneq, ftemp, feq;
  double totalForce;
  double coeff_force = 1.0-.5*omega;

  /*North side*/
  y0 = Dy-1;
  fNorth = 0.0;
  for(int x=0;x<Dx;x++)
    {
      /*Compute local macro. fields*/
      rho_ = 0.0; ux = 0.0; uy = 0.0;
      for(int k=0;k<9;k++)
	{
	  ftemp = f[IDX(x,y0,k)];
	  rho_ += ftemp;
	  ux += ftemp*c[k][0];
	  uy += ftemp*c[k][1];
	}
#ifdef _INCOMP
      uy /= rho0;
      ux /= rho0;
#else
      uy /= rho_;
      ux /= rho_;
#endif
      /*Compute tensor Pi1*/
      Pi_xy = 0.0;
      u2 = -1.5*(ux*ux + uy*uy);
      for(int k=0;k<9;k++)
	{
	  eu = c[k][0]*ux + c[k][1]*uy;
	  eueu = 4.5*eu*eu;
#ifdef _INCOMP
	  feq = w[k]*(rho_ +/*rho0**/(3.0*eu+eueu+u2));
#else
	  feq = w[k]*rho_*(1.0+3.0*eu+eueu+u2);
#endif
	  fneq = f[IDX(x,y0,k)] - feq;
	  Pi_xy += fneq*c[k][0]*c[k][1];
	}
      fNorth += + coeff_force*Pi_xy;
    }

  /*South side*/
  y0 = 0;
  fSouth = 0.0;
  for(int x=0;x<Dx;x++)
    {
      rho_ = 0.0; ux = 0.0; uy = 0.0;
      /*Compute local macro. fields*/
      for(int k=0;k<9;k++)
	{
	  ftemp = f[IDX(x,y0,k)];
	  rho_ += ftemp;
	  ux += ftemp*c[k][0];
	  uy += ftemp*c[k][1];
	}
#ifdef _INCOMP
      uy /= rho0;
      ux /= rho0;
#else
      uy /= rho_;
      ux /= rho_;
#endif
      /*Compute tensor Pi1*/
      Pi_xy = 0.0;
      u2 = -1.5*(ux*ux + uy*uy);
      for(int k=0;k<9;k++)
	{
	  eu = c[k][0]*ux + c[k][1]*uy;
	  eueu = 4.5*eu*eu;
#ifdef _INCOMP
	  feq = w[k]*(rho_ +/*rho0**/(3.0*eu+eueu+u2));
#else
	  feq = w[k]*rho_*(1.0+3.0*eu+eueu+u2);
#endif
	  fneq = f[IDX(x,y0,k)] - feq;
	  Pi_xy += fneq*c[k][0]*c[k][1];
	}
      fSouth += - coeff_force*Pi_xy;
    }

  totalForce = fNorth + fSouth;
  return totalForce;
}
