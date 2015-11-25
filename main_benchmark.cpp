#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <malloc.h>
#include <unistd.h>
#include <sys/time.h>

#ifndef __global__
#define __global__
#include "global.h"
#endif


#include "initialize_lattice_arrays.h"
#include "streamCollCompute.h"
#include "boundaryConditions.h"
#include "write_vtk.h"

using namespace std;

int main()
{
  /*Parameters for LB simulation*/
  int nbOfTimeSteps, Lx, Ly;
  double tau, beta;
  double *fin, *fout, *temp, *rho, *ux, *uy;
  double Ma;   //Mach number
  const double w[9]={4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};
  
  /*Reads input file*/
      ifstream input_file("input.datin");
      input_file >> nbOfTimeSteps;
      input_file >> Lx; Ly = Lx;
      input_file >> tau;
      input_file >> Ma;
      input_file.close();
      
  /* --- Compute or define other parameters --- */
      int Dy = 4*Ly + 1, Dx = 2*(Dy-1) + 1;
      int xmin = (Dx-1)/2; int xmax = xmin + Lx;
      int ymin = (Dy-1)/2 - Ly/2; int ymax = ymin + Ly;
      double cs = 1./sqrt(3); double rho0 = 1.0;
      double u0 = cs*cs*Ma; double uxSum, uxMean;
      double nu = 1./3.*(tau-0.5);
      double omega = 1.0/tau;
      beta = 8*nu*u0/((Dy-1)/2)/((Dy-1)/2);
  
  /* ---- | Allocate populations and fields | --- */

      fin = (double *) memalign(getpagesize(), Dx*Dy*9*sizeof(double));
      fout = (double *) memalign(getpagesize(), Dx*Dy*9*sizeof(double));
      rho = (double *) memalign(getpagesize(), Dx*Dy*sizeof(double));
      ux = (double *) memalign(getpagesize(), Dx*Dy*sizeof(double));
      uy = (double *) memalign(getpagesize(), Dx*Dy*sizeof(double));

   /* --- Initialize pops to equilibrium value --- */
      initializePopulations(fin, Dx, Dy);
      initializeFields(rho, ux, uy, Dx, Dy);
      int dummy2 = 0;


      /*Define timer variables*/
      struct timeval start, end;
      gettimeofday(&start,NULL);
  /* --- START LBM ---*/
      int tt=0;
	  for (int lbTimeStepCount=0; lbTimeStepCount<nbOfTimeSteps;lbTimeStepCount++)
	    {
	      streamingAndCollision_POSIX(fin, fout, rho, ux, uy, beta, tau, Dx, Dy);
	      computeDomainNoSlipWalls_BB(fout, fin, Dx, Dy);
	      computeSquareBounceBack_TEST(fout, fin, xmin, xmax, ymin, ymax, Dx, Dy);
	      /*Reset square nodes to equilibrium*/
	      for(int x=xmin+1;x<xmax;x++)
		{
		  for(int y=ymin+1;y<ymax;y++)
		    {
		      for(int k=0;k<9;k++)
			{
			  fout[IDX(x,y,k)] = w[k];
			}
		    }
		}
	      /*Swap populations*/
	      temp = fin;
	      fin = fout;
	      fout = temp;

	      /*Tick Timer*/

		
	    }
	  gettimeofday(&end,NULL);
	  cout << (end.tv_sec - start.tv_sec)*1e6 + (end.tv_usec - start.tv_usec) << endl;
	  
}





