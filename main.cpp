#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <malloc.h>
#include <unistd.h>

#ifndef __global__
#define __global__
#include "global.h"
#endif


#include "initialize_lattice_arrays.h"
#include "streamCollCompute.h"
#include "boundaryConditions.h"

using namespace std;

int main()
{
  /*Parameters for LB simulation*/
  int nbOfTimeSteps, Lx, Ly;
  int facquVtk, facquRe, facquForce;
  double tau, beta;
  double *fin, *fout, *temp, *rho, *ux, *uy;
  double Ma;   //Mach number
  string folderName, inputPopsFileName;
  const double w[9]={4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};
  
  /*Reads input file*/
      ifstream input_file("input.datin");
      input_file >> nbOfTimeSteps;
      input_file >> Lx; Ly = Lx;
      input_file >> tau;
      input_file >> Ma;
      input_file >> folderName;
      input_file >> inputPopsFileName;
      input_file >> facquVtk;
      input_file >> facquRe;
      input_file >> facquForce;
      input_file.close();
      
  /* --- Compute or define other parameters --- */
      int Dy = 4*Ly + 1, Dx = 2*(Dy-1) + 1;
      int xmin = (Dx-1)/2; int xmax = xmin + Lx;
      int ymin = (Dy-1)/2 - Ly/2; int ymax = ymin + Ly;
      double cs = 1./sqrt(3); double rho0 = 1.0;
      double u0 = cs*cs*Ma;
      double nu = 1./3.*(tau-0.5);
      double omega = 1.0/tau;
      beta = 8*nu*u0/((Dy-1)/2)/((Dy-1)/2);
  
  /* --- | Create folder for storing data | ---  */
      string instru = "mkdir " + folderName;
      system(instru.c_str());
      instru = "mkdir " + folderName + "/vtk_fluid/";
      system(instru.c_str());
      
  /* --- | Create parameters file | --- */
      string openParamFile = folderName + "/parameters.datout";
      ofstream param;
      param.open(openParamFile.c_str());
      param << "Number of timesteps : " << nbOfTimeSteps << endl;
      param << "L : "  << Lx << endl;
      param << "Dx : " << Dx << endl;
      param << "Dy : " << Dy << endl;
      param << "tau : " << tau << endl;
      param << "beta : " << beta << endl;
      param.close();


  /* ---- | Allocate populations and fields | --- */

      fin = (double *) memalign(getpagesize(), Dx*Dy*9*sizeof(double));
      fout = (double *) memalign(getpagesize(), Dx*Dy*9*sizeof(double));
      rho = (double *) memalign(getpagesize(), Dx*Dy*sizeof(double));
      ux = (double *) memalign(getpagesize(), Dx*Dy*sizeof(double));
      uy = (double *) memalign(getpagesize(), Dx*Dy*sizeof(double));

   /* --- Initialize pops to equilibrium value --- */
      initializePopulations(fin, Dx, Dy);
      initializeFields(rho, ux, uy, Dx, Dy);
   
  /* --- START LBM ---*/

	  for (int lbTimeStepCount=0; lbTimeStepCount<nbOfTimeSteps;lbTimeStepCount++)
	    {
	      cout << "Streaming and collison" <<endl;
	      streamingAndCollision_POSIX(fout, fin, rho, ux, uy, beta, tau, Dx, Dy);
	      cout << "walls" << endl;
	      computeDomainNoSlipWalls_BB(fout, fin, Dx, Dy);
	      cout << "square" << endl;
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

	    }
}



