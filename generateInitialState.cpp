#include <iostream>
#include <cstdlib>
#include <cmath>
#include<fstream>

using namespace std;

int generateInitialState(double *f, string path_to_folder, int N)
{

  srand(time(NULL));
  
  double *fprim;
  double sum = 0.0;
  
  int nbFiles = 5;

  ifstream fileList, popFile; string path_to_file, fileName;
  int error = 0; double random_prefactor;
  
  path_to_file = path_to_folder + "popfiles_list.dat";
  fileList.open(path_to_file.c_str());
  if(fileList.is_open()){error = 0;}
  else{error = 1;}

  fprim = new double[N];
  for (int nn=0;nn<nbFiles;nn++)
    {
      random_prefactor = rand() / (double) RAND_MAX;
      sum += random_prefactor;
      fileList >> fileName;
      path_to_file = path_to_folder + fileName;
      cout << "Opening " << path_to_file << endl;
      popFile.open(path_to_file.c_str(), ios::binary);
      if(popFile.is_open())
	{
	  popFile.read((char*)&fprim[0], N*sizeof(double));
	  for(int i=0;i<N;i++)
	    {
	      f[i] = f[i] + random_prefactor*fprim[i];
	    }
	}
      else
	{
	  error = 1;
	}
      popFile.close();
    }

  delete[] fprim;

  for(int i =0 ; i<N ; i++)
    {
      f[i] = f[i]/sum;
    }

  return error;
}
