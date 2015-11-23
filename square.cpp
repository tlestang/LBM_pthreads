#ifndef __global__
#define __global__
#include "global.h"
#endif

void computeSquareBounceBack_TEST(double *fout, double *fin, int xmin, int xmax, int ymin, int ymax, int Dx, int Dy)
{
  for(int x=xmin+1;x<xmax;x++)
    {
      fout[IDX(x,ymin,7)] = fin[IDX(x,ymin,5)];
      fout[IDX(x,ymin,4)] = fin[IDX(x,ymin,2)];
      fout[IDX(x,ymin,8)] = fin[IDX(x,ymin,6)];
      /*
      fout[x][ymin][7] = fin[x][ymin][5];
      fout[x][ymin][4] = fin[x][ymin][2];
      fout[x][ymin][8] = fin[x][ymin][6];
      */

      fout[IDX(x,ymax,5)] = fin[IDX(x,ymax,7)];
      fout[IDX(x,ymax,2)] = fin[IDX(x,ymax,4)];
      fout[IDX(x,ymax,6)] = fin[IDX(x,ymax,8)];
      /*
      fout[x][ymax][5] = fin[x][ymax][7];
      fout[x][ymax][2] = fin[x][ymax][4];
      fout[x][ymax][6] = fin[x][ymax][8];
      */
    }
  for(int y=ymin+1;y<ymax;y++)
    {
      fout[IDX(xmin,y,6)] = fin[IDX(xmin,y,8)];
      fout[IDX(xmin,y,3)] = fin[IDX(xmin,y,1)];
      fout[IDX(xmin,y,7)] = fin[IDX(xmin,y,5)];
      /*
      fout[xmin][y][6] = fin[xmin][y][8];
      fout[xmin][y][3] = fin[xmin][y][1];
      fout[xmin][y][7] = fin[xmin][y][5];
      */
      fout[IDX(xmax,y,8)] = fin[IDX(xmax,y,6)];
      fout[IDX(xmax,y,1)] = fin[IDX(xmax,y,3)];
      fout[IDX(xmax,y,5)] = fin[IDX(xmax,y,7)];
      /*
      fout[xmax][y][8] = fin[xmax][y][6];
      fout[xmax][y][1] = fin[xmax][y][3];
      fout[xmax][y][5] = fin[xmax][y][7];
      */
    }
  fout[IDX(xmin,ymin,7)] = fin[IDX(xmin,ymin,5)];
  fout[IDX(xmin,ymax,6)] = fin[IDX(xmin,ymax,8)];
  fout[IDX(xmax,ymax,5)] = fin[IDX(xmax,ymax,7)];
  fout[IDX(xmax,ymin,8)] = fin[IDX(xmax,ymin,6)];
  /*
  fout[xmin][ymin][7] = fin[xmin][ymin][5];
  fout[xmin][ymax][6] = fin[xmin][ymax][8];
  fout[xmax][ymax][5] = fin[xmax][ymax][7];
  fout[xmax][ymin][8] = fin[xmax][ymin][6];
  */
}

      
