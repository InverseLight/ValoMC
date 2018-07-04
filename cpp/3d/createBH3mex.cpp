
#include <string>
#define _USE_MATH_DEFINES
#define VALOMC_MEX
#include <cmath>
#include <limits>
#include <inttypes.h>
#include <string>

#include "mex.h"
#include "Array.hpp"
#include "ArrayMEX.hpp"
#include "MC3D.hpp"

void mexFunction(int nlhs, mxArray **plhs, int nrhs, const mxArray **prhs)
{
  
  if ((nrhs != 1) || (nlhs != 1))
  {
    mexErrMsgTxt("Syntax:\n [HN] = createBH3mex(H)\n");
  }
  Array<int_fast64_t> H;
  Array<double> r;
  Convert_mxArray(prhs[0], H);
  // Set parameters to MC
  MC3D MC;
  MC.H = H;
  int_fast64_t maxi = 0;
  for(int ii = 0; ii < MC.H.N; ii++) {if(maxi < MC.H[ii]) {maxi=MC.H[ii];} MC.H[ii]=MC.H[ii]-1;}
  r.resize(maxi);
  MC.r = r;
  MC.BuildNeighbourhoods();
  r.destroy();
  Array<int_fast64_t> HNo;
  int num_no_neighbor = 0;
  for(int ii = 0; ii < MC.HN.Nx; ii++) {
     if(MC.HN(ii, 0) == INT_FAST64_MAX) {
        num_no_neighbor++;
     }
     if(MC.HN(ii, 1) == INT_FAST64_MAX) {
        num_no_neighbor++;
     }
     if(MC.HN(ii, 2) == INT_FAST64_MAX) {
        num_no_neighbor++;
     }
     if(MC.HN(ii, 3) == INT_FAST64_MAX) {
        num_no_neighbor++;
     }
  }
  Convert_mxArray(&plhs[0], HNo, num_no_neighbor, 3);

  int bhi = 0;
  for(int ii = 0; ii < MC.HN.Nx; ii++) {
     if(MC.HN(ii, 0) == INT_FAST64_MAX) {
        HNo(bhi, 0) = MC.H(ii, 0)+1;
        HNo(bhi, 1) = MC.H(ii, 1)+1;
        HNo(bhi, 2) = MC.H(ii, 2)+1;
        bhi++;
     }
     if(MC.HN(ii, 1) == INT_FAST64_MAX) {
        HNo(bhi, 0) = MC.H(ii, 0)+1;
        HNo(bhi, 1) = MC.H(ii, 1)+1;
        HNo(bhi, 2) = MC.H(ii, 3)+1;
        bhi++;
     }
     if(MC.HN(ii, 2) == INT_FAST64_MAX) {
        HNo(bhi, 0) = MC.H(ii, 0)+1;
        HNo(bhi, 1) = MC.H(ii, 2)+1;
        HNo(bhi, 2) = MC.H(ii, 3)+1;
        bhi++;
     }
     if(MC.HN(ii, 3) == INT_FAST64_MAX) {
        HNo(bhi, 0) = MC.H(ii, 1)+1;
        HNo(bhi, 1) = MC.H(ii, 2)+1;
        HNo(bhi, 2) = MC.H(ii, 3)+1;
        bhi++;
     }
  }
}
