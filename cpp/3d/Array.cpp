#include <stdio.h>

#include "Array.hpp"


// Print array
void PrintArrLin(Array<int> &arr){
  int ii;
  for(ii = 0; ii < arr.N; ii++){
    printf(" %d\n", arr[ii]);
  }
}

// Print 2D array
void PrintArr(Array<int> &arr){
  int ii, jj;
  for(ii = 0; ii < arr.Nx; ii++){
    for(jj = 0; jj < arr.Ny; jj++){
      printf(" %d", arr(ii, jj));
    }
    printf("\n");
  }
}



// Test the functionality of Array
int main(int argc, char **argv){

  Array<int> arr;
  int ii, jj, kk;

  arr.resize(7, 3);
  kk = 0;
  for(ii = 0; ii < 7; ii++)
    for(jj = 0; jj < 3; jj++) 
      arr(ii, jj) = kk++;

  printf("Original array (%ld x %ld x %ld = %ld elements):\n", arr.Nx, arr.Ny, arr.Nz, arr.N);
  PrintArr(arr);
  printf("\n");

  printf("Original array with linear print (%ld x %ld x %ld = %ld elements):\n", arr.Nx, arr.Ny, arr.Nz, arr.N);
  PrintArrLin(arr);
  printf("\n");

  
  Array<int> cparr(arr);
  printf("Copy array (%ld x %ld x %ld = %ld elements):\n", cparr.Nx, cparr.Ny, cparr.Nz, cparr.N);
  PrintArr(cparr);
  printf("\n"); 

  return(0);
}
