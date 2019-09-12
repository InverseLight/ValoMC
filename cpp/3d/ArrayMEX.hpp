#ifndef __ARRAYMEX_HPP__
#define __ARRAYMEX_HPP__

#include <typeinfo>
#include "Array.hpp"
#include "mex.h"



// Convert MATLABs mxArray inp into Array class outp
template <typename T> void Convert_mxArray(const mxArray *inp, Array<T> &outp){
  // Check that the mxArray type is supported & matches T
  mxClassID id = mxGetClassID(inp);
  int idok = 0;
  if( (id == mxCHAR_CLASS) && (typeid(T) == typeid(char)) ) idok++;
  else if( (id == mxDOUBLE_CLASS) && (typeid(T) == typeid(double)) ) idok++;
  else if( (id == mxSINGLE_CLASS) && (typeid(T) == typeid(float)) ) idok++;
  else if( (id == mxINT8_CLASS) && (typeid(T) == typeid(char)) ) idok++;
  else if( (id == mxUINT8_CLASS) && (typeid(T) == typeid(unsigned char)) ) idok++;
  else if( (id == mxINT16_CLASS) && (typeid(T) == typeid(short)) ) idok++;
  else if( (id == mxUINT16_CLASS) && (typeid(T) == typeid(unsigned short)) ) idok++;
  else if( (id == mxINT32_CLASS) && (typeid(T) == typeid(int)) ) idok++;
  else if( (id == mxUINT32_CLASS) && (typeid(T) == typeid(unsigned int)) ) idok++;
  else if( (id == mxINT64_CLASS) && (typeid(T) == typeid(int_fast64_t)) ) idok++;
  else if( (id == mxUINT64_CLASS) && (typeid(T) == typeid(uint_fast64_t)) ) idok++;

  if(!idok) mexErrMsgTxt("mxArray type error! Incorrect type received!\n");

  // Check that the mxArray is of supported size
  mwSize ndim = mxGetNumberOfDimensions(inp);
  if(ndim != 2) mexErrMsgTxt("Only 2D-mxArrays are supported!\n");

  // Parse data to Array
  outp.data = (T *) mxGetData(inp);
  outp.Nx = mxGetM(inp);
  outp.Ny = mxGetN(inp);
  outp.Nz = 1;
  outp.Nxy = outp.Nx * outp.Ny;
  outp.N = outp.Nx * outp.Ny * outp.Nz;
  outp.IsRef = 1;
  outp.rank = 2;
}

// Convert MATLABs complex mxArray inp into Array class outp and ioutp which correspond to real and imaginary data
template <typename T> void Convert_mxArray(const mxArray *inp, Array<T> &outp, Array<T> &ioutp){
  // Check that the mxArray type is supported & matches T
  mxClassID id = mxGetClassID(inp);
  int idok = 0;
  if( (id == mxCHAR_CLASS) && (typeid(T) == typeid(char)) ) idok++;
  else if( (id == mxDOUBLE_CLASS) && (typeid(T) == typeid(double)) ) idok++;
  else if( (id == mxSINGLE_CLASS) && (typeid(T) == typeid(float)) ) idok++;
  else if( (id == mxINT8_CLASS) && (typeid(T) == typeid(char)) ) idok++;
  else if( (id == mxUINT8_CLASS) && (typeid(T) == typeid(unsigned char)) ) idok++;
  else if( (id == mxINT16_CLASS) && (typeid(T) == typeid(short)) ) idok++;
  else if( (id == mxUINT16_CLASS) && (typeid(T) == typeid(unsigned short)) ) idok++;
  else if( (id == mxINT32_CLASS) && (typeid(T) == typeid(int)) ) idok++;
  else if( (id == mxUINT32_CLASS) && (typeid(T) == typeid(unsigned int)) ) idok++;
  else if( (id == mxINT64_CLASS) && (typeid(T) == typeid(int_fast64_t)) ) idok++;
  else if( (id == mxUINT64_CLASS) && (typeid(T) == typeid(uint_fast64_t)) ) idok++;
  if(!idok) mexErrMsgTxt("mxArray type error! Incorrect type received!\n");

  // Check that the mxArray is of supported size
  mwSize ndim = mxGetNumberOfDimensions(inp);
  if(ndim != 2) mexErrMsgTxt("Only 2D-mxArrays are supported!\n");

  // Parse data to Arrays
  outp.data = (T *) mxGetData(inp);
  ioutp.data = (T *) mxGetImagData(inp);
  if(outp.data != NULL){
    outp.Nx = ioutp.Nx = mxGetM(inp);
    outp.Ny = ioutp.Ny = mxGetN(inp);
    outp.Nz = ioutp.Nz = 1;
    outp.Nxy = ioutp.Nxy = outp.Nx * outp.Ny;
    outp.N = ioutp.N = outp.Nx * outp.Ny * outp.Nz;
    outp.IsRef = ioutp.IsRef = 1;
    outp.rank = ioutp.rank = 2;
  }
  else{
    outp.Nx = outp.Ny = outp.Nz = outp.Nxy = outp.N = outp.IsRef = outp.rank = 0;
    ioutp.Nx = ioutp.Ny = ioutp.Nz = ioutp.Nxy = ioutp.N = ioutp.IsRef = ioutp.rank = 0;
  }
}


// Create output array of size Nx, Ny content of which can be modified through arr
template <typename T> void Convert_mxArray(mxArray **mx, Array<T> &arr, long Nx, long Ny){
  mxClassID id;
  if( typeid(T) == typeid(char) ) id = mxCHAR_CLASS;
  else if( typeid(T) == typeid(double) ) id = mxDOUBLE_CLASS;
  else if( typeid(T) == typeid(float) ) id = mxSINGLE_CLASS;
  else if( typeid(T) == typeid(char) ) id = mxINT8_CLASS;
  else if( typeid(T) == typeid(unsigned char) ) id = mxUINT8_CLASS;
  else if( typeid(T) == typeid(short) ) id = mxINT16_CLASS;
  else if( typeid(T) == typeid(unsigned short) ) id = mxUINT16_CLASS;
  else if( typeid(T) == typeid(int) ) id = mxINT32_CLASS;
  else if( typeid(T) == typeid(unsigned int) ) id = mxUINT32_CLASS;
//  else if( typeid(T) == typeid(long) ) id = mxINT64_CLASS;
//  else if( typeid(T) == typeid(unsigned long) ) id = mxUINT64_CLASS;
  else if( typeid(T) == typeid(int_fast64_t) ) id = mxINT64_CLASS;
  else if( typeid(T) == typeid(uint_fast64_t) ) id = mxUINT64_CLASS;
//  else if( (id == mxINT64_CLASS) && (typeid(T) == typeid(int_fast64_t)) ) idok++;
//  else if( (id == mxUINT64_CLASS) && (typeid(T) == typeid(uint_fast64_t)) ) idok++
  else mexErrMsgTxt("Trying to initialize mxArray with unsupported type!\n");

  // Create array
  mwSize dims[2] = { (mwSize) Nx, (mwSize) Ny };
  *mx = mxCreateNumericArray(2, const_cast <const mwSize *> (dims), id, mxREAL);

  // Parse data to Array
  arr.data = (T *) mxGetData(*mx);
  arr.Nx = mxGetM(*mx);
  arr.Ny = mxGetN(*mx);
  arr.Nz = 1;
  arr.Nxy = arr.Nx * arr.Ny;
  arr.N = arr.Nx * arr.Ny * arr.Nz;
  arr.IsRef = 1;
  arr.rank = 2;
}




// Create complex output array of size Nx, Ny content of which can be modified through arr and iarr
template <typename T> void Convert_mxArray(mxArray **mx, Array<T> &arr, Array<T> &iarr, long Nx, long Ny){
  mxClassID id;
  if( typeid(T) == typeid(char) ) id = mxCHAR_CLASS;
  else if( typeid(T) == typeid(double) ) id = mxDOUBLE_CLASS;
  else if( typeid(T) == typeid(float) ) id = mxSINGLE_CLASS;
  else if( typeid(T) == typeid(char) ) id = mxINT8_CLASS;
  else if( typeid(T) == typeid(unsigned char) ) id = mxUINT8_CLASS;
  else if( typeid(T) == typeid(short) ) id = mxINT16_CLASS;
  else if( typeid(T) == typeid(unsigned short) ) id = mxUINT16_CLASS;
  else if( typeid(T) == typeid(int) ) id = mxINT32_CLASS;
  else if( typeid(T) == typeid(unsigned int) ) id = mxUINT32_CLASS;
  else if( typeid(T) == typeid(long) ) id = mxINT64_CLASS;
  else if( typeid(T) == typeid(unsigned long) ) id = mxUINT64_CLASS;
  else mexErrMsgTxt("Trying to initialize mxArray with unsupported type!\n");

  // Create array
  mwSize dims[2] = { (mwSize) Nx, (mwSize) Ny };
  *mx = mxCreateNumericArray(2, const_cast <const mwSize *> (dims), id, mxCOMPLEX);

  // Parse data to Arrays
  arr.data = (T *) mxGetData(*mx);
  iarr.data = (T *) mxGetImagData(*mx);
  arr.Nx = iarr.Nx = mxGetM(*mx);
  arr.Ny = iarr.Ny = mxGetN(*mx);
  arr.Nz = iarr.Nz = 1;
  arr.Nxy = iarr.Nxy = arr.Nx * arr.Ny;
  arr.N = iarr.N = arr.Nx * arr.Ny * arr.Nz;
  arr.IsRef = iarr.IsRef = 1;
  arr.rank = iarr.rank = 2;
}


#endif
