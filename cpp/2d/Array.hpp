#ifndef __ARRAY_HPP__
#define __ARRAY_HPP__
#include <inttypes.h>

// A simple 1, 2, 3D array class with fortran style indexing

template <class T> class Array{
public:
  Array(){
    IsRef = 0;
    rank = 0;
    Nx = Ny = Nz = N = 0;
    data = NULL;
  } 

  Array(Array<T> &ref){
    IsRef = 1;
    rank = ref.rank;
    Nx = ref.Nx; Ny = ref.Ny; Nz = ref.Nz; N = ref.N;
    data = ref.data;
  } 

  Array &operator=(const Array &ref){
    if(this != &ref){
      IsRef = 1;
      rank = ref.rank;
      Nx = ref.Nx; Ny = ref.Ny; Nz = ref.Nz; N = ref.N;
      data = ref.data;
    }
    return(*this);
  }

  T* resize(int_fast64_t _Nx, int_fast64_t _Ny, int_fast64_t _Nz){
    rank = 3;
    Nx = _Nx; Ny = _Ny; Nz = _Nz; N = Nx * Ny * Nz; Nxy = Nx * Ny;
    data = new T[N];
    return(data);
  }
  T* resize(int_fast64_t _Nx, int_fast64_t _Ny){
    rank = 2;
    Nx = _Nx; Ny = _Ny; Nz = 1; N = Nx * Ny * Nz; Nxy = Nx * Ny;
    data = new T[N];
    return(data);
  }
  T* resize(int_fast64_t _Nx){
    rank = 1;
    Nx = _Nx; Ny = 1; Nz = 1; N = Nx * Ny * Nz; Nxy = Nx * Ny;
    data = new T[N];
    return(data);
  }

  T& operator()(int_fast64_t ix, int_fast64_t iy, int_fast64_t iz){
    return( data[ ix + Nx * iy + Nxy * iz] );
  }
  T& operator()(int_fast64_t ix, int_fast64_t iy){
    return( data[ ix + Nx * iy] );
  }
  T& operator()(int_fast64_t ix){
    return( data[ ix ] );
  }
  T& operator[](int_fast64_t index){
    return( data[index] );
  }

  void destroy(){
    rank = Nx = Ny = Nz = Nxy = N = 0;
    if(!IsRef) if(data != NULL) delete[] data;
    data = NULL;
  }

  ~Array(){
    destroy();
  }

public:
  int IsRef;
  int_fast64_t rank, Nx, Ny, Nz, Nxy, N;
  T *data;
};




#endif
