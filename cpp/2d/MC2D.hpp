#ifndef __MC2D_HPP__
#define __MC2D_HPP__

// define USE_OMP prior to including MC3D.hpp to utilize OpenMP
// define USE_MPI prior to including MC3D.hpp to utilize MPI
#define _USE_MATH_DEFINES
#include <math.h>

#include <ctime>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <inttypes.h>
#include <vector>

#ifndef INT_FAST64_MAX
#define INT_FAST64_MAX __INT_FAST64_MAX__
#endif

#include "Array.hpp"
#include "../Errors.hpp"

#include "mt_rng.hpp"

#ifdef USE_OMP
#include <omp.h>
#endif

#ifdef USE_MPI
#include "mpi.h"
#include "ArrayMPI.hpp"
#endif

#define abs(a) ((a) >= 0.0 ? (a) : (-(a)))
#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))

double eps = std::numeric_limits<double>::epsilon();

// Check if ray R + t D intersects with a segment (U, V)
int RayRayIntersects(double R[2], double D[2], double U[2], double V[2], double *t);

// Structure to hold information of a photon-packet
typedef struct _Photon
{
  double pos[2], dir[2];
  int_fast64_t curel, nextel, curface, nextface;
  double weight;
  double phase;
} Photon;

void PrintPhoton(Photon *p)
{
  printf(" curel: %ld -> %ld, curface: %ld -> %ld\n", p->curel, p->nextel, p->curface, p->nextface);
  printf(" %e %e\n\n", p->weight, p->phase);
}

// Class for 2D Optical Monte Carlo
class MC2D
{
public:
  // Constructor, destructor & Assingment operator
  MC2D();
  ~MC2D();
  MC2D &operator=(const MC2D &ref);

  // Random number generation related
  void InitRand();
  double UnifClosed();   // [0, 1]
  double UnifOpen();     // ]0 ,1[
  double UnifHalfDown(); // ]0, 1]
  double UnifHalfUp();   // [0, 1[

  // Area of triangle el
  double ElementArea(int_fast64_t el);
  // Length of boundary element ib, or triangle el edge f
  double ElementLength(int_fast64_t ib);
  double ElementLength(int_fast64_t el, int_fast64_t f);

  // Normal of boundary element ib, or triangle el and edge f
  void Normal(int_fast64_t ib, double *normal);
  void Normal(int_fast64_t el, int_fast64_t f, double *normal);

  // Initializes the MC2D after all the problem definition parameters have been given
  // Ie. constructs missing parameters
  void Init();

  // Build neighbourhood for volumetric elements
  void BuildNeighbourhoods();

  // Build light sources based on boundary conditions
  void BuildLightSource();

  // Given position, direction, current element & current face photon is on,
  // will return nonzero if the photon will hit a face in element, and also the distance to the point of intersection
  int WhichFace(Photon *phot, double *dist);

  // Create, scatter, mirror & propagate a photon
  void CreatePhoton(Photon *phot);
  void ScatterPhoton(Photon *phot);
  void MirrorPhoton(Photon *phot, int_fast64_t ib);
  void MirrorPhoton(Photon *phot, int_fast64_t el, int_fast64_t f);
  int FresnelPhoton(Photon *phot);
  void PropagatePhoton(Photon *phot);
  inline void search_neighbor(std::vector<int_fast64_t> &neighborlist, int_fast64_t element);
  // Perform MonteCarlo computation
  void MonteCarlo(bool (*progress)(double) = NULL);

  // [AL] Check if the arrays seem valid
  void ErrorChecks();

public:
  // Geometry
  Array<int_fast64_t> H, HN, BH; // Topology, Neigbourhood, Boundary
  Array<double> r;               // Grid nodes

  // Material parameters for each element
  Array<double> mua, mus, g, n; // Absorption, scattering & scattering inhomogeneity, index of refraction
  Array<double> k, g2;          // Wave number = omega / c * n, square of g

  // Boundary definitions for each Boundary element
  Array<char> BCType;

  // BCType   =   a  -  Absorbing boundary condition
  //              l  -  Laser light source boundary condition
  //              i  -  Isotropic light source
  //              c  -  Cosinic light source
  //              g  -  Gaussian light source

  Array<char> BCLightDirectionType; // [AL] Sets the direction of the light source

  // BCLightDirectionType =   n - No direction provided, use surface normal
  //                          a - Absolute
  //                          r - Relative to surface normal (\vec i is aligned along surface
  //                                                          \vec j is aligned along normal)

  Array<double> BCIntensity; // [AL] Sets the intensity of the light source
  Array<double> BCLNormal;
  Array<double> BCn;
  Array<double> GaussianSigma;

  // BCLNormal = Main direction the light propagates towards
  // BCn = Index of refraction on the outside of the boundary
  // GaussianSigma = Standard deviation describing Gaussian light sources

  // Frequency and angular frequency of amplitude modulation of the light source
  double f, omega, phase0; // [AL] Phase0

  // Number of photons to compute
  int_fast64_t Nphoton;

  // Speed of light (mm / ps)
  double c0;

  // Calculatable parameters
  Array<double> ER, EI;     // Absorbed power density in the volumetric elements (real & imaginary)
  Array<double> EBR, EBI;   // Absorbed power density in the boundary elements (real & imaginary)
  Array<double> DEBR, DEBI; // Absorbed power density in the boundary elements (real & imaginary) weighted by the dot product
                            // of the direction of the photon packets and the boundary normal

  // Light source likelyhood & creation variables
  Array<int_fast64_t> LightSources;
  Array<int_fast64_t> LightSourcesMother;
  Array<double> LightSourcesCDF;

  // Model parameters
  double weight0, chance; // Weight when to commence the roulette & the chance of revitalizing the photon
  long loss;              // Number of photons lost

  // Thread-safe random number generator
  mt_rng rng;

  // Seed id
  unsigned long seed;

  // OpeMP & MPI related variables
  int threadcount;
  int rank, nodecount;
  int totalthreads;
};

// Constuctor, set some default values for Monte Carlo
MC2D::MC2D()
{
  c0 = 2.99792458e11;

  Nphoton = 1;
  f = omega = 0.0;
  weight0 = 0.001;
  chance = 0.1;
  loss = 0;

  seed = 5489UL;

  rank = 0;
  threadcount = nodecount = totalthreads = 1;
}

// Nothing needs to be done, Arrays will deallocate themselves when it's time
MC2D::~MC2D()
{
}

// Assingment operator:
//  This will copy references to geometrical and parametric variables
//  Only new variables in the left hand side will be ER/EI, EBR/EBI, VER/VEI
MC2D &MC2D::operator=(const MC2D &ref)
{
  if (this != &ref)
  {
    H = ref.H;
    HN = ref.HN;
    BH = ref.BH;

    r = ref.r;
    mua = ref.mua;
    mus = ref.mus;
    g = ref.g;
    n = ref.n;
    k = ref.k;
    g2 = ref.g2;

    BCType = ref.BCType;
    BCLNormal = ref.BCLNormal;
    BCIntensity = ref.BCIntensity;
    BCLightDirectionType = ref.BCLightDirectionType;

    BCn = ref.BCn;
    GaussianSigma = ref.GaussianSigma;
    f = ref.f;
    phase0 = ref.phase0; // [AL]
    omega = ref.omega;
    Nphoton = ref.Nphoton;
    c0 = ref.c0;

    ER.resize(ref.ER.N);
    EI.resize(ref.EI.N);
    EBR.resize(ref.EBR.N);
    EBI.resize(ref.EBI.N);
    DEBR.resize(ref.DEBR.N); // [AL]
    DEBI.resize(ref.DEBI.N); // [AL]

    long ii;

    for (ii = 0; ii < ER.N; ii++)
      ER[ii] = EI[ii] = 0.0;
    for (ii = 0; ii < EBR.N; ii++)
      EBR[ii] = EBI[ii] = 0.0;
    for (ii = 0; ii < DEBR.N; ii++) // [AL]
      DEBR[ii] = DEBI[ii] = 0.0;    // [AL]

    LightSources = ref.LightSources;
    LightSourcesMother = ref.LightSourcesMother;
    LightSourcesCDF = ref.LightSourcesCDF;

    weight0 = ref.weight0;
    chance = ref.chance;
    loss = ref.loss;

    threadcount = ref.threadcount;
    rank = ref.rank;
    nodecount = ref.nodecount;
    totalthreads = ref.totalthreads;

    seed = ref.seed;
    InitRand();
  }

  return (*this);
}

// Initialize random number generator
void MC2D::InitRand()
{
  rng.Seed(seed);
}

// Draw random number on [0, 1]
double MC2D::UnifClosed()
{
  return (rng.drand_closed());
}

// Draw random number on ]0, 1[
double MC2D::UnifOpen()
{
  return (rng.drand_open());
}

// Draw random number ]0, 1]
double MC2D::UnifHalfDown()
{
  return (rng.drand_open_down());
}

// Draw random number [0, 1[
double MC2D::UnifHalfUp()
{
  return (rng.drand_open_up());
}

// Area of triangle el
double MC2D::ElementArea(int_fast64_t el)
{
  double ax = r(H(el, 1), 0) - r(H(el, 0), 0);
  double ay = r(H(el, 1), 1) - r(H(el, 0), 1);
  double bx = r(H(el, 2), 0) - r(H(el, 0), 0);
  double by = r(H(el, 2), 1) - r(H(el, 0), 1);
  return (0.5 * abs(ax * by - ay * bx));
}

// Length of boundary triangle ib
double MC2D::ElementLength(int_fast64_t ib)
{
  double len;
  len = sqrt(pow(r(BH(ib, 1), 0) - r(BH(ib, 0), 0), 2) + pow(r(BH(ib, 1), 1) - r(BH(ib, 0), 1), 2));
  return (len);
}

// Length of edge f of triangle el
double MC2D::ElementLength(int_fast64_t el, int_fast64_t f)
{
  double len;
  int i0, i1;
  if (f == 0)
  {
    i0 = 0;
    i1 = 1;
  }
  else if (f == 1)
  {
    i0 = 1;
    i1 = 2;
  }
  else if (f == 2)
  {
    i0 = 2;
    i1 = 0;
  }
  else
    return (0.0);
  len = sqrt(pow(r(H(el, i1), 0) - r(H(el, i0), 0), 2) + pow(r(H(el, i1), 1) - r(H(el, i0), 1), 2));
  return (len);
}

// Normal of boundary element ib
void MC2D::Normal(int_fast64_t ib, double *normal)
{
  double dx, dy, norm;
  dx = r(BH(ib, 1), 0) - r(BH(ib, 0), 0);
  dy = r(BH(ib, 1), 1) - r(BH(ib, 0), 1);
  norm = sqrt(dx * dx + dy * dy);
  normal[0] = dy / norm;
  normal[1] = -dx / norm;
}

// Normal of edge f on triangle el
void MC2D::Normal(int_fast64_t el, int_fast64_t f, double *normal)
{
  int i0, i1;
  if (f == 0)
  {
    i0 = 0;
    i1 = 1;
  }
  else if (f == 1)
  {
    i0 = 1;
    i1 = 2;
  }
  else if (f == 2)
  {
    i0 = 2;
    i1 = 0;
  }
  else
  {
    normal[0] = normal[1] = 0.0;
    return;
  }
  double dx, dy, norm;
  dx = r(H(el, i1), 0) - r(H(el, i0), 0);
  dy = r(H(el, i1), 1) - r(H(el, i0), 1);
  norm = sqrt(dx * dx + dy * dy);
  normal[0] = dy / norm;
  normal[1] = -dx / norm;
}

// Perform errorchecking and throw an error
void MC2D::ErrorChecks()
{
  /* SANITY CHECKS */
  // Check that
  // row size of H and g are equal
  // row size of H and mus are equal
  // row size of H and mua are equal
  // row size of H and n are equal
  // row size of BH and BCn are equal
  // row size of BH and BCType are equal
  // if given,  row size of BH and BCLNormal are equal
  // row size of BH and BCLightDirectionType are equal
  // H contains an index that cannot be found in r
  // BH contains an index that cannot be found in r

  if (g.Nx != H.Nx)
  {
    throw SIZE_MISMATCH_G;
    return;
  }

  // row size of H and mua are equal
  if (mua.Nx != H.Nx)
  {
    throw SIZE_MISMATCH_MUA;
    return;
  }

  // row size of H and mus are equal
  if (mus.Nx != H.Nx)
  {
    throw SIZE_MISMATCH_MUS;
    return;
  }

  // row size of H and n are equal
  if (n.Nx != H.Nx)
  {
    throw SIZE_MISMATCH_N;
    return;
  }

  // Sanity checks for BCn
  if (!BCn.N)
  {
    throw MISSING_BCN;
    return;
  }

  //row size of BH and BCn are equal
  if (BCn.Nx != BH.Nx)
  {
    throw SIZE_MISMATCH_BCN;
    return;
  }

  //row size of BH and BCType are equal
  if (BCType.Nx != BH.Nx)
  {
    throw SIZE_MISMATCH_BCTYPE;
    return;
  }

  // Make sure a light source exists
  bool has_lightsource = false;
  bool has_gaussian_lightsource = false;
  for (int ii = 0; ii < BCType.N; ii++)
  {
    if (BCType[ii] != 'a')
      has_lightsource = true;
    if (BCType[ii] == 'g')
      has_gaussian_lightsource = true;
  }

  if (!has_lightsource)
    throw NO_LIGHTSOURCE;

  if (has_gaussian_lightsource)
  {
    if (GaussianSigma.N != BH.Nx)
    {
      throw NO_GAUSSIAN_SIGMA;
    }
  }
  // Sanity checks for BCLNormal
  if (BCLNormal.N)
  {
    if (BCLNormal.Nx != BH.Nx)
    {
      throw SIZE_MISMATCH_LIGHT_DIRECTION;
      return;
    }
  }

  if (BCLightDirectionType.Nx)
  {
    if (BCLightDirectionType.Nx != BH.Nx)
    {
      throw SIZE_MISMATCH_LIGHT_DIRECTION_TYPE;
      return;
    }
  }

  // H contains an index that cannot be found in r

  for (int_fast64_t ii = 0; ii < H.N; ii++)
  {
    if (H[ii] < 0 || H[ii] >= r.Nx)
    {
      printf("Attempting to write entry %i to neighborhood (value %i), while it's size is %i \n", ii, H[ii], r.Nx);
      throw INCONSISTENT_H;
      return;
    }
  }
  // BH contains an index that cannot be found in r

  for (int_fast64_t ii = 0; ii < BH.N; ii++)
  {
    if (BH[ii] < 0 || BH[ii] >= r.Nx)
    {
      printf("Attempting to write entry %i to neighborhood (value %i), while it's size is %i \n", ii, BH[ii], r.Nx);
      throw INCONSISTENT_BH;
      return;
    }
  }
}

// Initialize Monte Carlo after geometry & material parameters have been assigned
// Under MPI also communicates relevant parameters to other computers and initializes
// mersenne twister with consequetive seed numbers
void MC2D::Init()
{
#ifdef USE_OMP
  threadcount = omp_get_max_threads();
#endif

#ifdef USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nodecount);

  MPI_Bcast(&seed, 1, MPI_LONG, 0, MPI_COMM_WORLD);
  seed += rank;

#ifdef USE_OMP
  MPI_Allreduce(&threadcount, &totalthreads, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif

  DistributeArray(H);
  DistributeArray(HN);
  DistributeArray(BH);
  DistributeArray(r);
  DistributeArray(mua);
  DistributeArray(mus);
  DistributeArray(g);
  DistributeArray(n);
  DistributeArray(BCType);
  DistributeArray(BCLNormal);
  DistributeArray(BCn);
  MPI_Bcast(&f, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&Nphoton, 1, MPI_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(&weight0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&chance, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  Nphoton /= nodecount;
#endif

  omega = 2.0 * M_PI * f;

  // Build neigborhood topology if required
  BuildNeighbourhoods();

  // Reserver memory for output variables
  int ii;
  ER.resize(H.Nx);
  EI.resize(H.Nx);
  for (ii = 0; ii < H.Nx; ii++)
    ER[ii] = EI[ii] = 0.0;
  EBR.resize(BH.Nx);
  EBI.resize(BH.Nx);

  DEBR.resize(BH.Nx); //[AL]
  DEBI.resize(BH.Nx); //[AL]

  for (ii = 0; ii < BH.Nx; ii++)
    EBR[ii] = EBI[ii] = 0.0;

  for (ii = 0; ii < BH.Nx; ii++) //[AL]
    DEBR[ii] = DEBI[ii] = 0.0;   //[AL]

  // Initialize BCIntensity to one if not given
  if (!BCIntensity.N)
  {
    BCIntensity.resize(BCType.N);
    int ii;
    for (ii = 0; ii < BCIntensity.N; ii++)
      BCIntensity[ii] = 1.0;
  }

  // [AL] Initialize BCLightDirectionType if not given
  if (!BCLightDirectionType.N)
  {
    int ii;
    BCLightDirectionType.resize(BCType.N);
    for (ii = 0; ii < BCLightDirectionType.N; ii++)
    {
      BCLightDirectionType[ii] = 'n';
    }
  }

  // Build lightsources
  BuildLightSource();

  // Compute wavenumber and g squared to speed computations
  k.resize(n.N);
  for (ii = 0; ii < n.N; ii++)
    k[ii] = omega / c0 * n[ii];

  g2.resize(g.N);

  for (ii = 0; ii < g.N; ii++)
    g2[ii] = pow(g[ii], 2);

  InitRand();

  // Normalize BCLNormal if given
  if (BCLNormal.N)
  {
    int ii;
    double norm;
    for (ii = 0; ii < BCLNormal.Nx; ii++)
    {
      norm = sqrt(pow(BCLNormal(ii, 0), 2) + pow(BCLNormal(ii, 1), 2));
      BCLNormal(ii, 0) /= norm;
      BCLNormal(ii, 1) /= norm;
    }
  }
  double n[2], t[2], norm;
  // [AL] Change BCLNormal coordinates to relative if needed
  for (ii = 0; ii < BCLightDirectionType.N; ii++)
  {
    Normal(ii, &n[0]);
    if (BCLightDirectionType[ii] == 'r')
    {
      t[0] = r(BH(ii, 1), 0) - r(BH(ii, 0), 0);
      t[1] = r(BH(ii, 1), 1) - r(BH(ii, 0), 1);
      norm = sqrt(t[0] * t[0] + t[1] * t[1]);
      t[0] /= norm;
      t[1] /= norm;
      double old_bclx = BCLNormal(ii, 0);
      double old_bcly = BCLNormal(ii, 1);
      BCLNormal(ii, 0) = old_bclx * t[0] + old_bcly * n[0];
      BCLNormal(ii, 1) = old_bclx * t[1] + old_bcly * n[1];
    }
  }

  return;
}

void MC2D::search_neighbor(std::vector<int_fast64_t> &neighborlist, int_fast64_t element)
{
#define NEW_NHOOD
#ifdef NEW_NHOOD
  for (int i = 0; i < neighborlist.size(); i++)
  {
    int_fast64_t neighbor = neighborlist[i];
    // Test for edge 0, 1, 2
    if (neighbor != element)
    {
      if (neighbor >= 0)
      {
      // If two vertices of a neighboring triangle
      // belongs to the current element, it must be a neighbor
#define ISPART(a) (H(element, a) == H(neighbor, 0) || H(element, a) == H(neighbor, 1) || H(element, a) == H(neighbor, 2))
        if (ISPART(0) && ISPART(1))
        {
          if (HN(element, 0) != INT_FAST64_MAX && HN(element, 0) != neighbor)
          {
            throw INCONSISTENT_MESH;
          }
          HN(element, 0) = neighbor;
        }
        if (ISPART(1) && ISPART(2))
        {
          if (HN(element, 1) != INT_FAST64_MAX && HN(element, 1) != neighbor)
          {
            throw INCONSISTENT_MESH;
          }
          HN(element, 1) = neighbor;
        }
        if (ISPART(2) && ISPART(0))
        {
          if (HN(element, 2) != INT_FAST64_MAX && HN(element, 2) != neighbor)
          {
            throw INCONSISTENT_MESH;
          }
          HN(element, 2) = neighbor;
        }
#undef ISPART
      }
      else
      {
      // the neighbor is a boundary element
#define ISPART2(a) (H(element, a) == BH(-(neighbor + 1), 0) || H(element, a) == BH(-(neighbor + 1), 1))
        if (ISPART2(0) && ISPART2(1))
        {
          // Something is already written to the neighborhood matrix
          // This should never happen
          if (HN(element, 0) != INT_FAST64_MAX && HN(element, 0) != neighbor)
          {
            throw INCONSISTENT_MESH;
          }
          HN(element, 0) = neighbor;
        }
        if (ISPART2(1) && ISPART2(2))
        {
          if (HN(element, 1) != INT_FAST64_MAX && HN(element, 1) != neighbor)
          {
            throw INCONSISTENT_MESH;
          }
          HN(element, 1) = neighbor;
        }
        if (ISPART2(2) && ISPART2(0))
        {
          if (HN(element, 2) != INT_FAST64_MAX && HN(element, 2) != neighbor)
          {
            throw INCONSISTENT_MESH;
          }
          HN(element, 2) = neighbor;
        }
#undef ISPART2
      }
    }
  }
#endif
}

// Build neigbourhood HN for topology H
void MC2D::BuildNeighbourhoods()
{
#ifdef NEW_NHOOD
  if (HN.N != H.N)
  {
    HN.resize(H.Nx, 3);
    std::vector<std::vector<int_fast64_t>> parent;
    parent.resize(r.Nx);
    // Build a vector (parent)
    // that contains all tetrahedrons shared by a vertex
    for (int i = 0; i < H.Nx; i++)
    {
      parent[H(i, 0)].push_back(i);
      parent[H(i, 1)].push_back(i);
      parent[H(i, 2)].push_back(i);
    }
    for (int i = 0; i < BH.Nx; i++)
    {
      parent[BH(i, 0)].push_back(-1 - i);
      parent[BH(i, 1)].push_back(-1 - i);
    }
    // Parent now contains all triangles and boundary elements for
    // each vertex in a vector

    int istart = 0, iend = H.Nx;
    for (int ii = istart; ii < iend; ii++)
    {

      HN(ii, 0) = INT_FAST64_MAX;
      HN(ii, 1) = INT_FAST64_MAX;
      HN(ii, 2) = INT_FAST64_MAX;

      // Search neighbors from the neighboring triangles instead of all triangles
      // Neighboring triangles are all triangles that share a vertex with the current triangle

      // give the neighborlist of each corner vertex as an input

      // neighbors of the first vertex
      search_neighbor(parent[H(ii, 0)], ii);
      // neighbors of the second vertex
      search_neighbor(parent[H(ii, 1)], ii);
      // neighbors of the third vertex
      search_neighbor(parent[H(ii, 2)], ii);

      if (HN(ii, 0) == INT_FAST64_MAX || HN(ii, 1) == INT_FAST64_MAX || HN(ii, 2) == INT_FAST64_MAX)
      {
        throw MISSING_BOUNDARY;
      }
    }
  }

#else
  if (HN.N != H.N)
  {
    HN.resize(H.Nx, 3);

    int_fast64_t istart = 0, iend = H.Nx;
    int_fast64_t ii, jj;
#ifdef USE_MPI
    istart = rank * H.Nx / nodecount;
    iend = (rank + 1) * H.Nx / nodecount;
    for (ii = 0; ii < HN.N; ii++)
      HN[ii] = 0;
#endif

#ifdef USE_OMP
#pragma omp parallel for private(ii, jj)
#endif
    for (ii = istart; ii < iend; ii++)
    {
      HN(ii, 0) = HN(ii, 1) = HN(ii, 2) = -1;

      for (jj = 0; jj < H.Nx; jj++)
      {
        if (ii == jj)
          continue;

#define ISIN(a, b) (H(ii, a) == H(jj, b))
#define ISPART(a) (ISIN(a, 0) || ISIN(a, 1) || ISIN(a, 2))
        // Test for Edge 0, 1, 2
        if (ISPART(0) && ISPART(1))
          HN(ii, 0) = jj;
        if (ISPART(1) && ISPART(2))
          HN(ii, 1) = jj;
        if (ISPART(2) && ISPART(0))
          HN(ii, 2) = jj;
#undef ISIN
#undef ISPART
      }

      // Fill -1's with appropriate boundary element from BH
      if ((HN(ii, 0) == -1) || (HN(ii, 1) == -1) || (HN(ii, 2) == -1))
      {
        for (jj = 0; jj < BH.Nx; jj++)
        {
#define ISIN2(a, b) (H(ii, a) == BH(jj, b))
#define ISPART2(a) (ISIN2(a, 0) || ISIN2(a, 1))
          if ((HN(ii, 0) == -1) && (ISPART2(0) && ISPART2(1)))
            HN(ii, 0) = -1 - jj;
          if ((HN(ii, 1) == -1) && (ISPART2(1) && ISPART2(2)))
            HN(ii, 1) = -1 - jj;
          if ((HN(ii, 2) == -1) && (ISPART2(2) && ISPART2(0)))
            HN(ii, 2) = -1 - jj;
#undef ISIN2
#undef ISPART2
        }
      }
    }

#ifdef USE_MPI
    // Distribute neigbourhood to other computers
    AllReduceArray(HN, MPI_SUM);
#endif
  }
#endif
#if 0 
  for(int i = 0; i < HN.Nx; i++) {
    printf("%i: %li %li %li\n", i, HN(i, 0), HN(i,1), HN(i,2));
  }
#endif
}

// Build light sources based on boundary conditions
//   LightSources will contain index to the boundary element in BH acting as lightsource
//   LightSourcesMother will contain index to volumetric element H for which BH is attached
//   LightSourcesCDF will be a cumulative/normalized sum of areas of all the lightsources,
//   this will ease randomizing the creation of photons
void MC2D::BuildLightSource()
{
  int_fast64_t ii, jj, kk, ib, NLightSource;
  if (LightSources.N != 0)
    return;
  // Loop over HN and determine the number of lightsources (boundary elements)
  NLightSource = 0;
  for (ii = 0; ii < HN.Nx; ii++)
  {
    for (jj = 0; jj < 3; jj++)
    {
      if (HN(ii, jj) <= -1)
      {
        ib = -1 - HN(ii, jj);
        if ((BCType[ib] == 'l') ||
            (BCType[ib] == 'i') ||
            (BCType[ib] == 'c') ||
            (BCType[ib] == 'g'))
          NLightSource++;
      }
    }
  }

  // Allocate space for light sources
  LightSources.resize(NLightSource);
  LightSourcesMother.resize(NLightSource);
  LightSourcesCDF.resize(NLightSource);

  // Compute length of each boundary element acting as a lightsource & assemble lightsource motherhood
  kk = 0;
  for (ii = 0; ii < HN.Nx; ii++)
  {
    for (jj = 0; jj < 3; jj++)
    {
      if (HN(ii, jj) <= -1)
      {
        ib = -1 - HN(ii, jj);
        if ((BCType[ib] == 'l') ||
            (BCType[ib] == 'i') ||
            (BCType[ib] == 'c') ||
            (BCType[ib] == 'g'))
        {
          LightSources[kk] = ib;
          LightSourcesMother[(int)kk] = ii;
          LightSourcesCDF[(int)kk] = ElementLength(ib);
          kk++;
        }
      }
    }
  }

  //  mexPrintf("%i\n", BCIntensity.N);

  for (ii = 0; ii < NLightSource; ii++)
    LightSourcesCDF[ii] *= BCIntensity[LightSources[ii]]; // [AL]

  // Compute cumsum of LightSourcesCDF and normalize -- Ie. form cumulated distribution function
  for (ii = 1; ii < NLightSource; ii++)
    LightSourcesCDF[ii] += LightSourcesCDF[ii - 1];
  for (ii = 0; ii < NLightSource; ii++)
    LightSourcesCDF[ii] /= LightSourcesCDF[NLightSource - 1];
}

// Determine which face a photon will exit a volumetric element from
int MC2D::WhichFace(Photon *phot, double *dist)
{
  // phot - photon under test
  // dist - distance the photon can travel before hitting the face
  //
  // Return values:
  // -1 if no hit
  // 0, 1, 2 for faces formed by (V0, V1), (V1, V2), (V2, V0) respectively

  double V0[2] = {r(H(phot->curel, 0), 0), r(H(phot->curel, 0), 1)};
  double V1[2] = {r(H(phot->curel, 1), 0), r(H(phot->curel, 1), 1)};
  double V2[2] = {r(H(phot->curel, 2), 0), r(H(phot->curel, 2), 1)};

  if (phot->curface != 0)
    if (RayRayIntersects(phot->pos, phot->dir, V0, V1, dist))
      if (*dist > 0.0)
      {
        phot->nextface = 0;
        phot->nextel = HN(phot->curel, phot->nextface);
        return (0);
      }

  if (phot->curface != 1)
    if (RayRayIntersects(phot->pos, phot->dir, V1, V2, dist))
      if (*dist > 0.0)
      {
        phot->nextface = 1;
        phot->nextel = HN(phot->curel, phot->nextface);
        return (1);
      }

  if (phot->curface != 2)
    if (RayRayIntersects(phot->pos, phot->dir, V2, V0, dist))
      if (*dist > 0.0)
      {
        phot->nextface = 2;
        phot->nextel = HN(phot->curel, phot->nextface);
        return (2);
      }

  return (-1);
}

// Create a new photon based on LightSources, LightSourcesMother and LighSourcesCDF
void MC2D::CreatePhoton(Photon *phot)
{
  double xi = UnifClosed();

  double n[2], t[2], norm;

  // Find boundary element index that will create this photon
  int_fast64_t ib;
  for (ib = 0; ib < LightSources.Nx; ib++)
    if (xi < LightSourcesCDF[ib])
      break;
  ib = ib > LightSources.Nx - 1 ? LightSources.Nx - 1 : ib;
  // Creator faces mother element
  phot->curel = LightSourcesMother[ib];
  // Creator face
  ib = LightSources[ib];
  if (-1 - HN(phot->curel, 0) == ib)
    phot->curface = 0;
  else if (-1 - HN(phot->curel, 1) == ib)
    phot->curface = 1;
  else if (-1 - HN(phot->curel, 2) == ib)
    phot->curface = 2;
  else
    phot->curface = -1;

  // Initial photon position uniformly distributed on the boundary element
  double w = UnifOpen();
  phot->pos[0] = r(BH(ib, 0), 0) + w * (r(BH(ib, 1), 0) - r(BH(ib, 0), 0));
  phot->pos[1] = r(BH(ib, 0), 1) + w * (r(BH(ib, 1), 1) - r(BH(ib, 0), 1));

  // Face normal
  Normal(ib, n);

  // Make sure that the normal points inside the element by checking dot product of normal & vector (phot->pos) to center of the element
  t[0] = (r(H(phot->curel, 0), 0) + r(H(phot->curel, 1), 0) + r(H(phot->curel, 2), 0)) / 3.0 - phot->pos[0];
  t[1] = (r(H(phot->curel, 0), 1) + r(H(phot->curel, 1), 1) + r(H(phot->curel, 2), 1)) / 3.0 - phot->pos[1];
  norm = sqrt(t[0] * t[0] + t[1] * t[1]);
  t[0] /= norm;
  t[1] /= norm;

  if (n[0] * t[0] + n[1] * t[1] < 0.0)
  {
    n[0] = -n[0];
    n[1] = -n[1];
  }

  // [AL] Changed the old if clause
  // If light direction type is is 'n' do not read anything from BCLNormal,
  // as the values might not be well defined
  if ((BCLightDirectionType[ib] == 'n'))
  {
    // No BCLNormal data -> Photons main propagation direction is towards inward normal
    // Select photons initial direction based on the boundary condition
    if ((BCType[ib] == 'l'))
    {
      // Laser -- Photons created in normal direction
      phot->dir[0] = n[0];
      phot->dir[1] = n[1];
    }
    else if ((BCType[ib] == 'i'))
    {
      // Isotropic -- Photons initial direction probality density is uniform on a circle
      double phi0, phi;
      phi0 = atan2(n[1], n[0]);
      phi = M_PI * (UnifOpen() - 0.5);
      phot->dir[0] = cos(phi0 + phi);
      phot->dir[1] = sin(phi0 + phi);
    }
    else if ((BCType[ib] == 'c'))
    {
      // Cosinic -- Directivity follows cosine pattern
      double phi0, phi;
      phi0 = atan2(n[1], n[0]);
      phi = asin(2.0 * UnifOpen() - 1.0);
      phot->dir[0] = cos(phi0 + phi);
      phot->dir[1] = sin(phi0 + phi);
    }
    else if ((BCType[ib] == 'g'))
    {
      // Gaussian -- Directivity pattern is Gaussian limited to incident angle ]-pi/2, pi/2[
      double phi0, /* u, v,*/ phi;
      phi0 = atan2(n[1], n[0]);
      // CDF of near-Gaussially distributed rng in range [-pi/2, pi/2] is
      //   (0.5 + 0.5 erf(x / sigma / sqrt(2)) - g) / (1 - 2 * h)
      // where
      //   h = 0.5 - 0.5 erf(M_PI / sigma / 2 / sqrt(2))
      double h = 0.5 - 0.5 * erf(M_PI / GaussianSigma[ib] / 2.0 / sqrt(2.0));
      double xi = UnifOpen();
      // Perform bisection search to invert CDF ie. find xi = CDF(phi)
      double a = -M_PI / 2, b = M_PI / 2, iter, CDF;
      for (iter = 0; iter < 16; iter++)
      {
        phi = (a + b) / 2.0;
        CDF = (0.5 + 0.5 * erf(phi / GaussianSigma[ib] / sqrt(2.0)) - h) / (1 - 2.0 * h);
        if (CDF - xi < 0.0)
          a = phi;
        else
          b = phi;
      }
      phot->dir[0] = cos(phi0 + phi);
      phot->dir[1] = sin(phi0 + phi);
    }
    //    mexPrintf("Photon without normal shot at %f %f\n", phot->dir[0],phot->dir[1]);
  }
  else
  {
    //   mexPrintf("Photon with a normal shot at %f %f\n", phot->dir[0],phot->dir[1]);

    // Select photons initial direction based on the boundary condition
    if ((BCType[ib] == 'l'))
    {
      // Laser -- Photons created in normal direction
      // Check normal type
      // direction is given in the mesh coordinates
      phot->dir[0] = BCLNormal(ib, 0);
      phot->dir[1] = BCLNormal(ib, 1);
      // printf("Launching photon %lf %lf\n", phot->dir[0], phot->dir[1]);
    }
    else if ((BCType[ib] == 'i'))
    {
      // Isotropic -- Photons initial direction probality density is uniform on a circle
      double phi0, phi;
      phi0 = atan2(BCLNormal(ib, 1), BCLNormal(ib, 0));
      phi = M_PI * (UnifOpen() - 0.5);
      phot->dir[0] = cos(phi0 + phi);
      phot->dir[1] = sin(phi0 + phi);
    }
    else if ((BCType[ib] == 'c'))
    {
      do
      {
        double phi0, phi;
        phi0 = atan2(BCLNormal(ib, 1), BCLNormal(ib, 0));
        phi = asin(2.0 * UnifOpen() - 1.0);
        phot->dir[0] = cos(phi0 + phi);
        phot->dir[1] = sin(phi0 + phi);

      } while (n[0] * phot->dir[0] + n[1] * phot->dir[1] <= 0.0);
    }
    else if ((BCType[ib] == 'g'))
    {
      // Gaussian -- Initial direction is draw from Gaussian distribution
      double phi0, /* u, v,*/ phi;
      phi0 = atan2(BCLNormal(ib, 1), BCLNormal(ib, 0));
      // CDF of near-Gaussially distributed rng in range [-pi/2, pi/2] is
      //   (0.5 + 0.5 erf(x / sigma / sqrt(2)) - g) / (1 - 2 * h)
      // where
      //   h = 0.5 - 0.5 erf(M_PI / sigma / 2 / sqrt(2))
      do
      {
        double h = 0.5 - 0.5 * erf(M_PI / GaussianSigma[ib] / 2.0 / sqrt(2.0));
        double xi = UnifOpen();
        // Perform bisection search to invert CDF ie. find xi = CDF(phi)
        double a = -M_PI / 2, b = M_PI / 2, iter, CDF;
        for (iter = 0; iter < 16; iter++)
        {
          phi = (a + b) / 2.0;
          CDF = (0.5 + 0.5 * erf(phi / GaussianSigma[ib] / sqrt(2.0)) - h) / (1 - 2.0 * h);
          if (CDF - xi < 0.0)
            a = phi;
          else
            b = phi;
        }
        phot->dir[0] = cos(phi0 + phi);
        phot->dir[1] = sin(phi0 + phi);
      } while (n[0] * phot->dir[0] + n[1] * phot->dir[1] <= 0.0);
    }
  }

  phot->nextel = -1;
  phot->nextface = -1;

  phot->weight = 1.0;
  phot->phase = phase0;
}

// Scatter a photon
void MC2D::ScatterPhoton(Photon *phot)
{
  if (g[phot->curel] != 0.0)
  {
    double xi = UnifHalfUp();
    double t = sqrt(pow(g2[phot->curel] + 1.0, 2) - 4.0 * g2[phot->curel]);
    double th = 2.0 * atan(tan(-xi * t * M_PI / (g2[phot->curel] - 1.0)) * (g2[phot->curel] - 2 * g[phot->curel] + 1.0) / t);
    double th0 = atan2(phot->dir[1], phot->dir[0]);
    phot->dir[0] = cos(th0 + th);
    phot->dir[1] = sin(th0 + th);
  }
  else
  {
    double th = 2.0 * M_PI * UnifHalfUp();
    phot->dir[0] = cos(th);
    phot->dir[1] = sin(th);
  }

  // This is to prevent RayRayIntersects from misbehaving after scattering event in the PropagatePhoton
  phot->curface = -1;
}

// Mirror photons propagation with respect to boundary element ib
void MC2D::MirrorPhoton(Photon *phot, int_fast64_t ib)
{
  double n[2], cdot;
  Normal(ib, n);
  cdot = n[0] * phot->dir[0] + n[1] * phot->dir[1];
  phot->dir[0] -= 2.0 * cdot * n[0];
  phot->dir[1] -= 2.0 * cdot * n[1];
}

// Mirror photon with respect to face f of element el
void MC2D::MirrorPhoton(Photon *phot, int_fast64_t el, int_fast64_t f)
{
  double n[2], cdot;
  Normal(el, f, n);
  cdot = n[0] * phot->dir[0] + n[1] * phot->dir[1];
  phot->dir[0] -= 2.0 * cdot * n[0];
  phot->dir[1] -= 2.0 * cdot * n[1];
}

// Fresnel transmission / reflection of a photon
int MC2D::FresnelPhoton(Photon *phot)
{
  // Likelyhood of reflection:
  //   R = 0.5 ( sin^2(theta_i - theta_t) / sin^2(theta_i + theta_t) + tan^2(theta_i - theta_t) / tan^2(theta_i + theta_t))
  //
  // For theta_i + theta_t < eps:
  //   R = ( (ni / nt - 1) / (ni / nt + 1) )^2
  // which is the limit as theta_i -> 0
  //
  // Transmission/Reflection of a incident direction di, through face with normal n:
  //   costhi = -n * di';
  //   costht = sqrt( 1 - (ni/nt)^2 * (1 - costhi^2) );
  //   dr = di + 2 * costhi * n;
  //   if(costhi > 0);
  //     dt = (ni/nt) * di + (ni/nt * costhi - costht) * n;
  //   else;
  //     dt = (ni/nt) * di + (ni/nt * costhi + costht) * n;
  //   end;

  // Normal of the tranmitting face
  double nor[3];
  Normal(phot->curel, phot->nextface, nor);
  double nipnt;

  // [AL] Made some minor modifications here
  // [AL] Do not change the direction of the photon if it escapes through the boundary
  if (phot->nextel < 0)
  {
    nipnt = n[phot->curel] / BCn[-1 - phot->nextel];

    double costhi = -(phot->dir[0] * nor[0] + phot->dir[1] * nor[1]);

    if (1.0 - pow(nipnt, 2) * (1.0 - pow(costhi, 2)) <= 0.0)
    {
      // Total reflection due to critical angle of Snell's law -- costht below would be complex
      phot->dir[0] += 2.0 * costhi * nor[0];
      phot->dir[1] += 2.0 * costhi * nor[1];
      phot->curface = phot->nextface;
      return (1);
    }
    double costht = sqrt(1.0 - pow(nipnt, 2) * (1.0 - pow(costhi, 2)));
    double thi;

    if (costhi > 0.0)
      thi = acos(costhi);
    else
      thi = acos(-costhi);

    double tht = acos(costht);
    double R;

    if (!(sin(thi + tht) > eps))
      R = pow((nipnt - 1.0) / (nipnt + 1.0), 2);
    else
      R = 0.5 * (pow(sin(thi - tht) / sin(thi + tht), 2) + pow(tan(thi - tht) / tan(thi + tht), 2));

    double xi = UnifClosed();

    if (xi <= R)
    {
      // Photon is reflected
      phot->dir[0] += 2.0 * costhi * nor[0];
      phot->dir[1] += 2.0 * costhi * nor[1];
      phot->curface = phot->nextface;

      return (1);
    }

    // Photon is transmitted
    return (0);
  }
  else
  {
    nipnt = n[phot->curel] / n[phot->nextel];

    double costhi = -(phot->dir[0] * nor[0] + phot->dir[1] * nor[1]);

    if (1.0 - pow(nipnt, 2) * (1.0 - pow(costhi, 2)) <= 0.0)
    {
      // Total reflection due to critical angle of Snell's law -- costht below would be complex
      phot->dir[0] += 2.0 * costhi * nor[0];
      phot->dir[1] += 2.0 * costhi * nor[1];
      phot->curface = phot->nextface;
      return (1);
    }
    double costht = sqrt(1.0 - pow(nipnt, 2) * (1.0 - pow(costhi, 2)));
    double thi;

    if (costhi > 0.0)
      thi = acos(costhi);
    else
      thi = acos(-costhi);

    double tht = acos(costht);
    double R;

    if (!(sin(thi + tht) > eps))
      R = pow((nipnt - 1.0) / (nipnt + 1.0), 2);
    else
      R = 0.5 * (pow(sin(thi - tht) / sin(thi + tht), 2) + pow(tan(thi - tht) / tan(thi + tht), 2));

    double xi = UnifClosed();

    if (xi <= R)
    {
      // Photon is reflected
      phot->dir[0] += 2.0 * costhi * nor[0];
      phot->dir[1] += 2.0 * costhi * nor[1];
      phot->curface = phot->nextface;

      return (1);
    }

    // Photon is transmitted - update propagation direction via Snell's law
    if (costhi > 0.0)
    {
      phot->dir[0] = nipnt * phot->dir[0] + (nipnt * costhi - costht) * nor[0];
      phot->dir[1] = nipnt * phot->dir[1] + (nipnt * costhi - costht) * nor[1];
    }
    else
    {
      phot->dir[0] = nipnt * phot->dir[0] + (nipnt * costhi + costht) * nor[0];
      phot->dir[1] = nipnt * phot->dir[1] + (nipnt * costhi + costht) * nor[1];
    }
    return (0);
  }
}

// Propagate a photon until it dies
void MC2D::PropagatePhoton(Photon *phot)
{
  double prop, dist, ds;
  int_fast64_t ib;

  // Propagate until the photon dies
  while (1)
  {
    // Draw the propagation distance
    prop = -log(UnifOpen()) / mus[phot->curel];

    // Propagate until the current propagation distance runs out (and a scattering will occur)
    while (1)
    {

      // Check through which face the photon will exit the current element
      if (WhichFace(phot, &dist) == -1)
      {
        loss++;
        return;
      }

      // Travel distance -- Either propagate to the boundary of the element, or to the end of the leap, whichever is closer
      ds = fmin(prop, dist);

      // Move photon
      phot->pos[0] += phot->dir[0] * ds;
      phot->pos[1] += phot->dir[1] * ds;

      // Upgrade element fluence
      if (omega <= 0.0)
      {
        // Unmodulated light
        if (mua[phot->curel] > 0.0)
        {
          ER[phot->curel] += (1.0 - exp(-mua[phot->curel] * ds)) * phot->weight;
        }
        else
        {
          ER[phot->curel] += phot->weight * ds;
        }
      }
      else
      {
        // Modulated light
        /*
 	    cls;
	    syms w0 mua k x ph0 s real;
	    % k = 0; ph0 = 0;
	    e = w0 * exp(-mua * x - j * (k * x + ph0));
	    g = int(e, x, 0, s);
	    syms a b real;
	    f = (a + i * b) / (mua + i * k);
	    % Change of element as photon passes it
	    pretty(simplify( real( g * (mua + i * k) ) ))
	    pretty(simplify( imag( g * (mua + i * k) ) ))
	    % Final divider / normalization
	    pretty( simplify( real(f) ) )
	    pretty( simplify( imag(f) ) )
	*/
        ER[phot->curel] += phot->weight * (cos(phot->phase) - cos(-phot->phase - k[phot->curel] * ds) * exp(-mua[phot->curel] * ds));
        EI[phot->curel] += phot->weight * (-sin(phot->phase) + sin(phot->phase + k[phot->curel] * ds) * exp(-mua[phot->curel] * ds));
        phot->phase += k[phot->curel] * ds;
      }

      // Upgrade photon weight
      phot->weight *= exp(-mua[phot->curel] * ds);

      // Photon has reached a situation where it has to be scattered
      prop -= ds;
      if (prop <= 0.0)
        break;

      // Otherwise the photon will continue to pass through the boundaries of the current element

      // Test for boundary conditions
      if (phot->nextel < 0)
      {
        // Boundary element index
        ib = -1 - phot->nextel;

        if ((BCType[ib] == 'm') || (BCType[ib] == 'L') || (BCType[ib] == 'I') || (BCType[ib] == 'C'))
        {
          // Mirror boundary condition -- Reflect the photon
          MirrorPhoton(phot, ib);
          phot->curface = phot->nextface;
          continue;
        }
        else
        {
          // Absorbing (a, l, i, c)
          // Check for mismatch between inner & outer index of refraction causes Fresnel transmission
          if (BCn[ib] > 0.0)
          {
            if (FresnelPhoton(phot))
            {
              continue;
            }
          }

          // Photon is transmitted

          double n[2], t[2];                      // [AL] Added normal calculation here for exitance
          Normal(phot->curel, phot->nextface, n); // [AL]

          // Make sure that the normal points outside of the element by checking dot product of
          // normal & vector (phot->pos) to center of the element

          t[0] = (r(H(phot->curel, 0), 0) + r(H(phot->curel, 1), 0) + r(H(phot->curel, 2), 0)) / 3.0 - phot->pos[0];
          t[1] = (r(H(phot->curel, 0), 1) + r(H(phot->curel, 1), 1) + r(H(phot->curel, 2), 1)) / 3.0 - phot->pos[1];

          if (n[0] * t[0] + n[1] * t[1] > 0.0)
          {
            n[0] = -n[0];
            n[1] = -n[1];
          }

          if (omega <= 0.0)
          {
            EBR[ib] += phot->weight;
            DEBR[ib] += phot->weight / (phot->dir[0] * n[0] + phot->dir[1] * n[1]);
          }
          else
          {
            EBR[ib] += phot->weight * cos(phot->phase);
            EBI[ib] -= phot->weight * sin(phot->phase);

            DEBR[ib] += phot->weight * cos(phot->phase) / (phot->dir[0] * n[0] + phot->dir[1] * n[1]); // AL
            DEBI[ib] -= phot->weight * sin(phot->phase) / (phot->dir[0] * n[0] + phot->dir[1] * n[1]); // AL
          }

          // Photon propagation will terminate
          return;
        }
      }

      // Upgrade remaining photon propagation lenght in case it is transmitted to different mus domain
      prop *= mus[phot->curel] / mus[phot->nextel];

      // Test transmission from vacuum -> scattering media
      if ((mus[phot->curel] <= 0.0) && (mus[phot->nextel] > 0.0))
      {
        // Draw new propagation distance -- otherwise photon might travel without scattering
        prop = -log(UnifOpen()) / mus[phot->nextel];
      }

      // Test for surival of the photon via roulette
      if (phot->weight < weight0)
      {
        if (UnifClosed() > chance)
          return;
        phot->weight /= chance;
      }

      // Fresnel transmission/reflection
      if (n[phot->curel] != n[phot->nextel])
      {
        if (FresnelPhoton(phot))
          continue;
      }

      // Update current face of the photon to that face which it will be on in the next element
      if (HN(phot->nextel, 0) == phot->curel)
        phot->curface = 0;
      else if (HN(phot->nextel, 1) == phot->curel)
        phot->curface = 1;
      else if (HN(phot->nextel, 2) == phot->curel)
        phot->curface = 2;
      else
      {
        loss++;
        return;
      }

      // Update current element of the photon
      phot->curel = phot->nextel;
    }

    // Scatter photon
    if (mus[phot->curel] > 0.0)
      ScatterPhoton(phot);
  }
}

// Run Monte Carlo
void MC2D::MonteCarlo(bool (*progress)(double))
{
#ifdef USE_OMP

  // OpenMP implementation

  // Spawn new MC2D classes with Nphoton' = Nphoton / ThreadCount, and initialize mt_rng seed
  long ii, jj, nthread = omp_get_max_threads(); // [AL]
  if (nthread <= 0)
  {
    printf("OpenMP initialization problem!");
    return;
  }
  MC2D *MCS = new MC2D[nthread];
  long *ticks = new long[nthread];
  bool abort_computation = false;

  // [AL] the progress bar gets an update after every TICK_VAL photons
#define TICK_VAL 1000

  for (ii = 0; ii < nthread; ii++)
  {
    MCS[ii] = *this;
    MCS[ii].Nphoton = Nphoton / nthread;
    MCS[ii].seed = MCS[ii].seed * totalthreads + ii;
    MCS[ii].InitRand();
    ticks[ii] = 0;
  }

  // [AL] if remainder of nphoton / nthread is non-zero, total photon count is not the same as Nphoton
  // therefore add the remaining photons to the last thread.
  long realnphot = 0;
  for (ii = 0; ii < nthread; ii++)
    realnphot += MCS[ii].Nphoton;
  MCS[nthread - 1].Nphoton += Nphoton - realnphot;
  // Compute Monte Carlo on each thread separetely
#pragma omp parallel
  {
    int_fast64_t iphoton, thread = omp_get_thread_num();
    Photon phot;
    for (iphoton = 1; iphoton <= MCS[thread].Nphoton; iphoton++)
    {
      if (iphoton % TICK_VAL == 0)
      {
        ticks[thread] = iphoton;
#pragma omp critical
        if (thread == 0)
        {
          int_fast64_t jj, csum = 0;
          {
            for (jj = 0; jj < nthread; jj++)
            {
              csum += ticks[jj];
            }
            if (!progress(100 * ((double)csum / (double)Nphoton)))
            {
              abort_computation = true;
            }
          }
        }
        if (abort_computation)
          break;
      }
      MCS[thread].CreatePhoton(&phot);
      MCS[thread].PropagatePhoton(&phot);
    }
  }
#ifdef VALOMC_MEX
  int_fast64_t csum = 0;
  for (jj = 0; jj < nthread; jj++)
  {
    csum += ticks[jj];
  }
  if (csum != Nphoton)
  {
    mexPrintf("WARNING: RUN WAS ABORTED OR PARALLEL ENVIRONMENT IS NOT WORKING CORRECTLY. IF YOU DID NOT ABORT THE RUN, PLEASE COMPILE AGAIN WITH OPENMP SUPPORT. \n");
  }
#endif
#pragma omp barrier

  // Sum up the results to first instance and delete MCS
  Nphoton = 0;
  loss = 0;
  for (jj = 0; jj < nthread; jj++)
  {
    Nphoton += MCS[jj].Nphoton;
    loss += MCS[jj].loss;
  }
  for (ii = 0; ii < H.Nx; ii++)
  {
    ER[ii] = EI[ii] = 0.0;
    for (jj = 0; jj < nthread; jj++)
    {
      ER[ii] += MCS[jj].ER[ii];
      EI[ii] += MCS[jj].EI[ii];
    }
  }
  for (ii = 0; ii < BH.Nx; ii++)
  {
    EBR[ii] = EBI[ii] = 0.0;
    for (jj = 0; jj < nthread; jj++)
    {
      EBR[ii] += MCS[jj].EBR[ii];
      EBI[ii] += MCS[jj].EBI[ii];
    }
  }

  for (ii = 0; ii < BH.Nx; ii++) // [AL]
  {
    DEBR[ii] = DEBI[ii] = 0.0;
    for (jj = 0; jj < nthread; jj++)
    {
      DEBR[ii] += MCS[jj].DEBR[ii];
      DEBI[ii] += MCS[jj].DEBI[ii];
    }
  }
  delete[] ticks;
  delete[] MCS;

#else

  // Single thread implementation
  long ii;

  //long itick = max(1, Nphoton / 100);
  long itick = 5000;
  double percentage = 0;
  long iphoton = 0;
  Photon phot;

  for (iphoton = 0; iphoton < Nphoton; iphoton++)
  {
    if ((iphoton % itick == 0) && (progress != NULL))
    {
      percentage = ((double)100.0 * iphoton / (double)Nphoton);
      if (!progress(percentage))
        break;
    }
    CreatePhoton(&phot);
    PropagatePhoton(&phot);
  }

#endif

#ifdef USE_MPI
  // Sum up computation from each computer & normalize
  AllReduceArray(ER, MPI_SUM);
  AllReduceArray(EI, MPI_SUM);
  AllReduceArray(EBR, MPI_SUM);
  AllReduceArray(EBI, MPI_SUM);
  AllReduceArray(DEBR, MPI_SUM);
  AllReduceArray(DEBI, MPI_SUM);

  // Sum up computed photons
  long tmplong;
  MPI_Allreduce(&Nphoton, &tmplong, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  Nphoton = tmplong;
  // Sum up lost photons
  MPI_Allreduce(&loss, &tmplong, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  loss = tmplong;
#endif

  // Normalize output variables
  if (omega <= 0.0)
  {
    for (ii = 0; ii < H.Nx; ii++)
    {
      if (mua[ii] > 0.0)
        ER[ii] /= mua[ii] * ElementArea(ii) * (double)Nphoton;
      else
        ER[ii] /= ElementArea(ii) * (double)Nphoton;
    }
    for (ii = 0; ii < BH.Nx; ii++)
      EBR[ii] /= (double)Nphoton * ElementLength(ii);
    for (ii = 0; ii < BH.Nx; ii++)
      DEBR[ii] /= (double)Nphoton * ElementLength(ii);
  }
  else
  {
    for (ii = 0; ii < H.Nx; ii++)
    {
      double a = ER[ii], b = EI[ii];
      ER[ii] = (b * k[ii] + a * mua[ii]) / (pow(k[ii], 2) + pow(mua[ii], 2)) / (double)Nphoton / ElementArea(ii);
      EI[ii] = -(a * k[ii] - b * mua[ii]) / (pow(k[ii], 2) + pow(mua[ii], 2)) / (double)Nphoton / ElementArea(ii);
    }

    for (ii = 0; ii < BH.Nx; ii++)
    {
      EBR[ii] /= (double)Nphoton * ElementLength(ii);
      EBI[ii] /= (double)Nphoton * ElementLength(ii);
      DEBR[ii] /= (double)Nphoton * ElementLength(ii);
      DEBI[ii] /= (double)Nphoton * ElementLength(ii);
    }
  }

  if (progress != NULL)
    progress(100);
}

// Check if ray R + t D intersects with a segment (U, V)
int RayRayIntersects(double R[2], double D[2], double U[2], double V[2], double *t)
{
  // [ dx, ux-vx ] [ t ] = [ ux-rx ]
  // [ dy, uy-vy ] [ s ]   [ uy-ry ]
  double det = D[0] * (U[1] - V[1]) - D[1] * (U[0] - V[0]);
  if ((-eps < det) && (det < eps))
    return (0);
  double s = (D[1] * (R[0] - U[0]) - D[0] * (R[1] - U[1])) / det;
  if ((s < 0.0) || (s > 1.0))
    return (0);
  *t = ((R[1] - U[1]) * (U[0] - V[0]) - (R[0] - U[0]) * (U[1] - V[1])) / det;
  return (1);
}

#endif
