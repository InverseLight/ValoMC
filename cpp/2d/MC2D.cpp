
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <time.h>

#include "../version.h"

// define USE_OMP to utilize threading with OpenMP 
//#define USE_OMP

#ifdef USE_OMP
#include <omp.h>
#endif

#include "Array.hpp"
#include "MC2D.hpp"


MC2D MC;


int LoadProblem(char *fin);
int SaveProblem(char *fout, int time);
bool Progress(double perc);



int main(int argc, char **argv) {

  printf("-----ValoMC-2D-----\n");
  printf("---Version %5s --\n", build_version);
  printf("---Revision%5s---\n", build_revision);
#ifdef USE_OMP
  printf("--OpenMP version --\n");
#endif
  printf("-------------------\n");
  printf("\n");

#ifdef USE_OMP
  printf("Using %d threads\n", omp_get_max_threads());
 //  double tstart = omp_get_wtime();
#endif

  MC.seed = (unsigned long) time(NULL);

  // Display help
  if( (argc < 3) ){
    printf("Use syntax:\n");
    printf(" MC2D inputfile outputfile\n");
    printf("\n");
    printf("Authors: Aki Pulkkinen, Aleksi Leino and Tanja Tarvainen (2018).\n");
    printf("The simulation code is originally written by Aki Pulkkinen.\n");
    printf("\n");
    return(0);
  }
  
  char *fin = argv[1];
  char *fout = argv[2];

  printf("Loading problem %s\n", fin);
  
  if(LoadProblem(fin)){
    printf("Error while loading!\n");
    return(1);
  }

  printf("Initializing MC2D\n");
   
  // [AL]
  try {
      MC.ErrorChecks();
      MC.Init();
  } catch (mcerror e) {
    printf("MC Initialization failed!: Reason %s\n", errorstring(e));
    return 1;
  }
 
  printf("Computing...\n");
  int tstart = (unsigned long) time(NULL);
  MC.MonteCarlo(Progress);
  int tend = (unsigned long) time(NULL);

  printf("Saving problem %s\n", fout);

  if(SaveProblem(fout, tend - tstart)) {
    printf("Error while saving!\n");
    return(1);
  }

  if(MC.loss) printf(" %ld photons lost during computation!\n", MC.loss);

//#ifdef USE_OMP
//  double tend = omp_get_wtime();

  printf("Computation took %i seconds\n", tend - tstart);
//#endif

  return(0);
}

int LoadProblem_TXT(char *fin){
  /*
    File Structure for text input:

      Ne Nb Nr Nphoton
      f phase0
      nx ny nz [AL]
      H
      BH
      r
      mua mus g n
      BCType
      [BCn]
      [BCLNormal] -- if provided, BCn has to be provided as well
      [BCLightDirectionType] // [AL]
      [BCIntensity] // [AL]
      [GaussianSigma] //[AL]
*/

  long ii;
  long Ne, Nb, Nr;
  int fsr;

  FILE *fp = fopen(fin, "r");
  if(fp == NULL) return(1);

  fsr=fscanf(fp, "%ld %ld %ld %ld\n", &Ne, &Nb, &Nr, &MC.Nphoton);
  fsr=fscanf(fp, "%lf %lf\n", &MC.f, &MC.phase0); // [AL]
  fsr=fscanf(fp, "%*i %*i\n"); // skip nx and ny [AL]
  // make negative phase0 positive by adding a multiple of 2*pi
  if(MC.phase0 < 0) {
    MC.phase0 += 2*M_PI*ceil(-MC.phase0 / (2*M_PI));
  }
  MC.H.resize(Ne, 3);
  MC.BH.resize(Nb, 2);
  MC.r.resize(Nr, 2);
  MC.BCType.resize(Nb);
  MC.BCLNormal.resize(Nb, 2);
  MC.BCn.resize(Nb);
  MC.mua.resize(Ne);
  MC.mus.resize(Ne);
  MC.g.resize(Ne);
  MC.n.resize(Ne);

  for(ii = 0; ii < Ne; ii++) fsr=fscanf(fp, "%ld %ld %ld\n", &MC.H(ii, 0), &MC.H(ii, 1), &MC.H(ii, 2));
  for(ii = 0; ii < Nb; ii++) fsr=fscanf(fp, "%ld %ld\n", &MC.BH(ii, 0), &MC.BH(ii, 1));
  for(ii = 0; ii < Nr; ii++) fsr=fscanf(fp, "%le %le\n", &MC.r(ii, 0), &MC.r(ii, 1));
  for(ii = 0; ii < Ne; ii++) fsr=fscanf(fp, "%le %le %le %le\n", &MC.mua(ii), &MC.mus(ii), &MC.g(ii), &MC.n(ii));
  for(ii = 0; ii < Nb; ii++) fsr=fscanf(fp, "%c\n", &MC.BCType(ii)); 

  // If the file still continues read BCn & BCLNormal
  if(!feof(fp)){
    for(ii = 0; ii < Nb; ii++) fsr=fscanf(fp, "%le\n", &MC.BCn(ii)); 
    if(!feof(fp)) for(ii = 0; ii < Nb; ii++) { 
       fsr=fscanf(fp, "%le %le\n", &MC.BCLNormal(ii, 0), &MC.BCLNormal(ii, 1));
    } else MC.BCLNormal.destroy();

    if(!feof(fp)) for(ii = 0; ii < Nb; ii++) fsr=fscanf(fp, "%c\n", &MC.BCLightDirectionType(ii));   
    else MC.BCLightDirectionType.destroy();

    if(!feof(fp)) for(ii = 0; ii < Nb; ii++) fsr=fscanf(fp, "%le\n", &MC.BCIntensity(ii));
    else MC.BCIntensity.destroy();

    if(!feof(fp)) for(ii = 0; ii < Nb; ii++) fsr=fscanf(fp, "%le\n", &MC.GaussianSigma(ii));  
    else MC.GaussianSigma.destroy();
  }
  else{
    MC.BCn.destroy();
    MC.BCLNormal.destroy();
    MC.BCLightDirectionType.destroy();
    MC.BCIntensity.destroy();
    MC.GaussianSigma.destroy();
  }
  
  fclose(fp);
  
  printf("Loaded:\n");
  printf("         H   (%ld x %ld)\n", MC.H.Nx, MC.H.Ny);
  printf("        BH   (%ld x %ld)\n", MC.BH.Nx, MC.BH.Ny);
  printf("         r   (%ld x %ld)\n", MC.r.Nx, MC.r.Ny);
  printf("    BCType   (%ld)\n", MC.BCType.Nx);
  printf("       BCn   (%ld)\n", MC.BCn.Nx);
  printf(" BCLNormal   (%ld x %ld)\n", MC.BCLNormal.Nx, MC.BCLNormal.Ny);
  printf("       mua   (%ld)\n", MC.mua.Nx);
  printf("       mus   (%ld)\n", MC.mus.Nx);
  printf("         g   (%ld)\n", MC.g.Nx);
  printf("         n   (%ld)\n", MC.n.Nx);
  printf("Constants:\n");
  printf("         f        (%e)\n", MC.f);
  printf("         phase0   (%e)\n", MC.phase0);  
  return(0);
}

int LoadProblem(char *fin){
  return( LoadProblem_TXT(fin) );
}

int SaveProblem_TXT(char *fout, int time){
  int ii;

  FILE *fp = fopen(fout, "w");
  if(fp == NULL) return(1);


  fprintf(fp, "%i 0\n", time);
  for(ii = 0; ii < MC.ER.Nx; ii++) fprintf(fp, "%e %e\n", MC.ER[ii], MC.EI[ii]);
  for(ii = 0; ii < MC.EBR.Nx; ii++) fprintf(fp, "%e %e\n", MC.EBR[ii], MC.EBI[ii]);  
  fclose(fp);

  return(0);
}

int SaveProblem(char *fin, int time){
  return( SaveProblem_TXT(fin, time) );
}

bool Progress(double perc){
  printf("  %f %%\r", perc);
  fflush(stdout);
  return true;
}

