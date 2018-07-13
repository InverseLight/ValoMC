#include "version.h"
#include <stdio.h>
#include <string.h>

void version_string(char *buf) {
  
  sprintf(buf,"--------------------------------------------\n");
  sprintf(buf+strlen(buf),"  Version:  %s\n", build_version);
  sprintf(buf+strlen(buf),"  Revision: %s\n", build_revision);
#ifdef USE_OMP
  sprintf(buf+strlen(buf),"  OpenMP enabled                     \n");
#else
  sprintf(buf+strlen(buf),"  OpenMP disabled, no parallelization\n");
#endif
#ifdef USE_OMP
  sprintf(buf+strlen(buf),"  Using %d threads\n", omp_get_max_threads());
#endif
  sprintf(buf+strlen(buf),"--------------------------------------------\n");
}

