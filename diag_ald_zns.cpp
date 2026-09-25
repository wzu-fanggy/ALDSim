/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */

#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "diag_ald_zns.h"
#include "app.h"
#include "app_ald_zns.h"
#include "comm_lattice.h"
#include "timer.h"
#include "error.h"
#include "memory.h"

using namespace SPPARKS_NS;


enum{VACANCY,S,SH,//2
SH2, ZnX2S, ZnX2SH, ZnX2SH2, // 6 DEZ adsorption 
ZnXS, ZnXSH, ZnS, ZnSH, Zn, ZnX, // 13 DEZ surface species
SH2Zn, SH2ZnX, SHZn, SHZnX, SZn, // 18 water pulse
Zn_i, ZnX_i, // 20 inert sites
QCM, SULFUR, ZINC, ADS_DEZ, HYDROGEN, DEZ, MEZ, LIGANDS, EVENTS,ONE,TWO,THREE
};       // same as DiagAld


/* ---------------------------------------------------------------------- */

DiagAldZns::DiagAldZns(SPPARKS *spk, int narg, char **arg) : Diag(spk,narg,arg)
{
  if (strcmp(app->style,"ald/zns") != 0)
    error->all(FLERR,"Diag style incompatible with app style");

  nlist = 0;

  int iarg = iarg_child;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"list") == 0) {
      nlist = narg - iarg - 1;
      list = new char*[nlist];
      int j = 0;
      for (int i = iarg+1; i < narg; i++) {
	int n = strlen(arg[i]) + 1;
	list[j] = new char[n];
	strcpy(list[j],arg[i]);
	j++;
      }
      iarg = narg;
    } else error->all(FLERR,"Diag_style aldZnS requires app_style aldZnS");                                      
  }

  if (nlist == 0) error->all(FLERR,"Diag_style aldZnS requires app_style aldZnS");
  which = new int[nlist];
  index = new int[nlist];
  ivector = new int[nlist];
}

/* ---------------------------------------------------------------------- */

DiagAldZns::~DiagAldZns()
{
  for (int i = 0; i < nlist; i++) delete [] list[i];
  delete [] list;
  delete [] which;
  delete [] index;
  delete [] ivector;
}

/* ---------------------------------------------------------------------- */

void DiagAldZns::init()
{
  appaldzns = (AppAldZns *) app;
  
  int none = appaldzns->none;
  int ntwo = appaldzns->ntwo;
  int nthree = appaldzns->nthree;
  for (int i = 0; i < nlist; i++) {
      if (strcmp(list[i],"S") == 0) which[i] = S;
      else if (strcmp(list[i],"SH") == 0) which[i] = SH;
      else if (strcmp(list[i],"SH2") == 0) which[i] = SH2;
      else if (strcmp(list[i],"VAC") == 0) which[i] = VACANCY;
      else if (strcmp(list[i],"HYDROGEN") == 0) which[i] = HYDROGEN;
      else if (strcmp(list[i],"ADS_DEZ") == 0) which[i] = ADS_DEZ;
      else if (strcmp(list[i],"QCM") == 0) which[i] = QCM;
      else if (strcmp(list[i],"SULFUR") == 0) which[i] = SULFUR;
      else if (strcmp(list[i],"ZINC") == 0) which[i] = ZINC;
      else if (strcmp(list[i],"DEZ") == 0) which[i] = DEZ;
      else if (strcmp(list[i],"MEZ") == 0) which[i] = MEZ;
      else if (strcmp(list[i],"LIGANDS") == 0) which[i] = LIGANDS;
      else if (strcmp(list[i],"events") == 0) which[i] = EVENTS;
      else if (strcmp(list[i],"ZnX2S") == 0) which[i] = ZnX2S;
      else if (strcmp(list[i],"ZnX2SH") == 0) which[i] = ZnX2SH;
      else if (strcmp(list[i],"ZnX2SH2") == 0) which[i] = ZnX2SH2;
      else if (strcmp(list[i],"ZnXS") == 0) which[i] = ZnXS;
      else if (strcmp(list[i],"ZnXSH") == 0) which[i] = ZnXSH;
      else if (strcmp(list[i],"ZnS") == 0) which[i] = ZnS;
      else if (strcmp(list[i],"ZnSH") == 0) which[i] = ZnSH;
      else if (strcmp(list[i],"Zn") == 0) which[i] = Zn;
      else if (strcmp(list[i],"ZnX") == 0) which[i] = ZnX;
      else if (strcmp(list[i],"SH2Zn") == 0) which[i] = SH2Zn;
      else if (strcmp(list[i],"SH2ZnX") == 0) which[i] = SH2ZnX;
      else if (strcmp(list[i],"SHZn") == 0) which[i] = SHZn;
      else if (strcmp(list[i],"SZn") == 0) which[i] = SZn;
      else if (strcmp(list[i],"SHZnX") == 0) which[i] = SHZnX;
      else if (strcmp(list[i],"Zn_i") == 0) which[i] = Zn_i;
      else if (strcmp(list[i],"ZnX_i") == 0) which[i] = ZnX_i;

    else if (list[i][0] == 's') {
      which[i] = ONE;
      int n = atoi(&list[i][1]);
      if (n < 1 || n > none) 
	error->all(FLERR,"Diag_style aldZnS requires app_style aldZnS");
      index[i] = n - 1;
    } else if (list[i][0] == 'd') {
      which[i] = TWO;
      int n = atoi(&list[i][1]);
      if (n < 1 || n > ntwo) 
	error->all(FLERR,"Diag_style aldZnS requires app_style aldZnS");
      index[i] = n - 1;
    } else if (list[i][0] == 'v') {
      which[i] = THREE;
      int n = atoi(&list[i][1]);
      if (n < 1 || n > nthree) 
	error->all(FLERR,"Diag_style aldZnS requires app_style aldZnS");
      index[i] = n - 1;
    } else error->all(FLERR,"Diag_style aldZnS requires app_style aldZnS");
  }

  siteflag = 1; 

  for (int i = 0; i < nlist; i++) ivector[i] = 0;
}

/* ---------------------------------------------------------------------- */

void DiagAldZns::compute()
{
  int sites[800],ivalue;
// here as well we have to consider some modification, generally it does not seem so difficult
  if (siteflag) {
    sites[S] = 0; sites[SH] = 0; sites[VACANCY] = 0; sites[SH2] = 0;
    sites[ZnX2S] = 0; sites[ZnX2SH] = 0;sites[ZnX2SH2] = 0; sites[ZnXS] = 0; sites[ZnXSH] = 0;   
    sites[ZnS] = 0; sites[ZnSH] = 0;sites[Zn] = 0; sites[ZnX] = 0; sites[SH2Zn] = 0; sites[SZn] = 0;   
    sites[SH2ZnX] = 0; sites[SHZn] = 0;sites[SHZnX] = 0; sites[Zn_i] = 0;sites[ZnX_i] = 0;


    int *element = appaldzns->element;
    int nlocal = appaldzns->nlocal;
    for (int i = 0; i < nlocal; i++) sites[element[i]]++;
  }

  for (int i = 0; i < nlist; i++) {
    if (which[i] == SH) ivalue = sites[SH];
    else if (which[i] == S) ivalue = sites[S];
    else if (which[i] == VACANCY) ivalue = sites[VACANCY];
    else if (which[i] == SH2) ivalue = sites[SH2];
    else if (which[i] == ZnX2S) ivalue = sites[ZnX2S];
    else if (which[i] == ZnX2SH) ivalue = sites[ZnX2SH];
    else if (which[i] == ZnX2SH2) ivalue = sites[ZnX2SH2];
    else if (which[i] == ZnXS) ivalue = sites[ZnXS];
    else if (which[i] == ZnXSH) ivalue = sites[ZnXSH];
    else if (which[i] == ZnS) ivalue = sites[ZnS];
    else if (which[i] == ZnSH) ivalue = sites[ZnSH];
    else if (which[i] == Zn) ivalue = sites[Zn];
    else if (which[i] == ZnX) ivalue = sites[ZnX];
    else if (which[i] == SH2Zn) ivalue = sites[SH2Zn];
    else if (which[i] == SH2ZnX) ivalue = sites[SH2ZnX];
    else if (which[i] == SHZn) ivalue = sites[SHZn];
    else if (which[i] == SHZnX) ivalue = sites[SHZnX];
    else if (which[i] == Zn_i) ivalue = sites[Zn_i];
    else if (which[i] == ZnX_i) ivalue = sites[ZnX_i];
    else if (which[i] == ADS_DEZ) ivalue = sites[ZnX2S] + sites[ZnX2SH] + sites[ZnX2SH2];
    else if (which[i] == HYDROGEN) ivalue = sites[SH] + 2*sites[SH2] + 2*sites[ZnX2SH2] + sites[ZnX2SH] + sites[ZnXSH] + sites[ZnSH] + 2*sites[SH2Zn] + sites[SHZn] + 2*sites[SH2ZnX] + sites[SHZnX];
    else if (which[i] == QCM) ivalue = 34.08*sites[SH2] + 33.07*sites[SH] + 32.06*sites[S] + 157.57*sites[ZnX2SH2] + 156.56*sites[ZnX2SH] + 155.55*sites[ZnX2S] + 127.51*sites[ZnXSH] + 126.50*sites[ZnXS] + 94.44*sites[ZnX] + 97.45*sites[ZnS] + 98.46*sites[ZnSH] + 65.39*sites[Zn] + 99.47*sites[SH2Zn] + 98.46*sites[SHZn] + 97.45*sites[SZn] + 128.52*sites[SH2ZnX] + 127.51*sites[SHZnX] + 65.39*sites[Zn_i] + 94.44*sites[ZnX_i];
    else if (which[i] == SULFUR) ivalue = sites[S] + sites[SH] + sites[SH2] + sites[ZnX2SH2] + sites[ZnX2SH] + sites[ZnX2S] + sites[ZnXS] + sites[ZnXSH] + sites[ZnS] + sites[ZnSH] + sites[SH2Zn] + sites[SHZn] + sites[SZn] + sites[SH2ZnX] + sites[SHZnX];
    else if (which[i] == ZINC) ivalue = sites[Zn] + sites[ZnX] + sites[ZnX2SH2] + sites[ZnX2SH] + sites[ZnX2S] + sites[ZnXS] + sites[ZnXSH] + sites[ZnS] + sites[ZnSH] + sites[SH2Zn] + sites[SHZn] + sites[SZn] + sites[SH2ZnX] + sites[SHZnX] + sites[Zn_i] + sites[ZnX_i];
    else if (which[i] == DEZ) ivalue = sites[ZnX2SH2] + sites[ZnX2SH] + sites[ZnX2S];
    else if (which[i] == MEZ) ivalue = sites[ZnXS] + sites[ZnXSH] + sites[ZnX] + sites[SH2ZnX] + sites[SHZnX] + sites[ZnX_i];
    else if (which[i] == LIGANDS) ivalue = 2 * (sites[ZnX2SH2] + sites[ZnX2SH] + sites[ZnX2S]) + sites[ZnXS] + sites[ZnXSH] + sites[ZnX] + sites[SH2ZnX] + sites[SHZnX] + sites[ZnX_i];
    else if (which[i] == EVENTS) ivalue = appaldzns->nevents;
    else if (which[i] == ONE) ivalue = appaldzns->scount[index[i]];
    else if (which[i] == TWO) ivalue = appaldzns->dcount[index[i]];
    else if (which[i] == THREE) ivalue = appaldzns->vcount[index[i]];
    
    MPI_Allreduce(&ivalue,&ivector[i],1,MPI_INT,MPI_SUM,world);
  }
}

/* ---------------------------------------------------------------------- */

void DiagAldZns::stats(char *str)
{
  for (int i = 0; i < nlist; i++) {
    sprintf(str,"%6d ",ivector[i]);
    str += strlen(str);
  }
}

/* ---------------------------------------------------------------------- */

void DiagAldZns::stats_header(char *str)
{
  for (int i = 0; i < nlist; i++) {
    sprintf(str,"%6s ",list[i]);
    str += strlen(str);
  }
}
