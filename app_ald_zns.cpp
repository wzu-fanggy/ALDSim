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

#include "math.h"
#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "app_ald_zns.h"
#include "solve.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS_NS;

enum{VACANCY,S,SH,// 3
SH2, ZnX2S, ZnX2SH, ZnX2SH2, // 7 DEZ adsorption 
ZnXS, ZnXSH, ZnS, ZnSH, Zn, ZnX, // 13 DEZ surface species
SH2Zn, SH2ZnX, SHZn, SHZnX, SZn}; // 18 water pulse



#define DELTAEVENT 100000
#define AUDIT_INTERVAL 1000000
#define MAX_HEALS 10000

// last executed event, kept for crash diagnostics (serial run)
static int g_last_rstyle = -1;
static int g_last_which = -1;
static int g_last_j = -2;
static int g_last_k = -2;
static double g_last_time = 0.0;

/* ---------------------------------------------------------------------- */

AppAldZns::AppAldZns(SPPARKS *spk, int narg, char **arg) : 
  AppLattice(spk,narg,arg)
{
  ninteger = 2;
  ndouble = 0;
  delpropensity = 1;
  delevent = 1;
  allow_kmc = 1;
  allow_rejection = 0;
  allow_masking = 0;
  

  create_arrays();

  if (narg != 1) error->all(FLERR,"Illegal app_style command");

  cycle = 0;
  pressureOn = 1;
  hello = 1;
  firsttime = 1;
  esites = NULL;
  echeck = NULL;
  msites = NULL;
  nmsites = 0;
  events = NULL;
  maxevent = 0;
  firstevent = NULL;
  audit_counter = 0;
  heal_count = 0;

  // reaction lists

  none = ntwo = nthree = 0;
  srate = drate = vrate = NULL;
  spropensity = dpropensity = vpropensity = NULL;
  sinput = soutput = NULL;
  dinput = doutput = NULL;
  vinput = voutput = NULL;
  comneigh = NULL;
  comneigh_size = 0;
  scount = dcount = vcount = NULL;
  sA = dA = vA = NULL;
  scoord = dcoord = vcoord = NULL;
  sexpon = dexpon = vexpon = NULL;
  spresson = dpresson = vpresson = NULL;
}

/* ---------------------------------------------------------------------- */

AppAldZns::~AppAldZns()
{
  delete [] esites;
  delete [] echeck;
  memory->sfree(msites);
  memory->sfree(events);
  memory->sfree(firstevent);
  memory->sfree(srate);
  memory->sfree(drate);
  memory->sfree(vrate);
  memory->sfree(spropensity);
  memory->sfree(dpropensity);
  memory->sfree(vpropensity);
  memory->sfree(sinput);
  memory->sfree(soutput);
  memory->sfree(dinput);
  memory->sfree(doutput);
  memory->sfree(vinput);
  memory->sfree(voutput);
  memory->sfree(comneigh);
  memory->sfree(scount);
  memory->sfree(dcount);
  memory->sfree(vcount);
  memory->sfree(sA);
  memory->sfree(dA);
  memory->sfree(vA);
  memory->sfree(scoord);
  memory->sfree(dcoord);
  memory->sfree(vcoord);
  memory->sfree(sexpon);
  memory->sfree(dexpon);
  memory->sfree(vexpon);
  memory->sfree(spresson);
  memory->sfree(dpresson);
  memory->sfree(vpresson);
}

/* ---------------------------------------------------------------------- */

void AppAldZns::input_app(char *command, int narg, char **arg)
{
  if (strcmp(command,"event") == 0) {
    if (narg < 1) error->all(FLERR,"Illegal event command E1");
    int rstyle = atoi(arg[0]);
    grow_reactions(rstyle);

    if (rstyle == 1) {
      if (narg != 9) error->all(FLERR,"Illegal event arg command");
//type I
      if (strcmp(arg[1],"S") == 0) sinput[none] = S;
      else if (strcmp(arg[1],"SH") == 0) sinput[none] = SH; 
      else if (strcmp(arg[1],"SH2") == 0) sinput[none] = SH2;
      else if (strcmp(arg[1],"Zn") == 0) sinput[none] = Zn;
      else if (strcmp(arg[1],"ZnX") == 0) sinput[none] = ZnX;
      else if (strcmp(arg[1],"ZnX2SH") == 0) sinput[none] = ZnX2SH;
      else if (strcmp(arg[1],"ZnX2SH2") == 0) sinput[none] = ZnX2SH2;
      else if (strcmp(arg[1],"ZnX2S") == 0) sinput[none] = ZnX2S;
      else if (strcmp(arg[1],"ZnXSH") == 0) sinput[none] = ZnXSH;
      else if (strcmp(arg[1],"ZnXS") == 0) sinput[none] = ZnXS;
      else if (strcmp(arg[1],"ZnSH") == 0) sinput[none] = ZnSH;
      else if (strcmp(arg[1],"ZnS") == 0) sinput[none] = ZnS;
      else if (strcmp(arg[1],"SHZn") == 0) sinput[none] = SHZn;
      else if (strcmp(arg[1],"SH2ZnX") == 0) sinput[none] = SH2ZnX;
      else if (strcmp(arg[1],"SHZnX") == 0) sinput[none] = SHZnX;
      else if (strcmp(arg[1],"SH2Zn") == 0) sinput[none] = SH2Zn;
      else if (strcmp(arg[1],"SZn") == 0) sinput[none] = SZn;


	
      else error->all(FLERR,"Illegal event arg1 command");

      if (strcmp(arg[2],"S") == 0) soutput[none] = S;
      else if (strcmp(arg[2],"SH") == 0) soutput[none] = SH;
      else if (strcmp(arg[2],"VAC") == 0) soutput[none] = VACANCY;
      else if (strcmp(arg[2],"SH2") == 0) soutput[none] = SH2;
      else if (strcmp(arg[2],"Zn") == 0) soutput[none] = Zn;
      else if (strcmp(arg[2],"ZnX") == 0) soutput[none] = ZnX;
      else if (strcmp(arg[2],"ZnX2SH") == 0) soutput[none] = ZnX2SH;
      else if (strcmp(arg[2],"ZnX2SH2") == 0) soutput[none] = ZnX2SH2;
      else if (strcmp(arg[2],"ZnX2S") == 0) soutput[none] = ZnX2S;
      else if (strcmp(arg[2],"ZnXSH") == 0) soutput[none] = ZnXSH;
      else if (strcmp(arg[2],"ZnXS") == 0) soutput[none] = ZnXS;
      else if (strcmp(arg[2],"ZnSH") == 0) soutput[none] = ZnSH;
      else if (strcmp(arg[2],"ZnS") == 0) soutput[none] = ZnS;
      else if (strcmp(arg[2],"SHZn") == 0) soutput[none] = SHZn;
      else if (strcmp(arg[2],"SH2ZnX") == 0) soutput[none] = SH2ZnX;
      else if (strcmp(arg[2],"SHZnX") == 0) soutput[none] = SHZnX;
      else if (strcmp(arg[2],"SH2Zn") == 0) soutput[none] = SH2Zn;
      else if (strcmp(arg[2],"SZn") == 0) soutput[none] = SZn;

      
      else error->all(FLERR,"Illegal event command E2");
      
      sA[none] = atof(arg[3]);
      if (sA[none] == 0.0) error->warning(FLERR,"Illegal coef during reading command");
      sexpon[none] = atoi(arg[4]);
      srate[none] = atof(arg[5]);
      scoord[none] = atoi(arg[6]);
      spresson[none] = atoi(arg[7]);

      none++;
      
//type II 
    } else if (rstyle == 2) {
      if (narg != 11) error->all(FLERR,"Illegal event command E3");

      if (strcmp(arg[1],"S") == 0) dinput[ntwo][0] = S;
      else if (strcmp(arg[1],"SH") == 0) dinput[ntwo][0] = SH;
      else if (strcmp(arg[1],"SH2") == 0) dinput[ntwo][0] = SH2;
      else if (strcmp(arg[1],"Zn") == 0) dinput[ntwo][0] = Zn;
      else if (strcmp(arg[1],"ZnX") == 0) dinput[ntwo][0] = ZnX;
      else if (strcmp(arg[1],"ZnX2SH") == 0) dinput[ntwo][0] = ZnX2SH;
      else if (strcmp(arg[1],"ZnX2SH2") == 0) dinput[ntwo][0] = ZnX2SH2;
      else if (strcmp(arg[1],"ZnX2S") == 0) dinput[ntwo][0] = ZnX2S;
      else if (strcmp(arg[1],"ZnXSH") == 0) dinput[ntwo][0] = ZnXSH;
      else if (strcmp(arg[1],"ZnXS") == 0) dinput[ntwo][0] = ZnXS;
      else if (strcmp(arg[1],"ZnSH") == 0) dinput[ntwo][0] = ZnSH;
      else if (strcmp(arg[1],"ZnS") == 0) dinput[ntwo][0] = ZnS;
      else if (strcmp(arg[1],"SHZn") == 0) dinput[ntwo][0] = SHZn;
      else if (strcmp(arg[1],"SH2ZnX") == 0) dinput[ntwo][0] = SH2ZnX;
      else if (strcmp(arg[1],"SHZnX") == 0) dinput[ntwo][0] = SHZnX;
      else if (strcmp(arg[1],"SH2Zn") == 0) dinput[ntwo][0] = SH2Zn;
      else if (strcmp(arg[1],"SZn") == 0) dinput[ntwo][0] = SZn;

	
      else error->all(FLERR,"Illegal event command E4");
      
      if (strcmp(arg[2],"S") == 0) doutput[ntwo][0] = S;
      else if (strcmp(arg[2],"SH") == 0) doutput[ntwo][0] = SH;
      else if (strcmp(arg[2],"SH2") == 0) doutput[ntwo][0] = SH2;
      else if (strcmp(arg[2],"Zn") == 0) doutput[ntwo][0] = Zn;
      else if (strcmp(arg[2],"ZnX") == 0) doutput[ntwo][0] = ZnX;
      else if (strcmp(arg[2],"ZnX2SH") == 0) doutput[ntwo][0] = ZnX2SH;
      else if (strcmp(arg[2],"ZnX2SH2") == 0) doutput[ntwo][0] = ZnX2SH2;
      else if (strcmp(arg[2],"ZnX2S") == 0) doutput[ntwo][0] = ZnX2S;
      else if (strcmp(arg[2],"ZnXSH") == 0) doutput[ntwo][0] = ZnXSH;
      else if (strcmp(arg[2],"ZnXS") == 0) doutput[ntwo][0] = ZnXS;
      else if (strcmp(arg[2],"ZnSH") == 0) doutput[ntwo][0] = ZnSH;
      else if (strcmp(arg[2],"ZnS") == 0) doutput[ntwo][0] = ZnS;
      else if (strcmp(arg[2],"SHZn") == 0) doutput[ntwo][0] = SHZn;
      else if (strcmp(arg[2],"SH2ZnX") == 0) doutput[ntwo][0] = SH2ZnX;
      else if (strcmp(arg[2],"SHZnX") == 0) doutput[ntwo][0] = SHZnX;
      else if (strcmp(arg[2],"SH2Zn") == 0) doutput[ntwo][0] = SH2Zn;
      else if (strcmp(arg[2],"SZn") == 0) doutput[ntwo][0] = SZn;

	
      else error->all(FLERR,"Illegal event command2");

      if (strcmp(arg[3],"SH") == 0) dinput[ntwo][1] = SH;
      else if (strcmp(arg[3],"S") == 0) dinput[ntwo][1] = S;
      else if (strcmp(arg[3],"SH2") == 0) dinput[ntwo][1] = SH2;
      else if (strcmp(arg[3],"Zn") == 0) dinput[ntwo][1] = Zn;
      else if (strcmp(arg[3],"ZnX") == 0) dinput[ntwo][1] = ZnX;
      else if (strcmp(arg[3],"ZnX2SH") == 0) dinput[ntwo][1] = ZnX2SH;
      else if (strcmp(arg[3],"ZnX2SH2") == 0) dinput[ntwo][1] = ZnX2SH2;
      else if (strcmp(arg[3],"ZnX2S") == 0) dinput[ntwo][1] = ZnX2S;
      else if (strcmp(arg[3],"ZnXSH") == 0) dinput[ntwo][1] = ZnXSH;
      else if (strcmp(arg[3],"ZnXS") == 0) dinput[ntwo][1] = ZnXS;
      else if (strcmp(arg[3],"ZnSH") == 0) dinput[ntwo][1] = ZnSH;
      else if (strcmp(arg[3],"ZnS") == 0) dinput[ntwo][1] = ZnS;
      else if (strcmp(arg[3],"SHZn") == 0) dinput[ntwo][1] = SHZn;
      else if (strcmp(arg[3],"SH2ZnX") == 0) dinput[ntwo][1] = SH2ZnX;
      else if (strcmp(arg[3],"SHZnX") == 0) dinput[ntwo][1] = SHZnX;
      else if (strcmp(arg[3],"SH2Zn") == 0) dinput[ntwo][1] = SH2Zn;
      else if (strcmp(arg[3],"SZn") == 0) dinput[ntwo][1] = SZn;

	
      else error->all(FLERR,"Illegal event command2");

      if (strcmp(arg[4],"S") == 0) doutput[ntwo][1] = S;
      else if (strcmp(arg[4],"SH") == 0) doutput[ntwo][1] = SH;
      else if (strcmp(arg[4],"SH2") == 0) doutput[ntwo][1] = SH2;
      else if (strcmp(arg[4],"Zn") == 0) doutput[ntwo][1] = Zn;
      else if (strcmp(arg[4],"ZnX") == 0) doutput[ntwo][1] = ZnX;
      else if (strcmp(arg[4],"ZnX2SH") == 0) doutput[ntwo][1] = ZnX2SH;
      else if (strcmp(arg[4],"ZnX2SH2") == 0) doutput[ntwo][1] = ZnX2SH2;
      else if (strcmp(arg[4],"ZnX2S") == 0) doutput[ntwo][1] = ZnX2S;
      else if (strcmp(arg[4],"ZnXSH") == 0) doutput[ntwo][1] = ZnXSH;
      else if (strcmp(arg[4],"ZnXS") == 0) doutput[ntwo][1] = ZnXS;
      else if (strcmp(arg[4],"ZnSH") == 0) doutput[ntwo][1] = ZnSH;
      else if (strcmp(arg[4],"ZnS") == 0) doutput[ntwo][1] = ZnS;
      else if (strcmp(arg[4],"SHZn") == 0) doutput[ntwo][1] = SHZn;
      else if (strcmp(arg[4],"SH2ZnX") == 0) doutput[ntwo][1] = SH2ZnX;
      else if (strcmp(arg[4],"SHZnX") == 0) doutput[ntwo][1] = SHZnX;
      else if (strcmp(arg[4],"SH2Zn") == 0) doutput[ntwo][1] = SH2Zn;
      else if (strcmp(arg[4],"SZn") == 0) doutput[ntwo][1] = SZn;
	
      else error->all(FLERR,"Illegal event command2");

      dA[ntwo] = atof(arg[5]);
      dexpon[ntwo] = atoi(arg[6]);
      if (dexpon[ntwo] != 0.0) error->warning(FLERR,"Illegal expon command2");
      drate[ntwo] = atof(arg[7]);
      dcoord[ntwo] = atoi(arg[8]);
      dpresson[ntwo] = atoi(arg[9]);
      ntwo++;
// type III
    }else if (rstyle == 3) {
      if (narg != 11) error->all(FLERR,"Illegal event command31");
 
      if (strcmp(arg[1],"S") == 0) vinput[nthree][0] = S;
      else if (strcmp(arg[1],"SH") == 0) vinput[nthree][0] = SH;
      else if (strcmp(arg[1],"SH2") == 0) vinput[nthree][0] = SH2;
      else if (strcmp(arg[1],"VAC") == 0) vinput[nthree][0] = VACANCY;
      else if (strcmp(arg[1],"Zn") == 0) vinput[nthree][0] = Zn;
      else if (strcmp(arg[1],"ZnX") == 0) vinput[nthree][0] = ZnX;
      else if (strcmp(arg[1],"ZnX2SH") == 0) vinput[nthree][0] = ZnX2SH;
      else if (strcmp(arg[1],"ZnX2SH2") == 0) vinput[nthree][0] = ZnX2SH2;
      else if (strcmp(arg[1],"ZnX2S") == 0) vinput[nthree][0] = ZnX2S;
      else if (strcmp(arg[1],"ZnXSH") == 0) vinput[nthree][0] = ZnXSH;
      else if (strcmp(arg[1],"ZnXS") == 0) vinput[nthree][0] = ZnXS;
      else if (strcmp(arg[1],"ZnSH") == 0) vinput[nthree][0] = ZnSH;
      else if (strcmp(arg[1],"ZnS") == 0) vinput[nthree][0] = ZnS;
      else if (strcmp(arg[1],"SHZn") == 0) vinput[nthree][0] = SHZn;
      else if (strcmp(arg[1],"SH2ZnX") == 0) vinput[nthree][0] = SH2ZnX;
      else if (strcmp(arg[1],"SHZnX") == 0) vinput[nthree][0] = SHZnX;
      else if (strcmp(arg[1],"SH2Zn") == 0) vinput[nthree][0] = SH2Zn;   
      else if (strcmp(arg[1],"SZn") == 0) vinput[nthree][0] = SZn;

      else error->all(FLERR,"Illegal event command32");

      if (strcmp(arg[2],"SH") == 0) voutput[nthree][0] = SH;
      else if (strcmp(arg[2],"S") == 0) voutput[nthree][0] = S;
      else if (strcmp(arg[2],"SH2") == 0) voutput[nthree][0] = SH2;
      else if (strcmp(arg[2],"VAC") == 0) voutput[nthree][0] = VACANCY;
      else if (strcmp(arg[2],"Zn") == 0) voutput[nthree][0] = Zn;
      else if (strcmp(arg[2],"ZnX") == 0) voutput[nthree][0] = ZnX;
      else if (strcmp(arg[2],"ZnX2SH") == 0) voutput[nthree][0] = ZnX2SH;
      else if (strcmp(arg[2],"ZnX2SH2") == 0) voutput[nthree][0] = ZnX2SH2;
      else if (strcmp(arg[2],"ZnX2S") == 0) voutput[nthree][0] = ZnX2S;
      else if (strcmp(arg[2],"ZnXSH") == 0) voutput[nthree][0] = ZnXSH;
      else if (strcmp(arg[2],"ZnXS") == 0) voutput[nthree][0] = ZnXS;
      else if (strcmp(arg[2],"ZnSH") == 0) voutput[nthree][0] = ZnSH;
      else if (strcmp(arg[2],"ZnS") == 0) voutput[nthree][0] = ZnS;
      else if (strcmp(arg[2],"SHZn") == 0) voutput[nthree][0] = SHZn;
      else if (strcmp(arg[2],"SH2ZnX") == 0) voutput[nthree][0] = SH2ZnX;
      else if (strcmp(arg[2],"SHZnX") == 0) voutput[nthree][0] = SHZnX;
      else if (strcmp(arg[2],"SH2Zn") == 0) voutput[nthree][0] = SH2Zn; 
      else if (strcmp(arg[2],"SZn") == 0) voutput[nthree][0] = SZn;

      else error->all(FLERR,"Illegal event command33");
      
      if (strcmp(arg[3],"VAC") == 0) vinput[nthree][1] = VACANCY;
      else if (strcmp(arg[3],"S") == 0) vinput[nthree][1] = S;
      else if (strcmp(arg[3],"SH2") == 0) vinput[nthree][1] = SH2;
      else if (strcmp(arg[3],"SH") == 0) vinput[nthree][1] = SH;
      else if (strcmp(arg[3],"Zn") == 0) vinput[nthree][1] = Zn;
      else if (strcmp(arg[3],"ZnX") == 0) vinput[nthree][1] = ZnX;
      else if (strcmp(arg[3],"ZnX2SH") == 0) vinput[nthree][1] = ZnX2SH;
      else if (strcmp(arg[3],"ZnX2SH2") == 0) vinput[nthree][1] = ZnX2SH2;
      else if (strcmp(arg[3],"ZnX2S") == 0) vinput[nthree][1] = ZnX2S;
      else if (strcmp(arg[3],"ZnXSH") == 0) vinput[nthree][1] = ZnXSH;
      else if (strcmp(arg[3],"ZnXS") == 0) vinput[nthree][1] = ZnXS;
      else if (strcmp(arg[3],"ZnSH") == 0) vinput[nthree][1] = ZnSH;
      else if (strcmp(arg[3],"ZnS") == 0) vinput[nthree][1] = ZnS;
      else if (strcmp(arg[3],"SHZn") == 0) vinput[nthree][1] = SHZn;
      else if (strcmp(arg[3],"SZn") == 0) vinput[nthree][1] = SZn;
      else if (strcmp(arg[3],"SH2ZnX") == 0) vinput[nthree][1] = SH2ZnX;
      else if (strcmp(arg[3],"SHZnX") == 0) vinput[nthree][1] = SHZnX;
      else if (strcmp(arg[3],"SH2Zn") == 0) vinput[nthree][1] = SH2Zn;   

    
      else error->all(FLERR,"Illegal event command34");

      if (strcmp(arg[4],"S") == 0) voutput[nthree][1] = S;
      else if (strcmp(arg[4],"SH") == 0) voutput[nthree][1] = SH;
      else if (strcmp(arg[4],"SH2") == 0) voutput[nthree][1] = SH2;
      else if (strcmp(arg[4],"VAC") == 0) voutput[nthree][1] = VACANCY;
      else if (strcmp(arg[4],"Zn") == 0) voutput[nthree][1] = Zn;
      else if (strcmp(arg[4],"ZnX") == 0) voutput[nthree][1] = ZnX;
      else if (strcmp(arg[4],"ZnX2SH") == 0) voutput[nthree][1] = ZnX2SH;
      else if (strcmp(arg[4],"ZnX2SH2") == 0) voutput[nthree][1] = ZnX2SH2;
      else if (strcmp(arg[4],"ZnX2S") == 0) voutput[nthree][1] = ZnX2S;
      else if (strcmp(arg[4],"ZnXSH") == 0) voutput[nthree][1] = ZnXSH;
      else if (strcmp(arg[4],"ZnXS") == 0) voutput[nthree][1] = ZnXS;
      else if (strcmp(arg[4],"ZnSH") == 0) voutput[nthree][1] = ZnSH;
      else if (strcmp(arg[4],"ZnS") == 0) voutput[nthree][1] = ZnS;
      else if (strcmp(arg[4],"SHZn") == 0) voutput[nthree][1] = SHZn;
      else if (strcmp(arg[4],"SH2ZnX") == 0) voutput[nthree][1] = SH2ZnX;
      else if (strcmp(arg[4],"SHZnX") == 0) voutput[nthree][1] = SHZnX;
      else if (strcmp(arg[4],"SH2Zn") == 0) voutput[nthree][1] = SH2Zn; 
      else if (strcmp(arg[4],"SZn") == 0) voutput[nthree][1] = SZn;

      else error->all(FLERR,"Illegal event command35");

      vA[nthree] = atof(arg[5]);
      vexpon[nthree] = atoi(arg[6]);
      if (vexpon[nthree] != 0.0) error->warning(FLERR,"Illegal vexpon command36");
      vrate[nthree] = atof(arg[7]);
      vcoord[nthree] = atoi(arg[8]);
      vpresson[nthree] = atoi(arg[9]);
      nthree++;

    } else error->all(FLERR,"Illegal event command37");
  } 
  else if (strcmp(command,"pulse_time") == 0) {
    if (narg != 2) error->all(FLERR,"Illegal pulse time");
      T1 = atof(arg[0]);
      T3 = atof(arg[1]);
  }
  else if (strcmp(command,"purge_time") == 0) {
    if (narg != 2) error->all(FLERR,"Illegal purge time");
      T2 = atof(arg[0]);
      T4 = atof(arg[1]);
  }else error->all(FLERR,"Unrecognized command38");
}

/* ----------------------------------------------------------------------
   set site value ptrs each time iarray/darray are reallocated
------------------------------------------------------------------------- */

void AppAldZns::grow_app()
{
  element = iarray[0];
  coord = iarray[1];
}

/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
------------------------------------------------------------------------- */

void AppAldZns::init_app()
{
  if (firsttime) {
    firsttime = 0;

    echeck = new int[nlocal];
    firstevent = (int *) memory->smalloc(nlocal*sizeof(int),"app:firstevent");
    //comneigh was defined to avoid double counting of common neighbor in site_propensity.
    // The second-neighbor scan can register up to ntwo*maxneigh*maxneigh distinct
    // (k,dpropensity) pairs per site, which exceeded nlocal in the ZnS 500K input
    // (201 type-II events, 4 max neighbors, 1280 sites) and silently corrupted memory.
    comneigh_size = 1 + ntwo * maxneigh * maxneigh;
    if (comneigh_size < nlocal) comneigh_size = nlocal;
    comneigh = memory->grow(comneigh,comneigh_size,2,"app/ald:comneigh");
    // esites must be large enough for 3 sites and their 1st neighbors
    
    esites = (int *) memory->smalloc(nlocal*sizeof(int),"app:esites");
    msites = (int *) memory->smalloc(nlocal*sizeof(int),"app:msites");
    //esites = new int[12*maxneigh]; 
  }
  // site validity

  int flag = 0;
  for (int i = 0; i < nlocal; i++) {
    if (coord[i] < -1 || coord[i] > 8) flag = 1;
    if (element[i] < VACANCY) flag = 1;
  }
  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall) error->all(FLERR,"One or more sites have invalid values");
}
/* ---------------------------------------------------------------------- */

void AppAldZns::setup_app()
{
  for (int i = 0; i < nlocal; i++) echeck[i] = 0;
  

  nevents = 0;
  for (int i = 0; i < nlocal; i++) firstevent[i] = -1;
  nmsites = 0;
  for (int i = 0; i < maxevent; i++) events[i].next = i+1;
  freeevent = 0;


  if (temperature == 0.0)
    error->all(FLERR,"Temperature cannot be 0.0 for app_ald");
  for (int m = 0; m < none; m++) {
    spropensity[m] = sA[m]*pow(temperature,sexpon[m])*exp(-srate[m]/temperature);
    scount[m] = 0;
  if (spropensity[m] == 0.0) error->warning(FLERR," spropensity cannot be 0.0 for app_ald");
  }
  for (int m = 0; m < ntwo; m++) {
    dpropensity[m] = dA[m]*pow(temperature,dexpon[m])*exp(-drate[m]/temperature);
    dcount[m] = 0;
  if (dpropensity[m] == 0.0) error->warning(FLERR,"dpropensity cannot be 0.0 for app_ald");
  }
  for (int m = 0; m < nthree; m++) {
    vpropensity[m] = vA[m]*pow(temperature,vexpon[m])*exp(-vrate[m]/temperature);
    vcount[m] = 0;
  if (vpropensity[m] == 0.0) error->warning(FLERR,"vpropensity cannot be 0.0 for app_ald");
  }
}

/* ----------------------------------------------------------------------
   compute energy of site
------------------------------------------------------------------------- */

double AppAldZns::site_energy(int i)
{
  return 0.0;
}

/* ----------------------------------------------------------------------
   KMC method
   compute total propensity of owned site summed over possible events
------------------------------------------------------------------------- */

double AppAldZns::site_propensity(int i)
{
  int j,k,m;


  clear_events(i);

  double proball = 0.0;


  //type I, check species of sites and consider possible events

  // count_coordS was added here to prevent adsorption of metal precursor
  // on the low coordinate SULFUR at sublayer

  for (m = 0; m < none; m++) {
    if (element[i] != sinput[m]) continue;
    if ((coord[i] == scoord[m] || scoord[m] == 0) && (spresson[m] == pressureOn || spresson[m] == 0)) {
	    add_event(i,1,m,spropensity[m],-1,-1);
	    proball += spropensity[m];
    }
  }

  // type II, check species of sites and second neighbor 
  // comneigh variable used to avoid double counting,
  // we consider more than one event between sites, therefore we need 2d comneigh array.

  int nextneib = 1;
  for (int jj = 0; jj < numneigh[i]; jj++) {
    j = neighbor[i][jj];
       for (int kk = 0; kk < numneigh[j]; kk++) {
            k = neighbor[j][kk];
            if (i == k) continue;
            for (m = 0; m < ntwo; m++) {
	    if ((element[i] == dinput [m][0] && element[k] == dinput [m][1]) && (dpresson[m] == pressureOn || dpresson[m] == 0) && (coord[i] == dcoord[m] || dcoord[m] == 0)) {
  	      comevent = 1;
  	      for (int ii = 0; ii < nextneib; ii++) {
  	      if ( comneigh[ii][0] == k && comneigh[ii][1] == dpropensity[m]) comevent = 0;
  	      }
  	      if (comevent){	
                add_event(i,2,m,dpropensity[m],-1,k);
                proball += dpropensity[m];
                if (nextneib >= comneigh_size)
                  error->all(FLERR,"comneigh overflow in app_ald_zns: nextneib exceeds comneigh_size");
  	        comneigh[nextneib][0] = k;
                comneigh[nextneib][1] = dpropensity[m];
  	        nextneib++;
  	    } 
          }
        }
      }
    }
    for (m = 0; m < nextneib; m++) comneigh[m][0] = comneigh[m][1] = 0;

  //type III, check species of sites and first neighbour 


  for (int jj = 0; jj < numneigh[i]; jj++) {
    j = neighbor[i][jj];
      for (m =0; m < nthree; m++) {
//	if (element[i] == vinput[m][0] && element[j] == vinput[m][1] && (coord[i] == vcoord[m] || vcoord[m] == 0) && (dpresson[m] == pressureOn || dpresson[m] == 0)) { // Zero cn does not allow reaction
	if (element[i] == vinput[m][0] && element[j] == vinput[m][1] && (coord[i] == vcoord[m] ) && (vpresson[m] == pressureOn || vpresson[m] == 0)) {
	add_event(i,3,m,vpropensity[m],j,-1);
	proball += vpropensity[m];
      }
    }
  }


  return proball;
}

/* ----------------------------------------------------------------------
   KMC method
   choose and perform an event for site
------------------------------------------------------------------------- */

void AppAldZns::site_event(int i, class RandomPark *random)
{
  int j,k,m,n,mm,jj;

  // periodic audit: catch event-list / solver-propensity desyncs before the
  // solver hits them (seen as "Site has no events" / "Corrupt event list" /
  // "Illegal execution event" crashes on the 500K run)
  if (++audit_counter >= AUDIT_INTERVAL) {
    audit_counter = 0;
    audit_consistency();
  }

  int elcoord = element[i];

  int isite2 = i2site[i];
  if (isite2 < 0 || firstevent[i] < 0 || propensity[isite2] <= 0.0) {
    dump_crash_context("Site has no events in site_event",i,-1,-1,-1,-1);
    heal_desync(i);
    return;
  }

  double threshhold = random->uniform() * propensity[isite2];
  double proball = 0.0;

  int ievent = firstevent[i];
  while (1) {
    proball += events[ievent].propensity;
    if (proball >= threshhold) break;
    ievent = events[ievent].next;
    if (ievent < 0) {
      dump_crash_context("Corrupt event list in site_event",i,-1,-1,-1,-1);
      heal_desync(i);
      return;
    }
  }

  int rstyle = events[ievent].style;
  int which = events[ievent].which;
  j = events[ievent].jpartner;
  k = events[ievent].kpartner;


  if (rstyle == 1) {
    element[i] = soutput[which];
    scount[which]++;
    } 
  else if (rstyle == 2 && j == -1) {
    element[i] = doutput[which][0];
    element[k] = doutput[which][1];
    dcount[which]++;
    }
  else if (rstyle == 3 && k == -1) {
    element[i] = voutput[which][0];
    element[j] = voutput[which][1];
    vcount[which]++;
    }
  else {
    dump_crash_context("Illegal execution event",i,rstyle,which,j,k);
    heal_desync(i);
    return;
  }

  // remember the event just executed, for crash diagnostics
  g_last_rstyle = rstyle;
  g_last_which = which;
  g_last_j = j;
  g_last_k = k;
  g_last_time = time;

  nmsites = 0;
  update_coord(elcoord,i,j,k,which);

  // sequence of ALD, 
  // 1 is metal pulse, 3 purge, 2 SULFUR pulse.
  if ((cycle+T1) > time ) {pressureOn = 1;}
  else if ((cycle+T1)<= time && time < (cycle+T1+T2)) {pressureOn = 3;}
  else if ((cycle+T1+T2) <= time && time < (cycle+T1+T2+T3)) {pressureOn = 2;}
  else if ((cycle+T1+T2+T3) <= time && time < (cycle+T1+T2+T3+T4)) {pressureOn = 3;}
  else {cycle += T1+T2+T3+T4; }

  

  int nsites = 0;
  int isite = i2site[i];
  
 

  propensity[isite] = site_propensity(i);
  esites[nsites++] = isite;
  echeck[isite] = 1;

  // go from site i to first and second neighbor in all type
  for (n = 0; n < numneigh[i]; n++) {
    m = neighbor[i][n];
    isite = i2site[m];
    if (isite >= 0 && echeck[isite] == 0) {
      propensity[isite] = site_propensity(m);
      esites[nsites++] = isite;
      echeck[isite] = 1;
    }
      for (jj = 0; jj< numneigh[m];jj++) {
        mm = neighbor[m][jj];
        isite = i2site[mm];
        if (isite >= 0 && echeck[isite] == 0) {
          propensity[isite] = site_propensity(mm);
          esites[nsites++] = isite;
          echeck[isite] = 1;
        }
     }
  }

  // go from site k to first and second neighbor in type II 
  if (rstyle == 2) {
    for (n = 0; n < numneigh[k]; n++) {
      m = neighbor[k][n];
      isite = i2site[m];
      if (isite >= 0 && echeck[isite] == 0) {
        propensity[isite] = site_propensity(m);
        esites[nsites++] = isite; 
        echeck[isite] = 1;
      }
      for (jj = 0; jj< numneigh[m];jj++) {
        mm = neighbor[m][jj];
        isite = i2site[mm];
        if (isite >= 0 && echeck[isite] == 0) {
          propensity[isite] = site_propensity(mm);
          esites[nsites++] = isite; 
          echeck[isite] = 1;
        }
      }
    }
  }

  // go from site j to first and second neighbor in type III 
  if (rstyle == 3) {
    for (n = 0; n < numneigh[j]; n++) {
      m = neighbor[j][n];
      isite = i2site[m];
      if (isite >= 0 && echeck[isite] == 0) {
        propensity[isite] = site_propensity(m);
        esites[nsites++] = isite;
        echeck[isite] = 1;
      }
      for (jj = 0; jj< numneigh[m];jj++) {
        mm = neighbor[m][jj];
        isite = i2site[mm];
        if (isite >= 0 && echeck[isite] == 0) {
          propensity[isite] = site_propensity(mm);
          esites[nsites++] = isite;
          echeck[isite] = 1;
        }
       }
    }
  }

  // merge sites whose coord was changed by put_mask/remove_mask in update_coord
  for (m = 0; m < nmsites; m++) {
    int mi = msites[m];
    int si = i2site[mi];
    if (si < 0) continue;
    if (echeck[si] == 0) {
      propensity[si] = site_propensity(mi);
      esites[nsites++] = si;
      echeck[si] = 1;
    }
  }

  solve->update(nsites,esites,propensity);
   // clear echeck array

  for (m = 0; m < nsites; m++)  {echeck[esites[m]] = 0; esites[m]=0;}
  nmsites = 0;
  
}

/* ----------------------------------------------------------------------
   clear all events out of list for site I
   add cleared events to free list
------------------------------------------------------------------------- */

void AppAldZns::clear_events(int i)
{
  int next;
  int index = firstevent[i];
  while (index >= 0) {
    next = events[index].next;
    events[index].next = freeevent;
    freeevent = index;
    nevents--;
    index = next;
  }
  firstevent[i] = -1;
}

/* ----------------------------------------------------------------------
   add an event to list for site I
   event = exchange with site J with probability = propensity
------------------------------------------------------------------------- */

void AppAldZns::add_event(int i, int rstyle, int which, double propensity,
			  int jpartner, int kpartner)
{
  if (nevents == maxevent) {
    maxevent += DELTAEVENT;
    events = 
      (Event *) memory->srealloc(events,maxevent*sizeof(Event),"app:events");
    for (int m = nevents; m < maxevent; m++) events[m].next = m+1;
    freeevent = nevents;
  }

  int next = events[freeevent].next;

  events[freeevent].style = rstyle;
  events[freeevent].which = which;
  events[freeevent].jpartner = jpartner;
  events[freeevent].kpartner = kpartner;
  events[freeevent].propensity = propensity;

  if ( propensity == 0 ) error->all(FLERR,"propensity in add_event wrong app ald");
  events[freeevent].next = firstevent[i];
  firstevent[i] = freeevent;
  freeevent = next;
  nevents++;
}

/* ----------------------------------------------------------------------
   grow list of stored reactions for single and double
------------------------------------------------------------------------- */

void AppAldZns::grow_reactions(int rstyle)
{
  if (rstyle == 1) {
    int n = none + 1;
    srate = (double *) 
      memory->srealloc(srate,n*sizeof(double),"app/ald:srate");
    spropensity = (double *) 
      memory->srealloc(spropensity,n*sizeof(double),"app/ald:spropensity");
    sinput = (int *) 
      memory->srealloc(sinput,n*sizeof(int),"app/ald:sinput");
    soutput = (int *) 
      memory->srealloc(soutput,n*sizeof(int),"app/ald:soutput");
    scount = (int *) 
      memory->srealloc(scount,n*sizeof(int),"app/ald:scount");
    sA = (double *) 
      memory->srealloc(sA,n*sizeof(double),"app/ald:sA");
    sexpon = (int *) 
      memory->srealloc(sexpon,n*sizeof(int),"app/ald:sexpon");
    scoord = (int *) 
      memory->srealloc(scoord,n*sizeof(int),"app/ald:scoord");
    spresson = (int *) 
      memory->srealloc(spresson,n*sizeof(int),"app/ald:spresson");

  } else if (rstyle == 2) {
    int n = ntwo + 1;
    drate = (double *) 
      memory->srealloc(drate,n*sizeof(double),"app/ald:drate");
    dpropensity = (double *) 
      memory->srealloc(dpropensity,n*sizeof(double),"app/ald:dpropensity");
    dinput = memory->grow(dinput,n,2,"app/ald:dinput");
    doutput = memory->grow(doutput,n,2,"app/ald:doutput");
    dcount = (int *) 
      memory->srealloc(dcount,n*sizeof(int),"app/ald:dcount");
    dA = (double *) 
      memory->srealloc(dA,n*sizeof(double),"app/ald:dA");
    dexpon = (int *) 
      memory->srealloc(dexpon,n*sizeof(int),"app/ald:dexpon");
    dcoord = (int *) 
      memory->srealloc(dcoord,n*sizeof(int),"app/ald:dcoord");
    dpresson = (int *) 
      memory->srealloc(dpresson,n*sizeof(int),"app/ald:dpresson");

  } else if (rstyle == 3) {
    int n = nthree + 1;
    vrate = (double *)
      memory->srealloc(vrate,n*sizeof(double),"app/ald:vrate");
    vpropensity = (double *)
      memory->srealloc(vpropensity,n*sizeof(double),"app/ald:vpropensity");
    vinput = memory->grow(vinput,n,2,"app/ald:vinput");
    voutput = memory->grow(voutput,n,2,"app/ald:voutput");
    vcount = (int *)
      memory->srealloc(vcount,n*sizeof(int),"app/ald:vcount");
    vA = (double *)
      memory->srealloc(vA,n*sizeof(double),"app/ald:vA");
    vexpon = (int *)
      memory->srealloc(vexpon,n*sizeof(int),"app/ald:vexpon");
    vcoord = (int *)
      memory->srealloc(vcoord,n*sizeof(int),"app/ald:vcoord");
    vpresson = (int *)
      memory->srealloc(vpresson,n*sizeof(int),"app/ald:vpresson");
  }
}

/* ----------------------------------------------------------------------
   update c.n. for Zn and S, put and remove mask for relative sites
------------------------------------------------------------------------- */
void AppAldZns::update_coord(int elcoord, int i, int j, int k, int which)
{
	if ((elcoord == S || elcoord == SH || elcoord == SH2) && (element[i] == ZnX2S || element[i] == ZnX2SH || element[i] == ZnX2SH2) && (j == -1 )) { // Adsorption of DEZ, event I
		coord[i]=coord[i]+1;
		put_mask(i);
	}
	else if ((elcoord == ZnX2S || elcoord == ZnX2SH || elcoord == ZnX2SH2) && (element[i] == S || element[i] == SH || element[i] == SH2) && (j == -1 )){ // Desorption of DEZ, event I
		coord[i]=coord[i]-1;
		remove_mask(i);
		count_coordS(i);
	}
    else if ((elcoord == ZnX2S || elcoord == ZnX2SH || elcoord == ZnX2SH2) && (element[i] == ZnXS || element[i] == ZnXSH ) && (element[k] == SH2 || element[k] == SH || element[k] == S ) && ( j == -1 ) && (k >= 0)){ // DEZ with H2S, event II
        remove_mask(i);
        put_mask(i);
    }
	else if ((j >= 0) && (elcoord == ZnX2S || elcoord == ZnX2SH || elcoord == ZnX2SH2) && (element[i] == ZnXS || element[i] == ZnXSH ) && ( element[j] == ZnX )){ // DEZ dissociation, event III
        remove_mask(i);
        put_mask(i);
		coord[j]=coord[j]+1;
		put_mask(j);
    }
    else if ((elcoord == SH2ZnX || elcoord == SHZnX ) && (element[i] == SH2Zn || element[i] == SHZn || element[i] == SZn) && (j == -1 ) ) { // MEZ with H2S, event I
        remove_mask(i);
        coord[i]=coord[i]-1;
    }
	else if ((j >= 0) && (elcoord == SH || elcoord == SH2) && (element[i] == SH || element[i] == S ) && ( element[j] == Zn ) ) { // MEZ with SH
        remove_mask(j);
        coord[j]=coord[j]-1;
    }
	else if ((j >= 0) && (elcoord == VACANCY) && (element[i] == ZnX || element[i] == Zn) && ( element[j] == S || element[j] == SH || element[j] == SH2 )  ) { // Zn densification
        if (element[i] == ZnX){
        	remove_mask(i, j); // Remove mask from previous S site
        }
        count_coord(j, i);
        if (element[i] == ZnX){
            put_mask(i); // Put mask on new Zn site
        }
	}
	else if ((j >= 0) && (elcoord == ZnX ) && (element[i] == VACANCY ) && ( element[j] == ZnXSH || element[j] == ZnXS) ) { // ZnX reverse densification
		remove_mask(i, j);
		count_coord(i, j);
		put_mask(j);
	}
	else if ((j >= 0) && (elcoord == Zn ) && (element[i] == VACANCY ) && ( element[j] == ZnSH || element[j] == ZnS) ) { // Zn reverse densification 
        count_coord(j, i);
        }
    else if ((j >= 0) && (elcoord == SH2Zn || elcoord == SHZn || elcoord == SZn) && ( element[i] == Zn ) && (element[j] == S || element[j] == SH || element[j] == SH2 )) { // SULFUR densification
        count_coord(i, j); // Reversed from normal ordering to avoid mixup in count_coord-function that expects i-> S site, j-> Zn site
        count_coordS(j);
    }
    else if ((j >= 0) && (elcoord == SH2ZnX || elcoord == SHZnX ) && ( element[i] == ZnX ) && (element[j] == SH || element[j] == SH2 ) && ( k == -1 )) { // SULFUR densification
        count_coord(i, j); // Reversed from normal ordering to avoid mixup in count_coord-function that expects i-> S site, j-> Zn site
    }
    else if ((elcoord == Zn || elcoord == ZnX ) && (element[i] == SH2Zn || element[i] == SH2ZnX ) && ( j == -1 )) { // Adsorption of H2S
        coord[i]=coord[i]+1;
    }
    else if ((elcoord == SH2Zn || elcoord == SH2ZnX ) && (element[i] == Zn || element[i] == ZnX ) && ( j == -1 )) { // Desorption of H2S
        coord[i]=coord[i]-1;
    }
    else if ((j >= 0) && (elcoord == ZnX ) && ( element[i] == ZnX && element[j] == VACANCY ) ) { // Desorption of H2S
        count_coord(j, i);
    }
    else if ((elcoord == SH2 ) && (element[i] == VACANCY ) && ( j == -1 )){ // Desorption of H2S
        count_coord(i, j);
    }
	else if ((j >= 0) && (elcoord == S || elcoord == SH || elcoord == SH2) && (element[i] == VACANCY) && (element[j] == SH2ZnX || element[j] == SHZnX || element[j] == SH2Zn || element[j] == SHZn || element[j] == SZn ) ) {// SULFUR reverse densification
		count_coord(i,j);
		coord[j]++;
	}
}
/* ----------------------------------------------------------------------
   put mask for affected sites
------------------------------------------------------------------------- */

void AppAldZns::put_mask(int i)
{
    int isite = i2site[i];
        if (isite < 0) return;
	int nsites = 0;
	esites[nsites++] = isite;
    echeck[isite] = 1;
// Add mask on the second neighbor (SULFUR) of the DEZ to block adsorption
	if (element[i] == ZnX2SH2 || element[i] == ZnX2SH || element[i] == ZnX2S ){
	  	for (int n = 0; n < numneigh[i]; n++) {
			int nn = neighbor[i][n];
			isite = i2site[nn];
			if (isite >= 0 && echeck[isite] == 0) { // Cover first neighbour Zn site
			    coord[isite]=coord[isite]-20; track_msite(nn);
			    esites[nsites++] = isite;
			    echeck[isite] = 1;
                        }
			for (int k = 0; k < numneigh[nn]; k++){
				int kk = neighbor[nn][k];
				isite = i2site[kk];
				if (isite >= 0 && echeck[isite] == 0) { // Cover second neighbour S site
					coord[isite]=coord[isite]-10; track_msite(kk);
					esites[nsites++] = isite;
					echeck[isite] = 1;
				}
				for (int m = 0; m < numneigh[kk]; m++) {
					int mm = neighbor[kk][m];
					isite = i2site[mm];
					if (isite >= 0 && echeck[isite] == 0) { // Cover third neighbour Zn site
					    coord[isite] = coord[isite]-10; track_msite(mm);
					    esites[nsites++] = isite;
					    echeck[isite] = 1;
                                        }
					for (int s = 0; s < numneigh[mm]; s++) {
						int ss = neighbor[mm][s];
						isite = i2site[ss];
						if (isite >= 0 && echeck[isite] == 0) { // Cover fourth neighbour S site
							coord[isite] = coord[isite]-10; track_msite(ss);
							esites[nsites++] = isite;
							echeck[isite] = 1;
						}
					}
				}
			}
	    }
    }

// Add mask to the second neighbor (SULFUR) of the DEZ to block adsorption (ZnX reverse densification)
    else if ( element[i] == ZnXSH || element[i] == ZnXS ){
        for (int n = 0; n < numneigh[i]; n++) {
            int nn = neighbor[i][n];
            isite = i2site[nn];
            if (isite >= 0 && echeck[isite] == 0) { // Cover first neighbour Zn site
                coord[isite]=coord[isite]-10; track_msite(nn);
                esites[nsites++] = isite;
                echeck[isite] = 1;
            }
            for (int k = 0; k < numneigh[nn]; k++){
                int kk = neighbor[nn][k];
                isite = i2site[kk];
                if (isite >= 0 && echeck[isite] == 0) { // Cover second neighbour S site
                    coord[isite]=coord[isite]-10; track_msite(kk);
                    esites[nsites++] = isite;
                    echeck[isite] = 1;
                }
                for (int m = 0; m < numneigh[kk]; m++) {
                    int mm = neighbor[kk][m];
                    isite = i2site[mm];
                    if (isite >= 0 && echeck[isite] == 0) { // Cover third neighbour Zn site
                        coord[isite] = coord[isite]-10; track_msite(mm);
                        esites[nsites++] = isite;
                        echeck[isite] = 1;
                    }
                    for (int s = 0; s < numneigh[mm]; s++) {
                        int ss = neighbor[mm][s];
                        isite = i2site[ss];
                        if (isite >= 0 && echeck[isite] == 0) { // Cover fourth neighbour S site
                            coord[isite] = coord[isite]-10; track_msite(ss);
                            esites[nsites++] = isite;
                            echeck[isite] = 1;
                        }
                    }
                }
            }
        }
    }	


// Add mask to the first neighbor (SULFUR) to block adsorption
	else if ( element[i] == ZnX ){
	  	for (int n = 0; n < numneigh[i]; n++) {
			int nn = neighbor[i][n];
			isite = i2site[nn];
			if (isite >= 0 && echeck[isite] == 0) { // Cover first neighbour S site
				coord[isite]=coord[isite]-10; track_msite(nn);
				esites[nsites++] = isite;
				echeck[isite] = 1;
			}
            for (int k = 0; k < numneigh[nn]; k++){
                int kk = neighbor[nn][k];
                isite = i2site[kk];
                if (isite >= 0 && echeck[isite] == 0) { // Cover second neighbour Zn site
                    coord[isite] = coord[isite]-10; track_msite(kk);
                    esites[nsites++] = isite;
                    echeck[isite] = 1;
                }
                for (int m = 0; m < numneigh[kk]; m++) {
                    int mm = neighbor[kk][m];
                    isite = i2site[mm];
                    if (isite >= 0 && echeck[isite] == 0) { // Cover third neighbour Zn site
                        coord[isite] = coord[isite]-10; track_msite(mm);
                        esites[nsites++] = isite;
                        echeck[isite] = 1;
                    }
                    for (int s = 0; s < numneigh[mm]; s++) {
                        int ss = neighbor[mm][s];
                        isite = i2site[ss];
                        if (isite >= 0 && echeck[isite] == 0) { // Cover fourth neighbour S site
                            coord[isite] = coord[isite]-10; track_msite(ss);
                            esites[nsites++] = isite;
                            echeck[isite] = 1;
                        }
                    }
                }
            }
        }
	}
        for (int m = 0; m < nsites; m++)  {echeck[esites[m]] = 0; esites[m]=0;}
}

/* ----------------------------------------------------------------------
   remove mask 
------------------------------------------------------------------------- */

void AppAldZns::remove_mask(int i, int j) // j flag for when Zn densification
{
    int isite = i2site[i];
        if (isite < 0) return;
	int nsites = 0;
	esites[nsites++] = isite;
    echeck[isite] = 1;
// Remove mask from SULFUR sites after desorption
	if ( element[i] == S || element[i] == SH || element[i] == SH2 || element[i] == ZnXS || element[i] == ZnXSH ){
	  	for (int n = 0; n < numneigh[i]; n++) {
			int nn = neighbor[i][n];
			isite = i2site[nn];
			if (isite >= 0 && echeck[isite] == 0) { // Remove first neighbour Zn site
                coord[isite]=coord[isite]+20; track_msite(nn);
                esites[nsites++] = isite;
                echeck[isite] = 1;
            }
			for (int k = 0; k < numneigh[nn]; k++){
				int kk = neighbor[nn][k];
				isite = i2site[kk];
				if (isite >= 0 && echeck[isite] == 0) { // Remove second neighbour S site
					coord[isite]=coord[isite]+10; track_msite(kk);
					esites[nsites++] = isite;
					echeck[isite] = 1;
				}
				for (int m = 0; m < numneigh[kk]; m++) {
					int mm = neighbor[kk][m];
					isite = i2site[mm];
					if (isite >= 0 && echeck[isite] == 0) {// Remove third neighbour Zn site
					    coord[isite] = coord[isite]+10; track_msite(mm);
					    esites[nsites++] = isite;
					    echeck[isite] = 1;
                    }
					for (int s = 0; s < numneigh[mm]; s++) {
						int ss = neighbor[mm][s];
						isite = i2site[ss];
						if (isite >= 0 && echeck[isite] == 0) {// Remove fourth neighbour S site
							coord[isite] = coord[isite]+10; track_msite(ss);
							esites[nsites++] = isite;
							echeck[isite] = 1;
						}
					}
				}
			}
        }
	}
	
// Remove mask from the SULFUR site after densification
	else if ( ( element[i] == ZnX && ( element[j] == S || element[j] == SH || element[j] == SH2 )) ){ 
	    echeck[i2site[i]] = 0;
	    for (int n = 0; n < numneigh[j]; n++) {
	        int nn = neighbor[j][n];
	        isite = i2site[nn];
	        if (isite >= 0 && echeck[isite] == 0 ) { // Remove first neighbour Zn site
	            coord[isite]=coord[isite]+10; track_msite(nn);
                esites[nsites++] = isite;
                echeck[isite] = 1;
            }
            for (int k = 0; k < numneigh[nn]; k++){
                int kk = neighbor[nn][k];
                if(kk != j){
                    isite = i2site[kk];
                    if (isite >= 0 && echeck[isite] == 0) { // Remove second neighbour S site
                        if(isite!=j){coord[isite] = coord[isite]+10; track_msite(kk);}
                        esites[nsites++] = isite;
                        echeck[isite] = 1;
                    }
                }
                for (int m = 0; m < numneigh[kk]; m++) {
					int mm = neighbor[kk][m];
					isite = i2site[mm];
					if (isite >= 0 && echeck[isite] == 0) { // Cover third neighbour Zn site
						coord[isite] = coord[isite]+10; track_msite(mm);
						esites[nsites++] = isite;
						echeck[isite] = 1;
					}
					for (int s = 0; s < numneigh[mm]; s++) {
						int ss = neighbor[mm][s];
						isite = i2site[ss];
						if (isite >= 0 && echeck[isite] == 0) { // Cover fourth neighbour S site
							coord[isite] = coord[isite]+10; track_msite(ss);
							esites[nsites++] = isite;
							echeck[isite] = 1;
						}
					}
				}
            }
        }
    }
	    
	
// Remove mask after second ligand has been removed
	else if ( element[i] == SZn || element[i] == SHZn || element[i] == SH2Zn ||  element[i] == ZnSH || element[i] == ZnS || element[i] == Zn ){
	  	for (int n = 0; n < numneigh[i]; n++) {
      	  	int nn = neighbor[i][n];
            isite = i2site[nn];
            if (isite >= 0 && echeck[isite] == 0 ) {
                    coord[isite]=coord[isite]+10; track_msite(nn);
                    esites[nsites++] = isite;
                    echeck[isite] = 1;
            }
            for (int k = 0; k < numneigh[nn]; k++){
                int kk = neighbor[nn][k];
                isite = i2site[kk];
                if (isite >= 0 && echeck[isite] == 0) {
                    coord[isite] = coord[isite]+10; track_msite(kk);
                    esites[nsites++] = isite;
                    echeck[isite] = 1;
                }
				for (int m = 0; m < numneigh[kk]; m++) {
					int mm = neighbor[kk][m];
					isite = i2site[mm];
					if (isite >= 0 && echeck[isite] == 0) { // Cover third neighbour Zn site
						coord[isite] = coord[isite]+10; track_msite(mm);
						esites[nsites++] = isite;
                        echeck[isite] = 1;
                    }
                    for (int s = 0; s < numneigh[mm]; s++) {
                        int ss = neighbor[mm][s];
                        isite = i2site[ss];
                        if (isite >= 0 && echeck[isite] == 0) { // Cover fourth neighbour S site
                            coord[isite] = coord[isite]+10; track_msite(ss);
                            esites[nsites++] = isite;
                            echeck[isite] = 1;
                        }
                    }
                }
            }
        }
	}
// Remove mask after reverse densification
    else if ( element[i] == VACANCY && ( element[j] == ZnXSH || element[j] == ZnXS )){
        for (int n = 0; n < numneigh[i]; n++) {
            int nn = neighbor[i][n];
            isite = i2site[nn];
            if (isite >= 0 && echeck[isite] == 0 ) {
                coord[isite]=coord[isite]+10; track_msite(nn);
                esites[nsites++] = isite;
                echeck[isite] = 1;
            }
            for (int k = 0; k < numneigh[nn]; k++){
                int kk = neighbor[nn][k];
                isite = i2site[kk];
                if (isite >= 0 && echeck[isite] == 0) {
                    coord[isite] = coord[isite]+10; track_msite(kk);
                    esites[nsites++] = isite;
                    echeck[isite] = 1;
                }
				for (int m = 0; m < numneigh[kk]; m++) {
					int mm = neighbor[kk][m];
					isite = i2site[mm];
					if (isite >= 0 && echeck[isite] == 0) { // Cover third neighbour Zn site
                        coord[isite] = coord[isite]+10; track_msite(mm);
                        esites[nsites++] = isite;
                        echeck[isite] = 1;
                    }
                    for (int s = 0; s < numneigh[mm]; s++) {
						int ss = neighbor[mm][s];
						isite = i2site[ss];
						if (isite >= 0 && echeck[isite] == 0) { // Cover fourth neighbour S site
                            coord[isite] = coord[isite]+10; track_msite(ss);
                            esites[nsites++] = isite;
                            echeck[isite] = 1;
                        }
                    }
                }
            }
        }
    }
      
	  	for (int m = 0; m < nsites; m++)  {echeck[esites[m]] = 0; esites[m]=0;}
}

/* ----------------------------------------------------------------------
   count c.n after densification
------------------------------------------------------------------------- */


void AppAldZns::count_coord(int i, int j) // i: SULFUR species, j: Zinc species (does not necessarily hold)
{
    if (i < 0) return;
// Densification of ZnXSH, ZnXS -> ZnX
    if (j >= 0 && ( element[i] == S || element[i] == SH || element[i] == SH2 ) &&  ( element[j] == ZnX )  ){
			coord[j]=coord[j]+1; // Add one because of X ligand
			for (int s = 0; s < numneigh[j]; s++){
				int nn = neighbor[j][s];
				if (element[nn] >= S && element[nn] <= ZnSH ) { // Check if neighbouring site is an SULFUR site
					coord[j]=coord[j] + 1;
					if (i != nn){ coord[nn]=coord[nn]+1; } // Careful not to change the cn of original site
				}
			}
    }
// Densification of ZnSH, ZnS -> Zn
    else if (j >= 0 && ( element[i] == S || element[i] == SH || element[i] == SH2 ) &&  ( element[j] == Zn )  ){
        for (int s = 0; s < numneigh[j]; s++){
                int nn = neighbor[j][s];
                if (element[nn] >= S && element[nn] <= ZnSH ) { // Check if neighbouring site is an SULFUR site
                        coord[j]=coord[j] + 1;
                        if (i != nn){ coord[nn]=coord[nn]+1; }
                }
        }
    }
// Densification of SULFUR-species
    else if (j >= 0 && ( element[i] == ZnX || element[i] == Zn ) &&  ( element[j] == S || element[j] == SH || element[j] == SH2 )  ){
        for (int s = 0; s < numneigh[j]; s++){
                int nn = neighbor[j][s];
                if ( Zn <= element[nn] && element[nn] <= SZn) { // Check if neighbouring site is an zinc site
                        coord[j]=coord[j] + 1;
                        if (i != nn){ coord[nn]=coord[nn]+1; 
//                            if(element[nn] == Zn and coord[nn] > 4){printf("i: %d %d %d j: %d %d %d nn: %d %d %d\n", i, element[i], coord[i], j, element[j], coord[j], nn, element[nn], coord[nn]);}
                        }
                }
        }
    }
// Reverse densification on ZnX
    else if (j >= 0 &&  element[i] == VACANCY  && (element[j] >= ZnXS && element[j] <= ZnSH ) ){
        if (j >= 0 &&  element[j] == ZnXS || element[j] == ZnXSH ){ coord[i]=coord[i] - 1; } // Remove the extra cn from ligand
        for (int s = 0; s < numneigh[i]; s++){
            int nn = neighbor[i][s];
            if ( element[nn] >= S && element[nn] <= ZnSH ) { // Check if neighbouring site is an SULFUR site
                coord[i]=coord[i] - 1;
                if (j != nn){ coord[nn]=coord[nn] - 1;
                }
            }
        }
    }
// Desorption of SH2, event 1 and 3
    else if ( element[i] == VACANCY  && ( j == -1 || element[j] == ZnX ) ){
        for (int s = 0; s < numneigh[i]; s++){
            int nn = neighbor[i][s];
            if ( element[nn] >= Zn && element[nn] <= SZn ) { // Check if neighbouring site is a zinc site
                coord[i]=coord[i] - 1;
                if (i != nn){ coord[nn]=coord[nn] - 1;}
            }
        }
    }
// Reverse densification of SH2 / SH / S
    else if (j >= 0 &&  element[i] == VACANCY  && ( element[j] == SH2ZnX || element[j] == SH2Zn || element[j] == SHZnX || element[j] == SHZn || element[j] == SZn) ){
        for (int s = 0; s < numneigh[i]; s++){
            int nn = neighbor[i][s];
            if ( element[nn] >= Zn && element[nn] <= SZn ) { // Check if neighbouring site is a zinc site
                coord[i]=coord[i] - 1;
                if (i != nn){ coord[nn]=coord[nn] - 1;}
            }
        }
    }

}

/* ----------------------------------------------------------------------
   remember a lattice site whose coord was changed by put_mask/remove_mask,
   so site_event can refresh its solver propensity (dedup by linear scan)
------------------------------------------------------------------------- */
void AppAldZns::track_msite(int ilat)
{
  for (int q = 0; q < nmsites; q++)
    if (msites[q] == ilat) return;
  if (nmsites >= nlocal)
    error->all(FLERR,"msites overflow in app_ald");
  msites[nmsites++] = ilat;
}

/* ----------------------------------------------------------------------
   count c.n of SULFUR before adsorption
------------------------------------------------------------------------- */
void AppAldZns::count_coordS(int i)
{
    int fullS = 0;
    int emptyS = 0;
    int totalS = 0;

    int isite = i2site[i];
    int nsites = 0;

	for (int m = 0; m < numneigh[i]; m++) {
		int mm = neighbor[i][m];
		for (int s = 0; s < numneigh[mm]; s++) {
			int ss = neighbor[mm][s];
			isite = i2site[ss];
			if (i==ss)  continue;
			if (isite >= 0 && echeck[isite] == 0) {
			  if ( element[ss] >= S && element[ss] <= ZnSH ) {fullS++;}
			  else if (element[ss] == VACANCY) {emptyS++;}
		          esites[nsites++] = isite;
		          echeck[isite] = 1;
		        }
    
		}
	}
   totalS = fullS+emptyS;
   if ( float(fullS) > 4*totalS/5 and coord[i] > -20) {coord[i] += -20;} // decrease the coord of the SULFUR site to render it inactive for adsorption
   for (int m = 0; m < nsites; m++)  {echeck[esites[m]] = 0; esites[m]=0;}
}

/* ----------------------------------------------------------------------
   periodic consistency audit
   finds sites whose solver propensity is > 0 while their event list is
   empty (or vice versa) and repairs them, so a single corrupt site cannot
   abort the run with "Site has no events" / "Corrupt event list"
------------------------------------------------------------------------- */

void AppAldZns::audit_consistency()
{
  int nbad1 = 0, nbad2 = 0, nbad3 = 0, nheal = 0;

  for (int x = 0; x < nlocal; x++) {
    int sx = i2site[x];
    if (sx < 0) continue;
    double p = propensity[sx];
    int fe = firstevent[x];
    if (p > 0.0 && fe < 0) nbad1++;
    else if (fe >= 0 && p == 0.0) nbad2++;
    if (fe >= maxevent) nbad3++;
  }

  if (nbad1 == 0 && nbad2 == 0 && nbad3 == 0) return;

  if (screen)
    fprintf(screen,
      "APPALD-AUDIT t=%.10g: %d prob>0/no-events, %d events/prob=0, %d bad-firstevent; repairing\n",
      time,nbad1,nbad2,nbad3);
  if (logfile)
    fprintf(logfile,
      "APPALD-AUDIT t=%.10g: %d prob>0/no-events, %d events/prob=0, %d bad-firstevent; repairing\n",
      time,nbad1,nbad2,nbad3);

  for (int x = 0; x < nlocal; x++) {
    int sx = i2site[x];
    if (sx < 0) continue;
    double p = propensity[sx];
    int fe = firstevent[x];
    if ((p > 0.0 && fe < 0) || fe >= maxevent) {
      heal_desync(x);
      nheal++;
    }
  }

  if (screen) fprintf(screen,"APPALD-AUDIT: repaired %d sites\n",nheal);
  if (logfile) fprintf(logfile,"APPALD-AUDIT: repaired %d sites\n",nheal);
}

/* ----------------------------------------------------------------------
   dump full context when the event machinery detects an inconsistent site
------------------------------------------------------------------------- */

void AppAldZns::dump_crash_context(const char *msg, int i,
                                int rstyle, int which, int jp, int kp)
{
  if (screen) {
    fprintf(screen,"\n=== AppALD event-list diagnostic: %s ===\n",msg);
    fprintf(screen,"time=%.10g pressureOn=%d cycle=%.10g T1=%.10g T2=%.10g T3=%.10g T4=%.10g\n",
            time,pressureOn,cycle,T1,T2,T3,T4);
    fprintf(screen,"last executed event: rstyle=%d which=%d j=%d k=%d at t=%.10g\n",
            g_last_rstyle,g_last_which,g_last_j,g_last_k,g_last_time);
    int isite = (i >= 0 && i < nlocal) ? i2site[i] : -1;
    fprintf(screen,"picked site i=%d i2site=%d element=%d coord=%d firstevent=%d",
            i,isite,(i>=0&&i<nlocal)?element[i]:-999,(i>=0&&i<nlocal)?coord[i]:-999,
            (i>=0&&i<nlocal)?firstevent[i]:-999);
    if (isite >= 0 && isite < nlocal)
      fprintf(screen," solver-prob=%.6g",propensity[isite]);
    fprintf(screen,"\n");
    if (rstyle >= 0)
      fprintf(screen,"picked event rstyle=%d which=%d j=%d k=%d\n",rstyle,which,jp,kp);
    fprintf(screen,"nevents=%d maxevent=%d solver-sum=%.6g num-active=%d\n",
            nevents,maxevent,
            solve ? solve->get_total_propensity() : -1.0,
            solve ? solve->get_num_active() : -1);
    int nbad1 = 0, nbad2 = 0, nbad3 = 0;
    for (int x = 0; x < nlocal; x++) {
      int sx = i2site[x];
      if (sx < 0) continue;
      double p = propensity[sx];
      int fe = firstevent[x];
      if (p > 0.0 && fe < 0) {
        if (nbad1 < 10)
          fprintf(screen,"  [prob>0,no-events] site=%d el=%d coord=%d prob=%.6g\n",
                  x,element[x],coord[x],p);
        nbad1++;
      }
      if (fe >= 0 && p == 0.0) {
        if (nbad2 < 10)
          fprintf(screen,"  [events,prob=0] site=%d el=%d coord=%d fe=%d\n",
                  x,element[x],coord[x],fe);
        nbad2++;
      }
      if (fe >= maxevent) {
        if (nbad3 < 10)
          fprintf(screen,"  [bad-firstevent] site=%d fe=%d maxevent=%d\n",
                  x,fe,maxevent);
        nbad3++;
      }
    }
    fprintf(screen,"summary: prob>0/no-events=%d events/prob=0=%d bad-firstevent=%d\n",
            nbad1,nbad2,nbad3);
  }
}

/* ----------------------------------------------------------------------
   repair an inconsistent site: free its (possibly corrupt) event list with
   a bounded walk, rebuild events via site_propensity, refresh solver entry
   returns 1 if the site could not be repaired (skip the pick)
------------------------------------------------------------------------- */

int AppAldZns::heal_desync(int i)
{
  if (++heal_count > MAX_HEALS) {
    dump_crash_context("too many event-list repairs, aborting",i,-1,-1,-1,-1);
    error->all(FLERR,"Excessive event-list corruption in app_ald_zns");
  }

  if (i < 0 || i >= nlocal) return 1;
  int isite = i2site[i];
  if (isite < 0) return 1;

  double p0 = propensity[isite];
  free_events_safe(i);
  propensity[isite] = site_propensity(i);
  int idx = isite;
  solve->update(1,&idx,propensity);

  if (screen)
    fprintf(screen,
      "APPALD-REPAIR t=%.10g site=%d el=%d coord=%d prob %.6g -> %.6g firstevent %d (heal #%d)\n",
      time,i,element[i],coord[i],p0,propensity[isite],firstevent[i],heal_count);
  if (logfile)
    fprintf(logfile,
      "APPALD-REPAIR t=%.10g site=%d el=%d coord=%d prob %.6g -> %.6g firstevent %d (heal #%d)\n",
      time,i,element[i],coord[i],p0,propensity[isite],firstevent[i],heal_count);
  return 0;
}

/* ----------------------------------------------------------------------
   bounded version of clear_events: a corrupted (cyclic) list must not
   hang the run, so stop after maxevent+1 steps and drop the tail
------------------------------------------------------------------------- */

void AppAldZns::free_events_safe(int i)
{
  int count = 0;
  int index = firstevent[i];
  while (index >= 0 && index < maxevent && count <= maxevent) {
    int next = events[index].next;
    events[index].next = freeevent;
    freeevent = index;
    if (nevents > 0) nevents--;
    index = next;
    count++;
  }
  firstevent[i] = -1;
}
