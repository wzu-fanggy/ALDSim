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

/* ----------------------------------------------------------------------
   ALD application of HfO2 was developed by 
   mahdi shirazi: m.shirazi@tue.nl, TU/e department of applied physics,
   Simon D. Elliott: simon.elliott@schrodinger.com, Schrodinger Materials Science.
   This application is a part of SPPARKS and authors retian the above term.
   See the manual-app-ald and examples folders for more information.
------------------------------------------------------------------------- */

#include "math.h"
#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "app_ald.h"
#include "solve.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS_NS;
enum{VACANCY,O,OH,//2
HfX4O,HfX4OH,HfHX4O,HfHX4OH,HfH2X4O,HfH2X4OH,HfH3X4O,HfH3X4OH,HfH4X4O,HfH4X4OH,//12
HfX3O,HfX3OH,HfHX3O,HfHX3OH,HfH2X3O,HfH2X3OH,HfH3X3O,HfH3X3OH,//20
HfX2O,HfX2OH,HfHX2O,HfHX2OH,HfH2X2O,HfH2X2OH,//26
HfX2,HfHX2,HfH2X2,//29
HfHX,HfX,Hf,//32
OH2HfX,OH2HfHX,OH2Hf,OHHfHX,//36
AlX3O,AlX3OH,AlHX2O,AlHX2OH,AlH2XO,AlH2XOH,//42
AlX2O,AlX2OH,AlHXO,AlHXOH,AlXO,AlXOH,AlOH,//49
AlX2,AlHX2,AlH2X2,//52
AlHX,AlX,Al,//55
OH2AlX,OH2AlHX,OH2Al,//58
OH2,Si};//60 both Hf and Al



/* ---------------------------------------------------------------------- */
/* map an input species name to its integer id, -1 if unknown             */

static int species_id(const char *name)
{
  if (strcmp(name,"VAC") == 0) return VACANCY;
  if (strcmp(name,"O") == 0) return O;
  if (strcmp(name,"OH") == 0) return OH;
  if (strcmp(name,"OH2") == 0) return OH2;
  if (strcmp(name,"AlX3O") == 0) return AlX3O;
  if (strcmp(name,"AlX3OH") == 0) return AlX3OH;
  if (strcmp(name,"AlHX2O") == 0) return AlHX2O;
  if (strcmp(name,"AlHX2OH") == 0) return AlHX2OH;
  if (strcmp(name,"AlH2XO") == 0) return AlH2XO;
  if (strcmp(name,"AlH2XOH") == 0) return AlH2XOH;
  if (strcmp(name,"AlX2O") == 0) return AlX2O;
  if (strcmp(name,"AlX2OH") == 0) return AlX2OH;
  if (strcmp(name,"AlHXO") == 0) return AlHXO;
  if (strcmp(name,"AlHXOH") == 0) return AlHXOH;
  if (strcmp(name,"AlXO") == 0) return AlXO;
  if (strcmp(name,"AlXOH") == 0) return AlXOH;
  if (strcmp(name,"AlOH") == 0) return AlOH;
  if (strcmp(name,"AlX2") == 0) return AlX2;
  if (strcmp(name,"AlHX2") == 0) return AlHX2;
  if (strcmp(name,"AlH2X2") == 0) return AlH2X2;
  if (strcmp(name,"AlHX") == 0) return AlHX;
  if (strcmp(name,"AlX") == 0) return AlX;
  if (strcmp(name,"Al") == 0) return Al;
  if (strcmp(name,"HfX4O") == 0) return HfX4O;
  if (strcmp(name,"HfX4OH") == 0) return HfX4OH;
  if (strcmp(name,"HfHX4O") == 0) return HfHX4O;
  if (strcmp(name,"HfHX4OH") == 0) return HfHX4OH;
  if (strcmp(name,"HfH2X4O") == 0) return HfH2X4O;
  if (strcmp(name,"HfH2X4OH") == 0) return HfH2X4OH;
  if (strcmp(name,"HfH3X4O") == 0) return HfH3X4O;
  if (strcmp(name,"HfH3X4OH") == 0) return HfH3X4OH;
  if (strcmp(name,"HfH4X4O") == 0) return HfH4X4O;
  if (strcmp(name,"HfH4X4OH") == 0) return HfH4X4OH;
  if (strcmp(name,"HfX3O") == 0) return HfX3O;
  if (strcmp(name,"HfX3OH") == 0) return HfX3OH;
  if (strcmp(name,"HfHX3O") == 0) return HfHX3O;
  if (strcmp(name,"HfHX3OH") == 0) return HfHX3OH;
  if (strcmp(name,"HfH2X3O") == 0) return HfH2X3O;
  if (strcmp(name,"HfH2X3OH") == 0) return HfH2X3OH;
  if (strcmp(name,"HfH3X3O") == 0) return HfH3X3O;
  if (strcmp(name,"HfH3X3OH") == 0) return HfH3X3OH;
  if (strcmp(name,"HfX2O") == 0) return HfX2O;
  if (strcmp(name,"HfX2OH") == 0) return HfX2OH;
  if (strcmp(name,"HfHX2O") == 0) return HfHX2O;
  if (strcmp(name,"HfHX2OH") == 0) return HfHX2OH;
  if (strcmp(name,"HfH2X2O") == 0) return HfH2X2O;
  if (strcmp(name,"HfH2X2OH") == 0) return HfH2X2OH;
  if (strcmp(name,"HfX2") == 0) return HfX2;
  if (strcmp(name,"HfHX2") == 0) return HfHX2;
  if (strcmp(name,"HfH2X2") == 0) return HfH2X2;
  if (strcmp(name,"HfHX") == 0) return HfHX;
  if (strcmp(name,"HfX") == 0) return HfX;
  if (strcmp(name,"Hf") == 0) return Hf;
  if (strcmp(name,"OH2HfX") == 0) return OH2HfX;
  if (strcmp(name,"OH2HfHX") == 0) return OH2HfHX;
  if (strcmp(name,"OH2Hf") == 0) return OH2Hf;
  if (strcmp(name,"OHHfHX") == 0) return OHHfHX;
  if (strcmp(name,"OH2AlX") == 0) return OH2AlX;
  if (strcmp(name,"OH2AlHX") == 0) return OH2AlHX;
  if (strcmp(name,"OH2Al") == 0) return OH2Al;
  if (strcmp(name,"Si") == 0) return Si;
  return -1;
}

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

AppAld::AppAld(SPPARKS *spk, int narg, char **arg) : 
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
  scount = dcount = vcount = NULL;
  sA = dA = vA = NULL;
  scoord = dcoord = vcoord = NULL;
  sexpon = dexpon = vexpon = NULL;
  spresson = dpresson = vpresson = NULL;
}

/* ---------------------------------------------------------------------- */

AppAld::~AppAld()
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

void AppAld::input_app(char *command, int narg, char **arg)
{
  if (strcmp(command,"event") == 0) {
    if (narg < 1) error->all(FLERR,"Illegal event command");
    int rstyle = atoi(arg[0]);
    grow_reactions(rstyle);

    if (rstyle == 1) {
      if (narg != 9) error->all(FLERR,"Illegal event arg command");
//type I
      //type I: single site, from -> to species
      int sid = species_id(arg[1]);
      if (sid < 0) error->all(FLERR,"Illegal event command");
      sinput[none] = sid;
      sid = species_id(arg[2]);
      if (sid < 0) error->all(FLERR,"Illegal event command");
      soutput[none] = sid;
            sA[none] = atof(arg[3]);
      if (sA[none] == 0.0) error->warning(FLERR,"Illegal coef during reading command");
      sexpon[none] = atoi(arg[4]);
      srate[none] = atof(arg[5]);
      scoord[none] = atoi(arg[6]);
      spresson[none] = atoi(arg[7]);

      none++;
      
//type II 
    } else if (rstyle == 2) {
      if (narg != 11) error->all(FLERR,"Illegal event command");

      //type II: site i and 2nd-neighbor k, both change species
      int sid = species_id(arg[1]);
      if (sid < 0) error->all(FLERR,"Illegal event command");
      dinput[ntwo][0] = sid;
      sid = species_id(arg[2]);
      if (sid < 0) error->all(FLERR,"Illegal event command");
      doutput[ntwo][0] = sid;
      sid = species_id(arg[3]);
      if (sid < 0) error->all(FLERR,"Illegal event command");
      dinput[ntwo][1] = sid;
      sid = species_id(arg[4]);
      if (sid < 0) error->all(FLERR,"Illegal event command");
      doutput[ntwo][1] = sid;
            dA[ntwo] = atof(arg[5]);
      dexpon[ntwo] = atoi(arg[6]);
      if (dexpon[ntwo] != 0.0) error->warning(FLERR,"Illegal expon command2");
      drate[ntwo] = atof(arg[7]);
      dcoord[ntwo] = atoi(arg[8]);
      dpresson[ntwo] = atoi(arg[9]);
      ntwo++;

    }else if (rstyle == 3) {
      if (narg != 11) error->all(FLERR,"Illegal event command31");
 
      //type III: site i and 1st-neighbor j, both change species
      int sid = species_id(arg[1]);
      if (sid < 0) error->all(FLERR,"Illegal event command32");
      vinput[nthree][0] = sid;
      sid = species_id(arg[2]);
      if (sid < 0) error->all(FLERR,"Illegal event command33");
      voutput[nthree][0] = sid;
      sid = species_id(arg[3]);
      if (sid < 0) error->all(FLERR,"Illegal event command34");
      vinput[nthree][1] = sid;
      sid = species_id(arg[4]);
      if (sid < 0) error->all(FLERR,"Illegal event command35");
      voutput[nthree][1] = sid;
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

void AppAld::grow_app()
{
  element = iarray[0];
  coord = iarray[1];
}

/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
------------------------------------------------------------------------- */

void AppAld::init_app()
{
  if (firsttime) {
    firsttime = 0;

    echeck = new int[nlocal];
    firstevent = (int *) memory->smalloc(nlocal*sizeof(int),"app:firstevent");
    //comneigh was defined to avoid double counting of common neighbor in site_propensity.
    // The number of distinct (k, dpropensity) pairs can exceed 12*maxneigh during the
    // H2O pulse when masked (coord-shifted) sites activate many type-II reactions,
    // which overflowed the old 12*maxneigh allocation. nlocal is a safe upper bound
    // (the dedup keys are site indices), with an explicit overflow guard below.
    comneigh = memory->grow(comneigh,nlocal,2,"app/ald:comneigh");
    // esites must hold every site up to 4th-neighbor distance of one site
    // (deep mask updates reach 4th neighbors); nlocal is a safe upper bound
    esites = (int *) memory->smalloc(nlocal*sizeof(int),"app:esites");
    // msites = lattice indices whose coord changed via put_mask/remove_mask
    // during update_coord; merged into the solver update list in site_event
    msites = (int *) memory->smalloc(nlocal*sizeof(int),"app:msites");
  }
  // site validity

  int flag = 0;
  for (int i = 0; i < nlocal; i++) {
    if (coord[i] < -1 || coord[i] > 8) flag = 1;
    if (element[i] < VACANCY || element[i] > Si) flag = 1;
  }
  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall) error->all(FLERR,"One or more sites have invalid values");
}
/* ---------------------------------------------------------------------- */

void AppAld::setup_app()
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

double AppAld::site_energy(int i)
{
  return 0.0;
}

/* ----------------------------------------------------------------------
   KMC method
   compute total propensity of owned site summed over possible events
------------------------------------------------------------------------- */

double AppAld::site_propensity(int i)
{
  int j,k,m;


  clear_events(i);

  double proball = 0.0;


  //type I, check species of sites and consider possible events

  // count_coordO was added here to prevent adsorption of HfX4 
  // on the low coordinate oxygen at sublayer

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
  	        if (nextneib >= nlocal)
  	          error->all(FLERR,"comneigh overflow in app_ald: nextneib exceeds nlocal");
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
	if (element[i] == vinput[m][0] && element[j] == vinput[m][1] && (coord[i] == vcoord[m] || vcoord[m] == 0) && (vpresson[m] == pressureOn || vpresson[m] == 0)) {
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

void AppAld::site_event(int i, class RandomPark *random)
{
  int j,k,m,n,mm,jj;

  // periodic audit: catch event-list / solver-propensity desyncs before the
  // solver hits them (seen as "Site has no events" / "Corrupt event list" /
  // "Illegal execution event" crashes on 400K/500K runs)
  if (++audit_counter >= AUDIT_INTERVAL) {
    audit_counter = 0;
    audit_consistency();
  }

  int elcoord = element[i];

  int isite = i2site[i];
  if (isite < 0 || firstevent[i] < 0 || propensity[isite] <= 0.0) {
    dump_crash_context("Site has no events in site_event",i,-1,-1,-1,-1);
    heal_desync(i);
    return;
  }

  double threshhold = random->uniform() * propensity[isite];
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
  update_coord(elcoord,i,j,which);

  // sequence of ALD, 
  // 1 is metal pulse, 3 purge, 2 oxygen pulse.
  if ((cycle+T1) > time ) {pressureOn = 1;}
  else if ((cycle+T1)<= time && time < (cycle+T1+T2)) {pressureOn = 3;}
  else if ((cycle+T1+T2) <= time && time < (cycle+T1+T2+T3)) {pressureOn = 2;}
  else if ((cycle+T1+T2+T3) <= time && time < (cycle+T1+T2+T3+T4)) {pressureOn = 3;}
  else {cycle += T1+T2+T3+T4; }

  

  int nsites = 0;
  isite = i2site[i];
  
 

  propensity[isite] = site_propensity(i);
  if (echeck[isite] == 0) {
    esites[nsites++] = isite;
    echeck[isite] = 1;
  }

  // go from site i to first and second neighbor in all type
  // for HfX4O, HfX4OH, and HfX2 go to the third and fourth neighbor to set mask
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
	//update mask up to fourth neighbor for AlX3O and AlX3OH
	if ((elcoord == O || elcoord == OH) && (element[i]== AlX3O || element[i]== AlX3OH || element[i]== HfX4O || element[i]== HfX4OH) && rstyle == 1) {

	   for (int ss = 0; ss< numneigh[mm]; ss++) {
	     int s = neighbor[mm][ss];
	     for (int ll = 0; ll< numneigh[s]; ll++) {
	       int l = neighbor[s][ll];
	       isite = i2site[l];
	       if (isite >= 0 && echeck[isite] == 0) {
			 propensity[isite] = site_propensity(l);
			 esites[nsites++] = isite;
			 echeck[isite] = 1;
	       }
	     }
	   }
	}
	if ((elcoord == AlX3OH || elcoord == AlX3O || elcoord == HfX4OH || elcoord == HfX4O) && (element[i]== OH  || element[i]== O) && rstyle == 1){

	   for (int ss = 0; ss< numneigh[mm]; ss++) {
	     int s = neighbor[mm][ss];
	     for (int ll = 0; ll< numneigh[s]; ll++) {
	       int l = neighbor[s][ll];
	       isite = i2site[l];
	       if (isite >= 0 && echeck[isite] == 0) {
			 propensity[isite] = site_propensity(l);
			 esites[nsites++] = isite;
			 echeck[isite] = 1;
	       }
	     }
	   }
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
        if ((elcoord == AlX2O || elcoord == AlHXO || elcoord == AlXO || elcoord == HfX2O || elcoord == HfHX2O || elcoord == HfH2X2O || elcoord == HfH4X4O) && element[i] == O ) {
	   for (int ss = 0; ss< numneigh[mm]; ss++) {
	     int s = neighbor[mm][ss];
	     isite = i2site[s];
	     if (isite >= 0 && echeck[isite] == 0) {
		 propensity[isite] = site_propensity(s);
		 esites[nsites++] = isite;
		 echeck[isite] = 1;
	     }
	   }
	}
        if ((elcoord == AlX2OH || elcoord == AlHXOH || elcoord == AlXOH || elcoord == HfX2OH || elcoord == HfHX2OH || elcoord == HfH2X2OH || elcoord == HfH4X4OH) && element[i] == OH ) {
	   for (int ss = 0; ss< numneigh[mm]; ss++) {
	     int s = neighbor[mm][ss];
	       isite = i2site[s];
	       if (isite >= 0 && echeck[isite] == 0) {
			 propensity[isite] = site_propensity(s);
			 esites[nsites++] = isite;
			 echeck[isite] = 1;
		    }
	     }
	   }
         if ((elcoord == AlX2 || elcoord == AlHX2 || elcoord == AlH2X2 || elcoord == HfX2 || elcoord == HfHX2 || elcoord == HfH2X2) && element[i] == VACANCY ) {
	   for (int ss = 0; ss< numneigh[mm]; ss++) {
	     int s = neighbor[mm][ss];
	     for (int ll = 0; ll< numneigh[s]; ll++) {
	       int l = neighbor[s][ll];
	       isite = i2site[l];
	       if (isite >= 0 && echeck[isite] == 0) {
			 propensity[isite] = site_propensity(l);
			 esites[nsites++] = isite;
			 echeck[isite] = 1;
	       }
	     }
	   }
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
   periodic consistency audit
   finds sites whose solver propensity is > 0 while their event list is
   empty (or vice versa) and repairs them, so a single corrupt site cannot
   abort the run with "Site has no events" / "Corrupt event list"
------------------------------------------------------------------------- */

void AppAld::audit_consistency()
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

void AppAld::dump_crash_context(const char *msg, int i,
                                int rstyle, int which, int jp, int kp)
{
  if (screen) {
    fprintf(screen,"\n=== AppAld event-list diagnostic: %s ===\n",msg);
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

int AppAld::heal_desync(int i)
{
  if (++heal_count > MAX_HEALS) {
    dump_crash_context("too many event-list repairs, aborting",i,-1,-1,-1,-1);
    error->all(FLERR,"Excessive event-list corruption in app_ald");
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

void AppAld::free_events_safe(int i)
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

/* ----------------------------------------------------------------------
   clear all events out of list for site I
   add cleared events to free list
------------------------------------------------------------------------- */

void AppAld::clear_events(int i)
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

void AppAld::add_event(int i, int rstyle, int which, double propensity,
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

void AppAld::grow_reactions(int rstyle)
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
   update c.n. for Hf and O, put and remove mask for relative sites
------------------------------------------------------------------------- */
void AppAld::update_coord(int elcoord, int i, int j, int which)
{
	if ((elcoord == O || elcoord == OH) && (element[i] == AlX3O || element[i] == AlX3OH)) {
		coord[i]++;
		put_mask(i);
	}
	else if ((elcoord == AlX3O || elcoord == AlX3OH) && (element[i] == O || element[i] == OH)){
		coord[i]--;
		remove_mask(i);
	}
	else if ((elcoord == AlX2 || elcoord == AlHX2 || elcoord == AlH2X2) && (element[i] == AlX || element[i] == AlHX || element[i] == Al)){
		remove_mask(i);
		coord[i]--;
		if (element[i] == Al) coord[i]=coord[i]-1;
	}
	else if ((elcoord == AlX || elcoord == AlHX) && (element[i] == AlHX2 || element[i] == AlX2 || element[i] == AlH2X2)){
		coord[i]++;
		put_mask(i);
	}

	else if (j >= 0 && elcoord == AlX2O && element[i] == O && element[j]==AlX2) {
		remove_mask(i);
		count_coord(i,j);
		put_mask(j);
	}
	else if (j >= 0 && elcoord == AlHXO && element[i] == O && element[j]==AlHX) {
		remove_mask(i);
		count_coord(i,j);
		put_mask(j);
	}
	else if (j >= 0 && elcoord == AlX2OH && element[i] == OH && element[j]==AlX2) {
		remove_mask(i);
		count_coord(i,j);
		put_mask(j);
	}
	else if (j >= 0 && elcoord == AlHXOH && element[i] == OH && element[j]==AlHX) {
		remove_mask(i);
		count_coord(i,j);
		put_mask(j);
	}
	else if (j >= 0 && elcoord == AlXO && element[i] == O && element[j]==AlX) {
		remove_mask(i);
		count_coord(i,j);
		put_mask(j);
	}
	else if (j >= 0 && elcoord == AlXOH && element[i] == OH && element[j]==AlX) {
		remove_mask(i);
		count_coord(i,j);
		put_mask(j);
	}
	else if (j >= 0 && (AlX2 <= elcoord && elcoord <= AlH2X2) && element[i] == VACANCY && (AlX2O <= element[j] && element[j] <= AlXOH)) {
		remove_mask(i);
		count_coord(i,j);
		put_mask(j);
	}
	else if (elcoord == AlHX && element[i] == Al){
		coord[i]--;
	}
	else if (elcoord == Al && element[i] == AlHX){
		coord[i]++;
	}
	else if ((elcoord == OH2AlHX || elcoord == OH2AlX) && element[i] == OH2Al){
		coord[i]--;
	}

		// HfO2 branches
	else if ((elcoord == O || elcoord == OH) && (element[i] == HfX4O || element[i] == HfX4OH)) {
		coord[i]++;
		put_mask(i);
	}
	else if ((elcoord == HfX4O || elcoord == HfX4OH) && (element[i] == O || element[i] == OH)){
		coord[i]--;
		remove_mask(i);
	}
	else if ((elcoord == HfX2 || elcoord == HfHX2 || elcoord == HfH2X2) && (element[i] == HfX || element[i] == HfHX || element[i] == Hf)){
		remove_mask(i);
		coord[i]--;
		if (element[i] == Hf) coord[i]=coord[i]-1;
	}
	else if ((elcoord == HfX || elcoord == HfHX) && (element[i] == HfHX2 || element[i] == HfX2 || element[i] == HfH2X2)){
		coord[i]++;
		put_mask(i);
	}

	else if (j >= 0 && elcoord == HfX2O && element[i] == O && element[j]==HfX2) {
		remove_mask(i);
		count_coord(i,j);
		put_mask(j);
	}
	else if (j >= 0 && elcoord == HfHX2O && element[i] == O && element[j]==HfHX2) {
		remove_mask(i);
		count_coord(i,j);
		put_mask(j);
	}
	else if (j >= 0 && elcoord == HfX2OH && element[i] == OH && element[j]==HfX2) {
		remove_mask(i);
		count_coord(i,j);
		put_mask(j);
	}
	else if (j >= 0 && elcoord == HfHX2OH && element[i] == OH && element[j]==HfHX2) {
		remove_mask(i);
		count_coord(i,j);
		put_mask(j);
	}
	else if (j >= 0 && elcoord == HfH2X2O && element[i] == O && element[j]==HfH2X2) {
		remove_mask(i);
		count_coord(i,j);
		put_mask(j);
	}
	else if (j >= 0 && elcoord == HfH2X2OH && element[i] == OH && element[j]==HfH2X2) {
		remove_mask(i);
		count_coord(i,j);
		put_mask(j);
	}
	else if (j >= 0 && elcoord == HfH4X4O && element[i] == O && element[j]==HfH2X2) {
		remove_mask(i);
		count_coord(i,j);
		put_mask(j);
	}
	else if (j >= 0 && elcoord == HfH4X4OH && element[i] == OH && element[j]==HfH2X2) {
		remove_mask(i);
		count_coord(i,j);
		put_mask(j);
	}
	else if (j >= 0 && (HfX2 <= elcoord && elcoord <= HfH2X2) && element[i] == VACANCY && (HfX2O <= element[j] && element[j] <= HfH2X2OH)) {
		remove_mask(i);
		count_coord(i,j);
		put_mask(j);
	}
	else if (elcoord == HfHX && element[i] == Hf){
		coord[i]--;
	}
	else if (elcoord == Hf && element[i] == HfHX){
		coord[i]++;
	}
	else if ((elcoord == OH2HfHX || elcoord == OH2HfX) && element[i] == OH2Hf){
		coord[i]--;
	}

	// densification of water molecule (Hf)
	else if (j >= 0 && (OH2HfX <= elcoord && elcoord <= OHHfHX) && (HfHX <= element[i] && element[i] <= Hf) && element[j]== OH2) {
		if((elcoord == OH2HfX || elcoord == OH2HfHX) && element[i] == Hf && element[j] == OH2) 
			{coord[i]--; }
		count_coord(i,j);
	}

	// the reverse of densification of water (Hf)
	else if (elcoord==OH2 && element[i]==VACANCY && (OH2HfX <= element[j] && element[j] <= OH2Hf) && j != -1) {
		count_coord(i,j);
	}

// densifiacation of water molecule
	else if (j >= 0 && (OH2AlX <= elcoord && elcoord <= OH2Al) && (AlHX <= element[i] && element[i] <= Al) && element[j]== OH2) {
		if((elcoord == OH2AlX || elcoord == OH2AlHX) && element[i] == Al && element[j] == OH2) 
			{coord[i]--; }
		count_coord(i,j);
	}

	// the reverse of densification of water
	else if (elcoord==OH2 && element[i]==VACANCY && (OH2AlX <= element[j] && element[j] <= OH2Al) && j != -1) {
		count_coord(i,j);
	}



        else if ((element[i]==OH || element[i]==O) && coord[i]==1 && pressureOn == 1) {
                count_coordO(i);
        }
}

/* ----------------------------------------------------------------------
   put mask for affected sites
------------------------------------------------------------------------- */
void AppAld::put_mask(int i)
{ 
        int isite = i2site[i];
        if (isite < 0) return;
	int nsites = 0;
	esites[nsites++] = isite;
        echeck[isite] = 1;
	if (element[i] == AlX3O || element[i] == AlX3OH || element[i] == HfX4O || element[i] == HfX4OH ){
	  	for (int n = 0; n < numneigh[i]; n++) {
			int nn = neighbor[i][n];
			for (int k = 0; k < numneigh[nn]; k++){
				int kk = neighbor[nn][k];
				isite = i2site[kk];
				if (isite >= 0 && echeck[isite] == 0) {
					coord[isite]=coord[isite]-10;
					esites[nsites++] = isite;
					echeck[isite] = 1;
					track_msite(kk);
				}
				for (int m = 0; m < numneigh[kk]; m++) {
					int mm = neighbor[kk][m];
					for (int s = 0; s < numneigh[mm]; s++) {
						int ss = neighbor[mm][s];
						isite = i2site[ss];
						if (isite >= 0 && echeck[isite] == 0) {
							coord[isite] = coord[isite]-10;
							esites[nsites++] = isite;
							echeck[isite] = 1;
							track_msite(ss);
						}
					}
				}
			}
	        }
	}
	else if (element[i] == AlX2 || element[i] == AlHX2 || element[i] == AlH2X2 || element[i] == HfX2 || element[i] == HfHX2 || element[i] == HfH2X2){
	  	for (int n = 0; n < numneigh[i]; n++) {
			int nn = neighbor[i][n];
			isite = i2site[nn];
			if (isite >= 0 && echeck[isite] == 0) {
				coord[isite]=coord[isite]-10;
				esites[nsites++] = isite;
				echeck[isite] = 1;
				track_msite(nn);
			}
			for (int k = 0; k < numneigh[nn]; k++){
				int kk = neighbor[nn][k];
				isite = i2site[kk];
				if (isite >= 0 && echeck[isite] == 0) {
					esites[nsites++] = isite;
					echeck[isite] = 1;
				}
				for (int m = 0; m < numneigh[kk]; m++) {
					int mm = neighbor[kk][m];
					isite = i2site[mm];
					if (isite >= 0 && echeck[isite] == 0) {
						coord[isite] = coord[isite]-10;
						esites[nsites++] = isite;
						echeck[isite] = 1;
						track_msite(mm);
					}
				}
			}
		}
	}
	else if ( (AlX2O <= element[i] && element[i] <= AlXOH) || (HfX2O <= element[i] && element[i] <= HfH2X2OH) ){
	  	for (int n = 0; n < numneigh[i]; n++) {
			int nn = neighbor[i][n];
			for (int k = 0; k < numneigh[nn]; k++){
				int kk = neighbor[nn][k];
				isite = i2site[kk];
				if (isite >= 0 && echeck[isite] == 0) {
					coord[isite]=coord[isite]-10;
					esites[nsites++] = isite;
					echeck[isite] = 1;
					track_msite(kk);
				}
				for (int m = 0; m < numneigh[kk]; m++) {
					int mm = neighbor[kk][m];
					for (int s = 0; s < numneigh[mm]; s++) {
						int ss = neighbor[mm][s];
						isite = i2site[ss];
						if (isite >= 0 && echeck[isite] == 0) {
							coord[isite] = coord[isite]-10;
							esites[nsites++] = isite;
							echeck[isite] = 1;
							track_msite(ss);
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
void AppAld::remove_mask(int i)
{ 
        int isite = i2site[i];
        if (isite < 0) return;
	int nsites = 0;
	esites[nsites++] = isite;
        echeck[isite] = 1;
	if (element[i] == AlX3O || element[i] == AlX3OH || element[i] == HfX4O || element[i] == HfX4OH ){
	  	for (int n = 0; n < numneigh[i]; n++) {
			int nn = neighbor[i][n];
			for (int k = 0; k < numneigh[nn]; k++){
				int kk = neighbor[nn][k];
				isite = i2site[kk];
				if (isite >= 0 && echeck[isite] == 0) {
					coord[isite]=coord[isite]+10;
					esites[nsites++] = isite;
					echeck[isite] = 1;
					track_msite(kk);
				}
				for (int m = 0; m < numneigh[kk]; m++) {
					int mm = neighbor[kk][m];
					for (int s = 0; s < numneigh[mm]; s++) {
						int ss = neighbor[mm][s];
						isite = i2site[ss];
						if (isite >= 0 && echeck[isite] == 0) {
							coord[isite] = coord[isite]+10;
							esites[nsites++] = isite;
							echeck[isite] = 1;
							track_msite(ss);
						}
					}
				}
			}
	        }
	}
	else if (element[i] == AlX2 || element[i] == AlHX2 || element[i] == AlH2X2 || element[i] == HfX2 || element[i] == HfHX2 || element[i] == HfH2X2){
	  	for (int n = 0; n < numneigh[i]; n++) {
			int nn = neighbor[i][n];
			isite = i2site[nn];
			if (isite >= 0 && echeck[isite] == 0) {
				coord[isite]=coord[isite]+10;
				esites[nsites++] = isite;
				echeck[isite] = 1;
				track_msite(nn);
			}
			for (int k = 0; k < numneigh[nn]; k++){
				int kk = neighbor[nn][k];
				isite = i2site[kk];
				if (isite >= 0 && echeck[isite] == 0) {
					esites[nsites++] = isite;
					echeck[isite] = 1;
				}
				for (int m = 0; m < numneigh[kk]; m++) {
					int mm = neighbor[kk][m];
					isite = i2site[mm];
					if (isite >= 0 && echeck[isite] == 0) {
						coord[isite] = coord[isite]+10;
						esites[nsites++] = isite;
						echeck[isite] = 1;
						track_msite(mm);
					}
				}
			}
		}
	}
	else if ( (AlX2O <= element[i] && element[i] <= AlXOH) || (HfX2O <= element[i] && element[i] <= HfH2X2OH) ){
	  	for (int n = 0; n < numneigh[i]; n++) {
			int nn = neighbor[i][n];
			for (int k = 0; k < numneigh[nn]; k++){
				int kk = neighbor[nn][k];
				isite = i2site[kk];
				if (isite >= 0 && echeck[isite] == 0) {
					coord[isite]=coord[isite]+10;
					esites[nsites++] = isite;
					echeck[isite] = 1;
					track_msite(kk);
				}
				for (int m = 0; m < numneigh[kk]; m++) {
					int mm = neighbor[kk][m];
					for (int s = 0; s < numneigh[mm]; s++) {
						int ss = neighbor[mm][s];
						isite = i2site[ss];
						if (isite >= 0 && echeck[isite] == 0) {
							coord[isite] = coord[isite]+10;
							esites[nsites++] = isite;
							echeck[isite] = 1;
							track_msite(ss);
						}
					}
				}
			}
	        }
	}
	else if (element[i] == VACANCY || element[i] == AlX || element[i] == AlHX || element[i] == Al || element[i] == HfX || element[i] == HfHX || element[i] == Hf){
	  	for (int n = 0; n < numneigh[i]; n++) {
			int nn = neighbor[i][n];
			isite = i2site[nn];
			if (isite >= 0 && echeck[isite] == 0) {
				coord[isite]=coord[isite]+10;
				esites[nsites++] = isite;
				echeck[isite] = 1;
				track_msite(nn);
			}
			for (int k = 0; k < numneigh[nn]; k++){
				int kk = neighbor[nn][k];
				isite = i2site[kk];
				if (isite >= 0 && echeck[isite] == 0) {
					esites[nsites++] = isite;
					echeck[isite] = 1;
				}
				for (int m = 0; m < numneigh[kk]; m++) {
					int mm = neighbor[kk][m];
					isite = i2site[mm];
					if (isite >= 0 && echeck[isite] == 0) {
						coord[isite] = coord[isite]+10;
						esites[nsites++] = isite;
						echeck[isite] = 1;
						track_msite(mm);
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
void AppAld::count_coord(int i, int j)
{
    if (j < 0 || i < 0) return;
    if ((element[i] == O || element[i] == OH) &&  ((AlX2 <= element[j] && element[j] <= AlH2X2) || (HfX2 <= element[j] && element[j] <= HfH2X2)) ){
	coord[j]=coord[j]+2;
	for (int s = 0; s < numneigh[j]; s++){
		int nn = neighbor[j][s];
		if (element[nn] == O || element[nn] == OH || element[nn] == OH2) {
			coord[j]++;
			if (i != nn) coord[nn]++;
		}
	}
    }
    else if (((AlHX <= element[i] && element[i] <= Al) || (HfHX <= element[i] && element[i] <= Hf)) && element[j] == OH2) {
	for (int s = 0; s < numneigh[j]; s++){
		int nn = neighbor[j][s];
		if ( (AlX2 <= element[nn] && element[nn] <= OH2Al) || (HfX2 <= element[nn] && element[nn] <= OHHfHX) ) {
			coord[j]++;
			coord[nn]++;
		}
	}
    }
    else if (((OH2AlX <= element[j] && element[j] <= OH2Al) || (OH2HfX <= element[j] && element[j] <= OH2Hf)) && element[i] == VACANCY) {
	for (int s = 0; s < numneigh[i]; s++){
		int nn = neighbor[i][s];
		if ( (AlX2 <= element[nn] && element[nn] <= OH2Al) || (HfX2 <= element[nn] && element[nn] <= OHHfHX) ) {
			coord[i]--;
			coord[nn]--;
		}
	}
    }
    else if (element[i] == VACANCY &&  ((AlX2O <= element[j] && element[j] <= AlXOH) || (HfX2O <= element[j] && element[j] <= HfH2X2OH))){
	coord[i]=coord[i]-3;
	for (int s = 0; s < numneigh[i]; s++){
		int nn = neighbor[i][s];
		if (element[nn] == O || element[nn] == OH || element[nn] == OH2) {
			coord[i]--;
			if (j != nn) coord[nn]--;
		}
	}
    }
}

/* ----------------------------------------------------------------------
   remember a lattice site whose coord was changed by put_mask/remove_mask,
   so site_event can refresh its solver propensity (dedup by linear scan)
------------------------------------------------------------------------- */
void AppAld::track_msite(int ilat)
{
  for (int q = 0; q < nmsites; q++)
    if (msites[q] == ilat) return;
  if (nmsites >= nlocal)
    error->all(FLERR,"msites overflow in app_ald");
  msites[nmsites++] = ilat;
}

/* ----------------------------------------------------------------------
   count c.n of oxygen before adsorption
------------------------------------------------------------------------- */
void AppAld::count_coordO(int i)
{
    int fullO = 0;
    int emptyO = 0;
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
			  if (element[ss] == O || element[ss] == OH || element[ss] == OH2) {fullO++;}
			  else if (element[ss] == VACANCY) {emptyO++;}
                          else {}
		          esites[nsites++] = isite;
		          echeck[isite] = 1;
		        }
    
		}
	}
   totalS = fullO+emptyO;
   if ( float(fullO) > 4*totalS/5 ) {coord[i]=2; }
   for (int m = 0; m < nsites; m++)  {echeck[esites[m]] = 0; esites[m]=0;}
}
