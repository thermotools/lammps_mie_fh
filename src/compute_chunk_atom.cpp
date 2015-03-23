/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "compute_chunk_atom.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "domain.h"
#include "region.h"
#include "lattice.h"
#include "modify.h"
#include "fix_store.h"
#include "comm.h"
#include "group.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"

#include <map>

using namespace LAMMPS_NS;

enum{BIN1D,BIN2D,BIN3D,TYPE,MOLECULE,COMPUTE,FIX,VARIABLE};
enum{LOWER,CENTER,UPPER,COORD};
enum{BOX,LATTICE,REDUCED};
enum{NODISCARD,MIXED,YESDISCARD};
enum{ONCE,NFREQ,EVERY};              // used in several files
enum{LIMITMAX,LIMITEXACT};

#define IDMAX 1024*1024
#define INVOKED_PERATOM 8

// allocate space for static class variable

ComputeChunkAtom *ComputeChunkAtom::cptr;

/* ---------------------------------------------------------------------- */

ComputeChunkAtom::ComputeChunkAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal compute chunk/atom command");

  peratom_flag = 1;
  size_peratom_cols = 0;
  create_attribute = 1;

  // chunk style and its args

  int iarg;

  binflag = 0;
  ncoord = 0;
  cfvid = NULL;

  if (strcmp(arg[3],"bin/1d") == 0) {
    binflag = 1;
    which = BIN1D;
    ncoord = 1;
    iarg = 4;
    readdim(narg,arg,iarg,0);
    iarg += 3;
  } else if (strcmp(arg[3],"bin/2d") == 0) {
    binflag = 1;
    which = BIN2D;
    ncoord = 2;
    iarg = 4;
    readdim(narg,arg,iarg,0);
    readdim(narg,arg,iarg+3,1);
    iarg += 6;
  } else if (strcmp(arg[3],"bin/3d") == 0) {
    binflag = 1;
    which = BIN3D;
    ncoord = 3;
    iarg = 4;
    readdim(narg,arg,iarg,0);
    readdim(narg,arg,iarg+3,1);
    readdim(narg,arg,iarg+6,2);
    iarg += 9;

  } else if (strcmp(arg[3],"type") == 0) {
    which = TYPE;
    iarg = 4;
  } else if (strcmp(arg[3],"molecule") == 0) {
    which = MOLECULE;
    iarg = 4;

  } else if (strstr(arg[3],"c_") == arg[3] || 
             strstr(arg[3],"f_") == arg[3] || 
             strstr(arg[3],"v_") == arg[3]) {
    if (arg[3][0] == 'c') which = COMPUTE;
    else if (arg[3][0] == 'f') which = FIX;
    else if (arg[3][0] == 'v') which = VARIABLE;
    iarg = 4;

    int n = strlen(arg[3]);
    char *suffix = new char[n];
    strcpy(suffix,&arg[3][2]);

    char *ptr = strchr(suffix,'[');
    if (ptr) {
      if (suffix[strlen(suffix)-1] != ']')
        error->all(FLERR,"Illegal fix ave/atom command");
      argindex = atoi(ptr+1);
      *ptr = '\0';
    } else argindex = 0;
    
    n = strlen(suffix) + 1;
    cfvid = new char[n];
    strcpy(cfvid,suffix);
    delete [] suffix;

  } else error->all(FLERR,"Illegal compute chunk/atom command");

  // optional args

  regionflag = 0;
  idregion = NULL;
  nchunksetflag = 0;
  nchunkflag = EVERY;
  limit = 0;
  limitstyle = LIMITMAX;
  limitfirst = 0;
  idsflag = EVERY;
  compress = 0;
  int discardsetflag = 0;
  discard = MIXED;
  minflag[0] = LOWER;
  minflag[1] = LOWER;
  minflag[2] = LOWER;
  maxflag[0] = UPPER;
  maxflag[1] = UPPER;
  maxflag[2] = UPPER;
  scaleflag = LATTICE;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute chunk/atom command");
      int iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1)
        error->all(FLERR,"Region ID for compute chunk/atom does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      idregion = new char[n];
      strcpy(idregion,arg[iarg+1]);
      regionflag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"nchunk") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute chunk/atom command");
      if (strcmp(arg[iarg+1],"once") == 0) nchunkflag = ONCE;
      else if (strcmp(arg[iarg+1],"every") == 0) nchunkflag = EVERY;
      else error->all(FLERR,"Illegal compute chunk/atom command");
      nchunksetflag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"limit") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute chunk/atom command");
      limit = force->inumeric(FLERR,arg[iarg+1]);
      if (limit < 0) error->all(FLERR,"Illegal compute chunk/atom command");
      if (limit && !compress) limitfirst = 1;
      iarg += 2;
      if (limit) {
        if (iarg+1 > narg) 
          error->all(FLERR,"Illegal compute chunk/atom command");
        if (strcmp(arg[iarg+1],"max") == 0) limitstyle = LIMITMAX;
        else if (strcmp(arg[iarg+1],"exact") == 0) limitstyle = LIMITEXACT;
        else error->all(FLERR,"Illegal compute chunk/atom command");
        iarg++;
      }
    } else if (strcmp(arg[iarg],"ids") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute chunk/atom command");
      if (strcmp(arg[iarg+1],"once") == 0) idsflag = ONCE;
      else if (strcmp(arg[iarg+1],"nfreq") == 0) idsflag = NFREQ;
      else if (strcmp(arg[iarg+1],"every") == 0) idsflag = EVERY;
      else error->all(FLERR,"Illegal compute chunk/atom command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"compress") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute chunk/atom command");
      else if (strcmp(arg[iarg+1],"no") == 0) compress = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) compress = 1;
      else error->all(FLERR,"Illegal compute chunk/atom command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"discard") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute chunk/atom command");
      if (strcmp(arg[iarg+1],"mixed") == 0) discard = MIXED;
      else if (strcmp(arg[iarg+1],"no") == 0) discard = NODISCARD;
      else if (strcmp(arg[iarg+1],"yes") == 0) discard = YESDISCARD;
      else error->all(FLERR,"Illegal compute chunk/atom command");
      discardsetflag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"bound") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal compute chunk/atom command");
      int idim;
      if (strcmp(arg[iarg+1],"x") == 0) idim = 0;
      else if (strcmp(arg[iarg+1],"y") == 0) idim = 1;
      else if (strcmp(arg[iarg+1],"z") == 0) idim = 2;
      else error->all(FLERR,"Illegal compute chunk/atom command");
      if (strcmp(arg[iarg+2],"lower") == 0) minflag[idim] = LOWER;
      else minflag[idim] = COORD;
      if (minflag[idim] == COORD) 
        minvalue[idim] = force->numeric(FLERR,arg[iarg+2]);
      if (strcmp(arg[iarg+3],"upper") == 0) maxflag[idim] = UPPER;
      else maxflag[idim] = COORD;
      if (maxflag[idim] == COORD) 
        maxvalue[idim] = force->numeric(FLERR,arg[iarg+3]);
      else error->all(FLERR,"Illegal compute chunk/atom command");
      iarg += 4;
    } else if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute chunk/atom command");
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = BOX;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = LATTICE;
      else if (strcmp(arg[iarg+1],"reduced") == 0) scaleflag = REDUCED;
      else error->all(FLERR,"Illegal compute chunk/atom command");
      iarg += 2;
    } else error->all(FLERR,"Illegal compute chunk/atom command");
  }

  // set nchunkflag and discard to default values if not explicitly set
  // for binning style, also check in init() if simulation box is static,
  //   which sets nchunkflag = ONCE

  if (!nchunksetflag) {
    if (binflag) {
      if (scaleflag == REDUCED) nchunkflag = ONCE;
      else nchunkflag = EVERY;
    }
    if (which == TYPE) nchunkflag = ONCE;
    if (which == MOLECULE) {
      if (regionflag) nchunkflag = EVERY;
      else nchunkflag = ONCE;
    }
    if (compress) nchunkflag = EVERY;
  }

  if (!discardsetflag) {
    if (binflag) discard = MIXED;
    else discard = YESDISCARD;
  }

  // error checks

  if (which == MOLECULE && !atom->molecule_flag) 
    error->all(FLERR,"Compute chunk/atom molecule for non-molecular system");

  if (!binflag && discard == MIXED)
    error->all(FLERR,"Compute chunk/atom without bins "
               "cannot use discard mixed");
  if (which == BIN1D && delta[0] <= 0.0) 
    error->all(FLERR,"Illegal compute chunk/atom command");
  if (which == BIN2D && (delta[0] <= 0.0 || delta[1] <= 0.0))
    error->all(FLERR,"Illegal compute chunk/atom command");
  if (which == BIN3D && 
      (delta[0] <= 0.0 || delta[1] <= 0.0 || delta[2] <= 0.0))
      error->all(FLERR,"Illegal compute chunk/atom command");

  if (which == COMPUTE) {
    int icompute = modify->find_compute(cfvid);
    if (icompute < 0)
      error->all(FLERR,"Compute ID for compute chunk/atom does not exist");
    if (modify->compute[icompute]->peratom_flag == 0)
      error->all(FLERR,
                 "Compute chunk/atom compute does not calculate "
                 "per-atom values");
    if (argindex == 0 &&
        modify->compute[icompute]->size_peratom_cols != 0)
      error->all(FLERR,"Compute chunk/atom compute does not "
                 "calculate a per-atom vector");
    if (argindex && modify->compute[icompute]->size_peratom_cols == 0)
      error->all(FLERR,"Compute chunk/atom compute does not "
                 "calculate a per-atom array");
    if (argindex &&
        argindex > modify->compute[icompute]->size_peratom_cols)
      error->all(FLERR,"Compute chunk/atom compute array is "
                 "accessed out-of-range");
  }

  if (which == FIX) {
    int ifix = modify->find_fix(cfvid);
    if (ifix < 0)
      error->all(FLERR,"Fix ID for compute chunk/atom does not exist");
    if (modify->fix[ifix]->peratom_flag == 0)
      error->all(FLERR,"Compute chunk/atom fix does not calculate "
                 "per-atom values");
    if (argindex == 0 && modify->fix[ifix]->size_peratom_cols != 0)
      error->all(FLERR,
                 "Compute chunk/atom fix does not calculate a per-atom vector");
    if (argindex && modify->fix[ifix]->size_peratom_cols == 0)
      error->all(FLERR,
                 "Compute chunk/atom fix does not calculate a per-atom array");
    if (argindex && argindex > modify->fix[ifix]->size_peratom_cols)
      error->all(FLERR,"Compute chunk/atom fix array is accessed out-of-range");
  }

  if (which == VARIABLE) {
    int ivariable = input->variable->find(cfvid);
    if (ivariable < 0)
      error->all(FLERR,"Variable name for compute chunk/atom does not exist");
    if (input->variable->atomstyle(ivariable) == 0)
      error->all(FLERR,"Compute chunk/atom variable is not "
                 "atom-style variable");
  }

  // setup scaling

  if (binflag) {
    if (domain->triclinic == 1 && scaleflag != REDUCED)
      error->all(FLERR,"Compute chunk/atom for triclinic boxes "
                 "requires units reduced");
  }

  if (scaleflag == LATTICE) {
    xscale = domain->lattice->xlattice;
    yscale = domain->lattice->ylattice;
    zscale = domain->lattice->zlattice;
  } else xscale = yscale = zscale = 1.0;

  // apply scaling factors

  if (binflag) {
    double scale;
    if (which == BIN1D) ndim = 1;
    if (which == BIN2D) ndim = 2;
    if (which == BIN3D) ndim = 3;
    for (int idim = 0; idim < ndim; idim++) {
      if (dim[idim] == 0) scale = xscale;
      else if (dim[idim] == 1) scale = yscale;
      else if (dim[idim] == 2) scale = zscale;
      delta[idim] *= scale;
      invdelta[idim] = 1.0/delta[idim];
      if (originflag[idim] == COORD) origin[idim] *= scale;
      if (minflag[idim] == COORD) minvalue[idim] *= scale;
      if (maxflag[idim] == COORD) maxvalue[idim] *= scale;
    }
  }

  // initialize chunk vector and per-chunk info

  nmax = 0;
  chunk = NULL;
  nmaxint = 0;
  ichunk = NULL;
  exclude = NULL;

  nchunk = 0;
  chunk_volume_scalar = 1.0;
  chunk_volume_vec = NULL;
  coord = NULL;
  chunkID = NULL;

  // computeflag = 1 if this compute might invoke another compute
  // during assign_chunk_ids()

  if (which == COMPUTE || which == FIX || which == VARIABLE) computeflag = 1;
  else computeflag = 0;

  // other initializations

  invoked_setup = -1;
  invoked_ichunk = -1;

  id_fix = NULL;
  fixstore = NULL;

  if (compress) hash = new std::map<tagint,int>();
  else hash = NULL;

  maxvar = 0;
  varatom = NULL;

  lockcount = 0;
  lockfix = NULL;

  if (which == MOLECULE) molcheck = 1;
  else molcheck = 0;
}

/* ---------------------------------------------------------------------- */

ComputeChunkAtom::~ComputeChunkAtom()
{
  // check nfix in case all fixes have already been deleted

  if (modify->nfix) modify->delete_fix(id_fix);
  delete [] id_fix;

  memory->destroy(chunk);
  memory->destroy(ichunk);
  memory->destroy(exclude);
  memory->destroy(chunk_volume_vec);
  memory->destroy(coord);
  memory->destroy(chunkID);

  delete [] idregion;
  delete [] cfvid;
  delete hash;

  memory->destroy(varatom);
}

/* ---------------------------------------------------------------------- */

void ComputeChunkAtom::init()
{
  // set and check validity of region

  if (regionflag) {
    int iregion = domain->find_region(idregion);
    if (iregion == -1)
      error->all(FLERR,"Region ID for compute chunk/atom does not exist");
    region = domain->regions[iregion];
  }

  // set compute,fix,variable

  if (which == COMPUTE) {
    int icompute = modify->find_compute(cfvid);
    if (icompute < 0)
      error->all(FLERR,"Compute ID for compute chunk/atom does not exist");
    cchunk = modify->compute[icompute];
  } else if (which == FIX) {
    int ifix = modify->find_fix(cfvid);
    if (ifix < 0)
      error->all(FLERR,"Fix ID for compute chunk/atom does not exist");
    fchunk = modify->fix[ifix];
  } else if (which == VARIABLE) {
    int ivariable = input->variable->find(cfvid);
    if (ivariable < 0)
      error->all(FLERR,"Variable name for compute chunk/atom does not exist");
    vchunk = ivariable;
  }

  // for style MOLECULE, check that no mol IDs exceed MAXSMALLINT
  // don't worry about group or optional region

  if (which == MOLECULE) {
    tagint *molecule = atom->molecule;
    int nlocal = atom->nlocal;
    tagint maxone = -1;
    for (int i = 0; i < nlocal; i++)
      if (molecule[i] > maxone) maxone = molecule[i];
    tagint maxall;
    MPI_Allreduce(&maxone,&maxall,1,MPI_LMP_TAGINT,MPI_MAX,world);
    if (maxall > MAXSMALLINT) 
      error->all(FLERR,"Molecule IDs too large for compute chunk/atom");
  }

  // for binning, if nchunkflag not already set, set it to ONCE or EVERY
  // depends on whether simulation box size is static or dynamic
  // reset invoked_setup if this is not first run and box just became static

  if (binflag && !nchunksetflag && !compress && scaleflag != REDUCED) {
    if (domain->box_change_size == 0) {
      if (nchunkflag == EVERY && invoked_setup >= 0) invoked_setup = -1;
      nchunkflag = ONCE;
    } else nchunkflag = EVERY;
  }

  // require nchunkflag = ONCE if idsflag = ONCE
  // b/c nchunk cannot change if chunk IDs are frozen
  // can't check until now since nchunkflag may have been adjusted in init()

  if (idsflag == ONCE && nchunkflag != ONCE)
    error->all(FLERR,"Compute chunk/atom ids once but nchunk is not once");

  // create/destroy fix STORE for persistent chunk IDs as needed
  // need to wait until init() so that fix ave/chunk command(s) are in place
  // they increment lockcount if they lock this compute
  // fixstore ID = compute-ID + COMPUTE_STORE, fix group = compute group
  // fixstore initializes all values to 0.0

  if (lockcount && !fixstore) {
    int n = strlen(id) + strlen("_COMPUTE_STORE") + 1;
    id_fix = new char[n];
    strcpy(id_fix,id);
    strcat(id_fix,"_COMPUTE_STORE");
  
    char **newarg = new char*[5];
    newarg[0] = id_fix;
    newarg[1] = group->names[igroup];
    newarg[2] = (char *) "STORE";
    newarg[3] = (char *) "1";
    newarg[4] = (char *) "1";
    modify->add_fix(5,newarg);
    fixstore = (FixStore *) modify->fix[modify->nfix-1];
    delete [] newarg;
  }

  if (!lockcount && fixstore) {
    delete fixstore;
    fixstore = NULL;
  }
}

/* ----------------------------------------------------------------------
   invoke setup_chunks and/or compute_ichunk if only done ONCE
   so that nchunks or chunk IDs are assigned when this compute was specified
     as opposed to first times compute_peratom() or compute_ichunk() is called
------------------------------------------------------------------------- */

void ComputeChunkAtom::setup()
{
  if (nchunkflag == ONCE) setup_chunks();
  if (idsflag == ONCE) compute_ichunk();
}

/* ----------------------------------------------------------------------
   only called by classes that use per-atom computes in standard way
     dump, variable, thermo output, other computes, etc
   not called by fix ave/chunk or compute chunk commands
     they invoke setup_chunks() and compute_ichunk() directly
------------------------------------------------------------------------- */

void ComputeChunkAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow floating point chunk vector if necessary

  if (atom->nlocal > nmax) {
    memory->destroy(chunk);
    nmax = atom->nmax;
    memory->create(chunk,nmax,"chunk/atom:chunk");
    vector_atom = chunk;
  }
  
  setup_chunks();
  compute_ichunk();

  // copy integer indices into floating-point chunk vector

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) chunk[i] = ichunk[i];
}

/* ----------------------------------------------------------------------
   set lock, so that nchunk will not change from startstep to stopstep
   called by fix ave/chunk for duration of its Nfreq epoch
   OK if called by multiple fix ave/chunk commands
     error if all callers do not have same duration
     last caller holds the lock, so it can also unlock
   lockstop can be positive for final step of finite-size time window
   or can be -1 for infinite-size time window
------------------------------------------------------------------------- */

void ComputeChunkAtom::lock(Fix *fixptr, bigint startstep, bigint stopstep)
{
  if (lockfix == NULL) {
    lockfix = fixptr;
    lockstart = startstep;
    lockstop = stopstep;
    return;
  }

  if (startstep != lockstart || stopstep != lockstop)
    error->all(FLERR,"Two fix ave commands using "
               "same compute chunk/atom command in incompatible ways");

  // set lock to last calling Fix, since it will be last to unlock()

  lockfix = fixptr;
}

/* ----------------------------------------------------------------------
   unset lock
   can only be done by fix ave/chunk command that holds the lock
------------------------------------------------------------------------- */

void ComputeChunkAtom::unlock(Fix *fixptr)
{
  if (fixptr != lockfix) return;
  lockfix = NULL;
}

/* ----------------------------------------------------------------------
   assign chunk IDs from 1 to Nchunk to every atom, or 0 if not in chunk
------------------------------------------------------------------------- */

void ComputeChunkAtom::compute_ichunk()
{
  int i;

  // skip if already done on this step

  if (invoked_ichunk == update->ntimestep) return;

  // if old IDs persist via storage in fixstore, then just retrieve them
  // yes if idsflag = ONCE, and already done once
  //   or if idsflag = NFREQ and lock is in place and are on later timestep
  // else proceed to recalculate per-atom chunk assignments

  int restore = 0;
  if (idsflag == ONCE && invoked_ichunk >= 0) restore = 1;
  if (idsflag == NFREQ && lockfix && update->ntimestep > lockstart) restore = 1;

  if (restore) {
    invoked_ichunk = update->ntimestep;
    double *vstore = fixstore->vstore;
    int nlocal = atom->nlocal;
    for (i = 0; i < nlocal; i++) ichunk[i] = static_cast<int> (vstore[i]);
    return;
  }

  invoked_ichunk = update->ntimestep;

  // assign chunk IDs to atoms
  // will exclude atoms not in group or in optional region
  // already invoked if this is same timestep as last setup_chunks()

  if (update->ntimestep > invoked_setup) assign_chunk_ids();

  // compress chunk IDs via hash of the original uncompressed IDs
  // also apply discard rule except for binning styles which already did

  int nlocal = atom->nlocal;

  if (compress) {
    if (binflag) {
      for (i = 0; i < nlocal; i++) {
        if (exclude[i]) continue;
        if (hash->find(ichunk[i]) == hash->end()) exclude[i] = 1;
        else ichunk[i] = hash->find(ichunk[i])->second;
      }
    } else if (discard == NODISCARD) {
      for (i = 0; i < nlocal; i++) {
        if (exclude[i]) continue;
        if (hash->find(ichunk[i]) == hash->end()) ichunk[i] = nchunk;
        else ichunk[i] = hash->find(ichunk[i])->second;
      }
    } else {
      for (i = 0; i < nlocal; i++) {
        if (exclude[i]) continue;
        if (hash->find(ichunk[i]) == hash->end()) exclude[i] = 1;
        else ichunk[i] = hash->find(ichunk[i])->second;
      }
    }

  // else if no compression apply discard rule by itself

  } else {
    if (discard == NODISCARD) {
      for (i = 0; i < nlocal; i++) {
        if (exclude[i]) continue;
        if (ichunk[i] < 1 || ichunk[i] > nchunk) ichunk[i] = nchunk;;
      }
    } else {
      for (i = 0; i < nlocal; i++) {
        if (exclude[i]) continue;
        if (ichunk[i] < 1 || ichunk[i] > nchunk) exclude[i] = 1;
      }
    }
  }

  // set ichunk = 0 for excluded atoms
  // this should set any ichunk values which have not yet been set

  for (i = 0; i < nlocal; i++)
    if (exclude[i]) ichunk[i] = 0;

  // if newly calculated IDs need to persist, store them in fixstore
  // yes if idsflag = ONCE or idsflag = NFREQ and lock is in place

  if (idsflag == ONCE || (idsflag == NFREQ && lockfix)) { 
    double *vstore = fixstore->vstore;
    int nlocal = atom->nlocal;
    for (int i = 0; i < nlocal; i++) vstore[i] = ichunk[i];
  }

  // one-time check if which = MOLECULE and 
  // any chunks do not contain all atoms in the molecule

  if (molcheck) {
    check_molecules();
    molcheck = 0;
  }
}

/* ----------------------------------------------------------------------
   setup chunks
   return nchunk = # of chunks
     all atoms will be assigned a chunk ID from 1 to Nchunk, or 0
   also setup any internal state needed to quickly assign atoms to chunks
   called from compute_peratom() and also directly from
     fix ave/chunk and compute chunk commands
------------------------------------------------------------------------- */

int ComputeChunkAtom::setup_chunks()
{
  if (invoked_setup == update->ntimestep) return nchunk;

  // check if setup needs to be done
  // no if lock is in place
  // no if nchunkflag = ONCE, and already done once
  // otherwise yes
  // even if no, check if need to re-compute bin volumes
  //   so that fix ave/chunk can do proper density normalization

  int flag = 0;
  if (lockfix) flag = 1;
  if (nchunkflag == ONCE && invoked_setup >= 0) flag = 1;

  if (flag) {
    if (binflag && scaleflag == REDUCED && domain->box_change_size)
      bin_volumes();
    return nchunk;
  }

  invoked_setup = update->ntimestep;

  // assign chunk IDs to atoms
  // will exclude atoms not in group or in optional region
  // for binning styles, need to setup bins and their volume first
  // IDs are needed to scan for max ID and for compress()

  if (binflag) {
    nchunk = setup_bins();
    bin_volumes();
  }
  assign_chunk_ids();

  // set nchunk for chunk styles other than binning
  // for styles other than TYPE, scan for max ID

  if (which == TYPE) nchunk = atom->ntypes;
  else if (!binflag) {
    
    int nlocal = atom->nlocal;
    int hi = -1;
    for (int i = 0; i < nlocal; i++) {
      if (exclude[i]) continue;
      if (ichunk[i] > hi) hi = ichunk[i];
    }
      
    MPI_Allreduce(&hi,&nchunk,1,MPI_INT,MPI_MAX,world);
    if (nchunk <= 0) nchunk = 1;
  }

  // apply limit setting as well as compression of chunks with no atoms
  // if limit is set, there are 3 cases:
  //   no compression, limit specified before compression, or vice versa

  if (limit && !binflag) {
    if (!compress) {
      if (limitstyle == LIMITMAX) nchunk = MIN(nchunk,limit);
      else if (limitstyle == LIMITEXACT) nchunk = limit;
    } else if (limitfirst) {
      nchunk = MIN(nchunk,limit);
    }
  }

  if (compress) compress_chunk_ids();

  if (limit && !binflag && compress) {
    if (limitstyle == LIMITMAX) nchunk = MIN(nchunk,limit);
    else if (limitstyle == LIMITEXACT) nchunk = limit;
  }

  return nchunk;
}

/* ----------------------------------------------------------------------
   assign chunk IDs for all atoms, via ichunk vector
   except excluded atoms, their chunk IDs are set to 0 later
   also set exclude vector to 0/1 for all atoms
     excluded atoms are those not in group or in optional region
   called from compute_ichunk() and setup_chunks()
------------------------------------------------------------------------- */

void ComputeChunkAtom::assign_chunk_ids()
{
  int i;

  // grow integer chunk index vector if necessary

  if (atom->nlocal > nmaxint) {
    memory->destroy(ichunk);
    memory->destroy(exclude);
    nmaxint = atom->nmax;
    memory->create(ichunk,nmaxint,"chunk/atom:ichunk");
    memory->create(exclude,nmaxint,"chunk/atom:exclude");
  }

  // update region if necessary

  if (regionflag) region->prematch();

  // exclude = 1 if atom is not assigned to a chunk
  // exclude atoms not in group or not in optional region

  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (regionflag) {
    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit && 
          region->match(x[i][0],x[i][1],x[i][2])) exclude[i] = 0;
      else exclude[i] = 1;
    }
  } else {
    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) exclude[i] = 0;
      else exclude[i] = 1;
    }
  }

  // set ichunk to style value for included atoms
  // binning styles apply discard rule, others do not yet

  if (binflag) {
    if (which == BIN1D) atom2bin1d();
    else if (which == BIN2D) atom2bin2d();
    else if (which == BIN3D) atom2bin3d();

  } else if (which == TYPE) {
    int *type = atom->type;
    for (i = 0; i < nlocal; i++) {
      if (exclude[i]) continue;
      ichunk[i] = type[i];
    }

  } else if (which == MOLECULE) {
    tagint *molecule = atom->molecule;
    for (i = 0; i < nlocal; i++) {
      if (exclude[i]) continue;
      ichunk[i] = static_cast<int> (molecule[i]);
    }

  } else if (which == COMPUTE) {
    if (!(cchunk->invoked_flag & INVOKED_PERATOM)) {
      cchunk->compute_peratom();
      cchunk->invoked_flag |= INVOKED_PERATOM;
    }

    if (argindex == 0) {
      double *vec = cchunk->vector_atom;
      for (i = 0; i < nlocal; i++) {
        if (exclude[i]) continue;
        ichunk[i] = static_cast<int> (vec[i]);
      }
    } else {
      double **array = cchunk->array_atom;
      int argm1 = argindex - 1;
      for (i = 0; i < nlocal; i++) {
        if (exclude[i]) continue;
        ichunk[i] = static_cast<int> (array[i][argm1]);
      }
    }

  } else if (which == FIX) {
    if (update->ntimestep % fchunk->peratom_freq)
      error->all(FLERR,"Fix used in compute chunk/atom not "
                 "computed at compatible time");

    if (argindex == 0) {
      double *vec = fchunk->vector_atom;
      for (i = 0; i < nlocal; i++) {
        if (exclude[i]) continue;
        ichunk[i] = static_cast<int> (vec[i]);
      }
    } else {
      double **array = fchunk->array_atom;
      int argm1 = argindex - 1;
      for (i = 0; i < nlocal; i++) {
        if (exclude[i]) continue;
        ichunk[i] = static_cast<int> (array[i][argm1]);
      }
    }

  } else if (which == VARIABLE) {
    if (nlocal > maxvar) {
      maxvar = atom->nmax;
      memory->destroy(varatom);
      memory->create(varatom,maxvar,"chunk/atom:varatom");
    }

    input->variable->compute_atom(vchunk,igroup,varatom,1,0);
    for (i = 0; i < nlocal; i++) {
      if (exclude[i]) continue;
      ichunk[i] = static_cast<int> (varatom[i]);
    }
  }
}

/* ----------------------------------------------------------------------
   compress chunk IDs currently assigned to atoms across all processors
     by removing those with no atoms assigned
   current assignment excludes atoms not in group or in optional region
   current Nchunk = max ID
   operation:
     use hash to store list of populated IDs that I own
     add new IDs to populated lists communicated from all other procs
     final hash has global list of populated ideas
   reset Nchunk = length of global list
   called by setup_chunks() when setting Nchunk
   remapping of chunk IDs to smaller Nchunk occurs later in compute_ichunk()
------------------------------------------------------------------------- */

void ComputeChunkAtom::compress_chunk_ids()
{
  hash->clear();

  // put my IDs into hash

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    if (exclude[i]) continue;
    if (hash->find(ichunk[i]) == hash->end()) (*hash)[ichunk[i]] = 0;
  }

  // n = # of my populated IDs
  // nall = n summed across all procs

  int n = hash->size();
  bigint nbone = n;
  bigint nball;
  MPI_Allreduce(&nbone,&nball,1,MPI_LMP_BIGINT,MPI_SUM,world);
                            
  // create my list of populated IDs

  int *list = NULL;
  memory->create(list,n,"chunk/atom:list");

  n = 0;
  std::map<tagint,int>::iterator pos;
  for (pos = hash->begin(); pos != hash->end(); ++pos)
    list[n++] = pos->first;

  // if nall < 1M, just allgather all ID lists on every proc
  // else perform ring comm
  // add IDs from all procs to my hash

  if (nball <= IDMAX) {

    // setup for allgatherv

    int nprocs = comm->nprocs;
    int nall = nball;
    int *recvcounts,*displs,*listall;
    memory->create(recvcounts,nprocs,"chunk/atom:recvcounts");
    memory->create(displs,nprocs,"chunk/atom:displs");
    memory->create(listall,nall,"chunk/atom:listall");

    MPI_Allgather(&n,1,MPI_INT,recvcounts,1,MPI_INT,world);

    displs[0] = 0;
    for (int iproc = 1; iproc < nprocs; iproc++)
      displs[iproc] = displs[iproc-1] + recvcounts[iproc-1];
    
    // allgatherv acquires list of populated IDs from all procs
    
    MPI_Allgatherv(list,n,MPI_INT,listall,recvcounts,displs,MPI_INT,world);

    // add all unique IDs in listall to my hash

    for (int i = 0; i < nall; i++)
      if (hash->find(listall[i]) == hash->end()) (*hash)[listall[i]] = 0;

    // clean up

    memory->destroy(recvcounts);
    memory->destroy(displs);
    memory->destroy(listall);

  } else {
    cptr = this;
    comm->ring(n,sizeof(int),list,1,idring,NULL,0);
  }

  memory->destroy(list);

  // nchunk = length of hash containing populated IDs from all procs

  nchunk = hash->size();

  // reset hash value of each original chunk ID to ordered index
  //   ordered index = new compressed chunk ID (1 to Nchunk)
  //   leverages fact that map stores keys in ascending order
  // also allocate and set chunkID = list of original chunk IDs
  //   used by fix ave/chunk and compute property/chunk

  memory->destroy(chunkID);
  memory->create(chunkID,nchunk,"chunk/atom:chunkID");

  n = 0;
  for (pos = hash->begin(); pos != hash->end(); ++pos) {
    chunkID[n] = pos->first;
    (*hash)[pos->first] = ++n;
  }
}

/* ----------------------------------------------------------------------
   callback from comm->ring()
   cbuf = list of N chunk IDs from another proc
   loop over the list, add each to my hash
   hash ends up storing all unique IDs across all procs
------------------------------------------------------------------------- */

void ComputeChunkAtom::idring(int n, char *cbuf)
{
  tagint *list = (tagint *) cbuf;
  std::map<tagint,int> *hash = cptr->hash;
  for (int i = 0; i < n; i++) (*hash)[list[i]] = 0;
}

/* ----------------------------------------------------------------------
   one-time check for which = MOLECULE to check
     if each chunk contains all atoms in the molecule
   issue warning if not
   note that this check is without regard to discard rule
   if discard == NODISCARD, there is no easy way to check that all
     atoms in an out-of-bounds molecule were added to a chunk,
     some could have been excluded by group or region, others not
------------------------------------------------------------------------- */

void ComputeChunkAtom::check_molecules()
{
  tagint *molecule = atom->molecule;
  int nlocal = atom->nlocal;

  int flag = 0;

  if (!compress) {
    for (int i = 0; i < nlocal; i++) {
      if (molecule[i] > 0 && molecule[i] <= nchunk && 
          ichunk[i] == 0) flag = 1;
    }
  } else {
    int molid;
    for (int i = 0; i < nlocal; i++) {
      molid = static_cast<int> (molecule[i]);
      if (hash->find(molid) != hash->end() && ichunk[i] == 0) flag = 1;
    }
  }
  
  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall && comm->me == 0)
    error->warning(FLERR,
                   "One or more chunks do not contain all atoms in molecule");
}

/* ----------------------------------------------------------------------
   setup spatial bins and their extent and coordinates
   return nbins = # of bins, will become # of chunks
   called from setup_chunks()
------------------------------------------------------------------------- */

int ComputeChunkAtom::setup_bins()
{
  int i,j,k,m,n,idim;
  double lo,hi,coord1,coord2;

  // lo = bin boundary immediately below boxlo or minvalue
  // hi = bin boundary immediately above boxhi or maxvalue
  // allocate and initialize arrays based on new bin count

  double binlo[3],binhi[3];
  double *prd;
  if (scaleflag == REDUCED) {
    binlo[0] = domain->boxlo_lamda[0];
    binlo[1] = domain->boxlo_lamda[1];
    binlo[2] = domain->boxlo_lamda[2];
    binhi[0] = domain->boxhi_lamda[0];
    binhi[1] = domain->boxhi_lamda[1];
    binhi[2] = domain->boxhi_lamda[2];
    prd = domain->prd_lamda;
  } else {
    binlo[0] = domain->boxlo[0];
    binlo[1] = domain->boxlo[1];
    binlo[2] = domain->boxlo[2];
    binhi[0] = domain->boxhi[0];
    binhi[1] = domain->boxhi[1];
    binhi[2] = domain->boxhi[2];
    prd = domain->prd;
  }

  if (minflag[0] == COORD) binlo[0] = minvalue[0];
  if (minflag[1] == COORD) binlo[1] = minvalue[1];
  if (minflag[2] == COORD) binlo[2] = minvalue[2];
  if (maxflag[0] == COORD) binhi[0] = maxvalue[0];
  if (maxflag[1] == COORD) binhi[1] = maxvalue[1];
  if (maxflag[2] == COORD) binhi[2] = maxvalue[2];

  int nbins = 1;

  for (m = 0; m < ndim; m++) {
    idim = dim[m];
    if (originflag[m] == LOWER) origin[m] = binlo[idim];
    else if (originflag[m] == UPPER) origin[m] = binhi[idim];
    else if (originflag[m] == CENTER)
      origin[m] = 0.5 * (binlo[idim] + binhi[idim]);

    if (origin[m] < binlo[idim]) {
      n = static_cast<int> ((binlo[idim] - origin[m]) * invdelta[m]);
      lo = origin[m] + n*delta[m];
    } else {
      n = static_cast<int> ((origin[m] - binlo[idim]) * invdelta[m]);
      lo = origin[m] - n*delta[m];
      if (lo > binlo[idim]) lo -= delta[m];
    }
    if (origin[m] < binhi[idim]) {
      n = static_cast<int> ((binhi[idim] - origin[m]) * invdelta[m]);
      hi = origin[m] + n*delta[m];
      if (hi < binhi[idim]) hi += delta[m];
    } else {
      n = static_cast<int> ((origin[m] - binhi[idim]) * invdelta[m]);
      hi = origin[m] - n*delta[m];
    }

    if (lo > hi) error->all(FLERR,"Invalid bin bounds in compute chunk/atom");

    offset[m] = lo;
    nlayers[m] = static_cast<int> ((hi-lo) * invdelta[m] + 0.5);
    nbins *= nlayers[m];
    chunk_volume_scalar *= delta[m]/prd[idim];
  }

  // allocate and set bin coordinates

  memory->destroy(coord);
  memory->create(coord,nbins,ndim,"chunk/atom:coord");

  if (ndim == 1) {
    for (i = 0; i < nlayers[0]; i++)
      coord[i][0] = offset[0] + (i+0.5)*delta[0];
  } else if (ndim == 2) {
    m = 0;
    for (i = 0; i < nlayers[0]; i++) {
      coord1 = offset[0] + (i+0.5)*delta[0];
      for (j = 0; j < nlayers[1]; j++) {
        coord[m][0] = coord1;
        coord[m][1] = offset[1] + (j+0.5)*delta[1];
        m++;
      }
    }
  } else if (ndim == 3) {
    m = 0;
    for (i = 0; i < nlayers[0]; i++) {
      coord1 = offset[0] + (i+0.5)*delta[0];
      for (j = 0; j < nlayers[1]; j++) {
        coord2 = offset[1] + (j+0.5)*delta[1];
        for (k = 0; k < nlayers[2]; k++) {
          coord[m][0] = coord1;
          coord[m][1] = coord2;
          coord[m][2] = offset[2] + (k+0.5)*delta[2];
          m++;
        }
      }
    }
  }

  return nbins;
}

/* ----------------------------------------------------------------------
   calculate chunk volumes = bin volumes
   scalar if all bins have same volume
   vector if per-bin volumes are different
------------------------------------------------------------------------- */

void ComputeChunkAtom::bin_volumes()
{
  if (which == BIN1D || which == BIN2D || which == BIN3D) {
    if (domain->dimension == 3)
      chunk_volume_scalar = domain->xprd * domain->yprd * domain->zprd;
    else chunk_volume_scalar = domain->xprd * domain->yprd;
    double *prd;
    if (scaleflag == REDUCED) prd = domain->prd_lamda;
    else prd = domain->prd;
    for (int m = 0; m < ndim; m++)
      chunk_volume_scalar *= delta[m]/prd[dim[m]];

  } else {
    memory->destroy(chunk_volume_vec);
    memory->create(chunk_volume_vec,nchunk,"chunk/atom:chunk_volume_vec");
    // fill in the vector values
  }
}

/* ----------------------------------------------------------------------
   assign each atom to a 1d spatial bin (layer)
------------------------------------------------------------------------- */

void ComputeChunkAtom::atom2bin1d()
{
  int i,ibin;
  double *boxlo,*boxhi,*prd;
  double xremap;

  double **x = atom->x;
  int nlocal = atom->nlocal;

  int idim = dim[0];
  int nlayer1m1 = nlayers[0] - 1;
  int periodicity = domain->periodicity[idim];

  if (periodicity) {
    if (scaleflag == REDUCED) {
      boxlo = domain->boxlo_lamda;
      boxhi = domain->boxhi_lamda;
      prd = domain->prd_lamda;
    } else {
      boxlo = domain->boxlo;
      boxhi = domain->boxhi;
      prd = domain->prd;
    }
  }

  // remap each atom's relevant coord back into box via PBC if necessary
  // if scaleflag = REDUCED, box coords -> lamda coords
  // apply discard rule

  if (scaleflag == REDUCED) domain->x2lamda(nlocal);

  for (i = 0; i < nlocal; i++) {
    if (exclude[i]) continue;

    xremap = x[i][idim];
    if (periodicity) {
      if (xremap < boxlo[idim]) xremap += prd[idim];
      if (xremap >= boxhi[idim]) xremap -= prd[idim];
    }
      
    ibin = static_cast<int> ((xremap - offset[0]) * invdelta[0]);
    if (xremap < offset[0]) ibin--;
      
    if (discard == MIXED) {
      if (!minflag[idim]) ibin = MAX(ibin,0);
      else if (ibin < 0) {
        exclude[i] = 1;
        continue;
      }
      if (!maxflag[idim]) ibin = MIN(ibin,nlayer1m1);
      else if (ibin > nlayer1m1) {
        exclude[i] = 1;
        continue;
      }
    } else if (discard == NODISCARD) {
      ibin = MAX(ibin,0);
      ibin = MIN(ibin,nlayer1m1);
    } else if (ibin < 0 || ibin > nlayer1m1) {
      exclude[i] = 1;
      continue;
    }
      
    ichunk[i] = ibin+1;
  }

  if (scaleflag == REDUCED) domain->lamda2x(nlocal);
}

/* ----------------------------------------------------------------------
   assign each atom to a 2d spatial bin (pencil)
------------------------------------------------------------------------- */

void ComputeChunkAtom::atom2bin2d()
{
  int i,ibin,i1bin,i2bin;
  double *boxlo,*boxhi,*prd;
  double xremap,yremap;

  double **x = atom->x;
  int nlocal = atom->nlocal;

  int idim = dim[0];
  int jdim = dim[1];
  int nlayer1m1 = nlayers[0] - 1;
  int nlayer2m1 = nlayers[1] - 1;
  int *periodicity = domain->periodicity;

  if (periodicity[idim] || periodicity[jdim]) {
    if (scaleflag == REDUCED) {
      boxlo = domain->boxlo_lamda;
      boxhi = domain->boxhi_lamda;
      prd = domain->prd_lamda;
    } else {
      boxlo = domain->boxlo;
      boxhi = domain->boxhi;
      prd = domain->prd;
    }
  }

  // remap each atom's relevant coord back into box via PBC if necessary
  // if scaleflag = REDUCED, box coords -> lamda coords
  // apply discard rule

  if (scaleflag == REDUCED) domain->x2lamda(nlocal);

  for (i = 0; i < nlocal; i++) {
    if (exclude[i]) continue;

    xremap = x[i][idim];
    if (periodicity[idim]) {
      if (xremap < boxlo[idim]) xremap += prd[idim];
      if (xremap >= boxhi[idim]) xremap -= prd[idim];
    }
      
    i1bin = static_cast<int> ((xremap - offset[0]) * invdelta[0]);
    if (xremap < offset[0]) i1bin--;
      
    if (discard == MIXED) {
      if (!minflag[idim]) i1bin = MAX(i1bin,0);
      else if (i1bin < 0) {
        exclude[i] = 1;
        continue;
      }
      if (!maxflag[idim]) i1bin = MIN(i1bin,nlayer1m1);
      else if (i1bin > nlayer1m1) {
        exclude[i] = 1;
        continue;
      }
    } else if (discard == NODISCARD) {
      i1bin = MAX(i1bin,0);
      i1bin = MIN(i1bin,nlayer1m1);
    } else if (i1bin < 0 || i1bin > nlayer1m1) {
      exclude[i] = 1;
      continue;
    }

    yremap = x[i][jdim];
    if (periodicity[jdim]) {
      if (yremap < boxlo[jdim]) yremap += prd[jdim];
      if (yremap >= boxhi[jdim]) yremap -= prd[jdim];
    }

    i2bin = static_cast<int> ((yremap - offset[1]) * invdelta[1]);
    if (yremap < offset[1]) i2bin--;

    if (discard == MIXED) {
      if (!minflag[jdim]) i2bin = MAX(i2bin,0);
      else if (i2bin < 0) {
        exclude[i] = 1;
        continue;
      }
      if (!maxflag[jdim]) i2bin = MIN(i2bin,nlayer2m1);
      else if (i2bin > nlayer1m1) {
        exclude[i] = 1;
        continue;
      }
    } else if (discard == NODISCARD) {
      i2bin = MAX(i2bin,0);
      i2bin = MIN(i2bin,nlayer2m1);
    } else if (i2bin < 0 || i2bin > nlayer2m1) {
      exclude[i] = 1;
      continue;
    }

    ibin = i1bin*nlayers[1] + i2bin;
    ichunk[i] = ibin+1;
  }

  if (scaleflag == REDUCED) domain->lamda2x(nlocal);
}

/* ----------------------------------------------------------------------
   assign each atom to a 3d spatial bin (brick)
------------------------------------------------------------------------- */

void ComputeChunkAtom::atom2bin3d()
{
  int i,ibin,i1bin,i2bin,i3bin;
  double *boxlo,*boxhi,*prd;
  double xremap,yremap,zremap;

  double **x = atom->x;
  int nlocal = atom->nlocal;

  int idim = dim[0];
  int jdim = dim[1];
  int kdim = dim[2];
  int nlayer1m1 = nlayers[0] - 1;
  int nlayer2m1 = nlayers[1] - 1;
  int nlayer3m1 = nlayers[2] - 1;
  int *periodicity = domain->periodicity;

  if (periodicity[idim] || periodicity[jdim] || periodicity[kdim]) {
    if (scaleflag == REDUCED) {
      boxlo = domain->boxlo_lamda;
      boxhi = domain->boxhi_lamda;
      prd = domain->prd_lamda;
    } else {
      boxlo = domain->boxlo;
      boxhi = domain->boxhi;
      prd = domain->prd;
    }
  }

  // remap each atom's relevant coord back into box via PBC if necessary
  // if scaleflag = REDUCED, box coords -> lamda coords
  // apply discard rule

  if (scaleflag == REDUCED) domain->x2lamda(nlocal);

  for (i = 0; i < nlocal; i++) {
    if (exclude[i]) continue;

    xremap = x[i][idim];
    if (periodicity[idim]) {
      if (xremap < boxlo[idim]) xremap += prd[idim];
      if (xremap >= boxhi[idim]) xremap -= prd[idim];
    }

    i1bin = static_cast<int> ((xremap - offset[0]) * invdelta[0]);
    if (xremap < offset[0]) i1bin--;
    
    if (discard == MIXED) {
      if (!minflag[idim]) i1bin = MAX(i1bin,0);
      else if (i1bin < 0) {
        exclude[i] = 1;
        continue;
      }
      if (!maxflag[idim]) i1bin = MIN(i1bin,nlayer1m1);
      else if (i1bin > nlayer1m1) {
        exclude[i] = 1;
        continue;
      }
    } else if (discard == NODISCARD) {
      i1bin = MAX(i1bin,0);
      i1bin = MIN(i1bin,nlayer1m1);
    } else if (i1bin < 0 || i1bin > nlayer1m1) {
      exclude[i] = 1;
      continue;
    }
    
    yremap = x[i][jdim];
    if (periodicity[jdim]) {
      if (yremap < boxlo[jdim]) yremap += prd[jdim];
      if (yremap >= boxhi[jdim]) yremap -= prd[jdim];
    }

    i2bin = static_cast<int> ((yremap - offset[1]) * invdelta[1]);
    if (yremap < offset[1]) i2bin--;

    if (discard == MIXED) {
      if (!minflag[jdim]) i2bin = MAX(i2bin,0);
      else if (i2bin < 0) {
        exclude[i] = 1;
        continue;
      }
      if (!maxflag[jdim]) i2bin = MIN(i2bin,nlayer2m1);
      else if (i2bin > nlayer1m1) {
        exclude[i] = 1;
        continue;
      }
    } else if (discard == NODISCARD) {
      i2bin = MAX(i2bin,0);
      i2bin = MIN(i2bin,nlayer2m1);
    } else if (i2bin < 0 || i2bin > nlayer2m1) {
      exclude[i] = 1;
      continue;
    }

    zremap = x[i][kdim];
    if (periodicity[kdim]) {
      if (zremap < boxlo[kdim]) zremap += prd[kdim];
      if (zremap >= boxhi[kdim]) zremap -= prd[kdim];
    }

    i3bin = static_cast<int> ((zremap - offset[2]) * invdelta[2]);
    if (zremap < offset[2]) i3bin--;

    if (discard == MIXED) {
      if (!minflag[kdim]) i3bin = MAX(i3bin,0);
      else if (i3bin < 0) {
        exclude[i] = 1;
        continue;
      }
      if (!maxflag[kdim]) i3bin = MIN(i3bin,nlayer3m1);
      else if (i3bin > nlayer3m1) {
        exclude[i] = 1;
        continue;
      }
    } else if (discard == NODISCARD) {
      i3bin = MAX(i3bin,0);
      i3bin = MIN(i3bin,nlayer3m1);
    } else if (i3bin < 0 || i3bin > nlayer3m1) {
      exclude[i] = 1;
      continue;
    }
    
    ibin = i1bin*nlayers[1]*nlayers[2] + i2bin*nlayers[2] + i3bin;
    ichunk[i] = ibin+1;
  }

  if (scaleflag == REDUCED) domain->lamda2x(nlocal);
}

/* ----------------------------------------------------------------------
   process args for one dimension of binning info
------------------------------------------------------------------------- */

void ComputeChunkAtom::readdim(int narg, char **arg, int iarg, int idim)
{
  if (narg < iarg+3) error->all(FLERR,"Illegal compute chunk/atom command");
  if (strcmp(arg[iarg],"x") == 0) dim[idim] = 0;
  else if (strcmp(arg[iarg],"y") == 0) dim[idim] = 1;
  else if (strcmp(arg[iarg],"z") == 0) dim[idim] = 2;

  if (dim[idim] == 2 && domain->dimension == 2)
    error->all(FLERR,"Cannot use compute chunk/atom bin z for 2d model");

  if (strcmp(arg[iarg+1],"lower") == 0) originflag[idim] = LOWER;
  else if (strcmp(arg[iarg+1],"center") == 0) originflag[idim] = CENTER;
  else if (strcmp(arg[iarg+1],"upper") == 0) originflag[idim] = UPPER;
  else originflag[idim] = COORD;
  if (originflag[idim] == COORD) 
    origin[idim] = force->numeric(FLERR,arg[iarg+1]);
  
  delta[idim] = force->numeric(FLERR,arg[iarg+2]);
}

/* ----------------------------------------------------------------------
   initialize one atom's storage values, called when atom is created
   just set chunkID to 0 for new atom
------------------------------------------------------------------------- */

void ComputeChunkAtom::set_arrays(int i)
{
  if (!fixstore) return;
  double *vstore = fixstore->vstore;
  vstore[i] = 0.0;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays and per-chunk arrays
   note: nchunk is actually 0 until first call
------------------------------------------------------------------------- */

double ComputeChunkAtom::memory_usage()
{
  double bytes = 2*nmaxint * sizeof(int);          // ichunk,exclude
  bytes += nmax * sizeof(double);                  // chunk
  bytes += ncoord*nchunk * sizeof(double);         // coord
  if (compress) bytes += nchunk * sizeof(int);     // chunkID
  return bytes;
}
