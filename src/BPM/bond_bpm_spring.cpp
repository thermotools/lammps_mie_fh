// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "bond_bpm_spring.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix_bond_history.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"

#define EPSILON 1e-10

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

BondBPMSpring::BondBPMSpring(LAMMPS *lmp) : BondBPM(lmp)
{
  partial_flag = 1;
}

/* ---------------------------------------------------------------------- */

BondBPMSpring::~BondBPMSpring()
{
  if (fix_bond_history) modify->delete_fix("BOND_HISTORY_BPM_SPRING");

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(k);
    memory->destroy(ecrit);
    memory->destroy(gamma);
  }
}

/* ----------------------------------------------------------------------
  Store data for a single bond - if bond added after LAMMPS init (e.g. pour)
------------------------------------------------------------------------- */

double BondBPMSpring::store_bond(int n,int i,int j)
{
  int m,k;
  double delx, dely, delz, r;
  double **x = atom->x;
  double **bondstore = fix_bond_history->bondstore;
  tagint *tag = atom->tag;

  delx = x[i][0] - x[j][0];
  dely = x[i][1] - x[j][1];
  delz = x[i][2] - x[j][2];

  r = sqrt(delx*delx + dely*dely + delz*delz);
  bondstore[n][0] = r;

  if (i < atom->nlocal) {
    for (m = 0; m < atom->num_bond[i]; m ++) {
      if (atom->bond_atom[i][m] == tag[j]) {
        fix_bond_history->update_atom_value(i, m, 0, r);
      }
    }
  }

  if (j < atom->nlocal) {
    for (m = 0; m < atom->num_bond[j]; m ++) {
      if (atom->bond_atom[j][m] == tag[i]) {
        fix_bond_history->update_atom_value(j, m, 0, r);
      }
    }
  }

  return r;
}

/* ----------------------------------------------------------------------
  Store data for all bonds called once
------------------------------------------------------------------------- */

void BondBPMSpring::store_data()
{
  int i, j, m, type;
  double delx, dely, delz, r;
  double **x = atom->x;
  int **bond_type = atom->bond_type;

  for (i = 0; i < atom->nlocal; i ++) {
    for (m = 0; m < atom->num_bond[i]; m ++) {
      type = bond_type[i][m];

      //Skip if bond was turned off
      if (type < 0) continue;

      // map to find index n
      j = atom->map(atom->bond_atom[i][m]);
      if (j == -1) error->one(FLERR, "Atom missing in BPM bond");

      delx = x[i][0] - x[j][0];
      dely = x[i][1] - x[j][1];
      delz = x[i][2] - x[j][2];

      // Get closest image in case bonded with ghost
      domain->minimum_image(delx, dely, delz);
      r = sqrt(delx*delx + dely*dely + delz*delz);

      fix_bond_history->update_atom_value(i, m, 0, r);
    }
  }

  fix_bond_history->post_neighbor();
}

/* ---------------------------------------------------------------------- */

void BondBPMSpring::compute(int eflag, int vflag)
{

  if (! fix_bond_history->stored_flag) {
    fix_bond_history->stored_flag = true;
    store_data();
  }

  int i1,i2,itmp,m,n,type,itype,jtype;
  double delx, dely, delz, delvx, delvy, delvz;
  double e, rsq, r, r0, rinv, smooth, fbond, ebond, dot;

  ev_init(eflag,vflag);

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  tagint *tag = atom->tag;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  double **bondstore = fix_bond_history->bondstore;

  for (n = 0; n < nbondlist; n++) {

    // skip bond if already broken
    if (bondlist[n][2] <= 0) continue;

    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    type = bondlist[n][2];
    r0 = bondstore[n][0];

    // Ensure pair is always ordered to ensure numerical operations
    // are identical to minimize the possibility that a bond straddling
    // an mpi grid (newton off) doesn't break on one proc but not the other
    if (tag[i2] < tag[i1]) {
      itmp = i1;
      i1 = i2;
      i2 = itmp;
    }

    // If bond hasn't been set - should be initialized to zero
    if (r0 < EPSILON || std::isnan(r0))
      r0 = store_bond(n,i1,i2);

    delx = x[i1][0] - x[i2][0];
    dely = x[i1][1] - x[i2][1];
    delz = x[i1][2] - x[i2][2];

    rsq = delx*delx + dely*dely + delz*delz;
    r = sqrt(rsq);
    e = (r - r0)/r0;

    if (fabs(e) > ecrit[type]) {
      bondlist[n][2] = 0;
      process_broken(i1, i2);
      continue;
    }

    rinv = 1.0/r;
    fbond = k[type]*(r0-r);
    if (eflag) ebond = -0.5*fbond*(r0-r);

    delvx = v[i1][0] - v[i2][0];
    delvy = v[i1][1] - v[i2][1];
    delvz = v[i1][2] - v[i2][2];
    dot = delx*delvx + dely*delvy + delz*delvz;
    fbond -= gamma[type]*dot*rinv;

    smooth = (r-r0)/(r0*ecrit[type]);
    smooth *= smooth;
    smooth *= smooth;
    smooth *= smooth;
    smooth = 1 - smooth;
    
    fbond *= rinv*smooth;

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += delx*fbond;
      f[i1][1] += dely*fbond;
      f[i1][2] += delz*fbond;
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] -= delx*fbond;
      f[i2][1] -= dely*fbond;
      f[i2][2] -= delz*fbond;
    }

    if (evflag) ev_tally(i1,i2,nlocal,newton_bond,ebond,fbond,delx,dely,delz);
  }
}

/* ---------------------------------------------------------------------- */

void BondBPMSpring::allocate()
{
  allocated = 1;
  int n = atom->nbondtypes;

  memory->create(k,n+1,"bond:k");
  memory->create(ecrit,n+1,"bond:ecrit");
  memory->create(gamma,n+1,"bond:gamma");

  memory->create(setflag,n+1,"bond:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void BondBPMSpring::coeff(int narg, char **arg)
{
  if (narg != 4) error->all(FLERR,"Incorrect args for bond coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  utils::bounds(FLERR,arg[0],1,atom->nbondtypes,ilo,ihi,error);

  double k_one = utils::numeric(FLERR,arg[1],false,lmp);
  double ecrit_one = utils::numeric(FLERR,arg[2],false,lmp);
  double gamma_one = utils::numeric(FLERR,arg[3],false,lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    ecrit[i] = ecrit_one;
    gamma[i] = gamma_one;
    setflag[i] = 1;
    count++;

    if (1.0+ecrit[i] > max_stretch) max_stretch = 1.0+ecrit[i];
  }

  if (count == 0) error->all(FLERR,"Incorrect args for bond coefficients");
}

/* ----------------------------------------------------------------------
   check if pair defined and special_bond settings are valid
------------------------------------------------------------------------- */

void BondBPMSpring::init_style()
{
  BondBPM::init_style();

  if (comm->ghost_velocity == 0)
    error->all(FLERR,"Bond bpm/spring requires ghost atoms store velocity");

  if (!fix_bond_history)
    fix_bond_history = (FixBondHistory *) modify->add_fix(
         "BOND_HISTORY_BPM_SPRING all BOND_HISTORY 0 1");
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void BondBPMSpring::write_restart(FILE *fp)
{
  fwrite(&k[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&ecrit[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&gamma[1],sizeof(double),atom->nbondtypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void BondBPMSpring::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    utils::sfread(FLERR,&k[1],sizeof(double),atom->nbondtypes,fp,nullptr,error);
    utils::sfread(FLERR,&ecrit[1],sizeof(double),atom->nbondtypes,fp,nullptr,error);
    utils::sfread(FLERR,&gamma[1],sizeof(double),atom->nbondtypes,fp,nullptr,error);
  }
  MPI_Bcast(&k[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&ecrit[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&gamma[1],atom->nbondtypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void BondBPMSpring::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nbondtypes; i++)
    fprintf(fp,"%d %g %g %g\n", i,k[i],ecrit[i],gamma[i]);
}

/* ---------------------------------------------------------------------- */

double BondBPMSpring::single(int type, double rsq, int i, int j,
                           double &fforce)
{
  // Not yet enabled
  if (type <= 0) return 0.0;

  //double r0;
  //for (int n = 0; n < atom->num_bond[i]; n ++) {
  //  if (atom->bond_atom[i][n] == atom->tag[j]) {
  //    r0 = fix_bond_history->get_atom_value(i, n, 0);
  //  }
  //}
}
