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

/* ----------------------------------------------------------------------
   Contributing authors: Paul Crozier (SNL)
                         Carolyn Phillips (University of Michigan)
------------------------------------------------------------------------- */

#include "fix_ttm.h"

#include <cmath>
#include <cstring>
#include "atom.h"
#include "force.h"
#include "update.h"
#include "domain.h"
#include "respa.h"
#include "comm.h"
#include "random_mars.h"
#include "memory.h"
#include "error.h"

#include "tokenizer.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define MAXLINE 1024
//#define OFFSET 16384    // to avoid outside-of-box atoms being rounded incorrectly
//#define SHIFT 0.5       // 0.5 for nearest grid point, 0.0 for lower-left grid point

#define OFFSET 0    // to avoid outside-of-box atoms being rounded incorrectly
#define SHIFT 0.0       // 0.5 for nearest grid point, 0.0 for lower-left grid point

/* ---------------------------------------------------------------------- */

FixTTM::FixTTM(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  random(nullptr), 
  gfactor1(nullptr), gfactor2(nullptr), ratio(nullptr), flangevin(nullptr),
  T_electron(nullptr), T_electron_old(nullptr), 
  net_energy_transfer(nullptr), net_energy_transfer_all(nullptr)
{
  if (narg != 14) error->all(FLERR,"Illegal fix ttm command");

  vector_flag = 1;
  size_vector = 2;
  global_freq = 1;
  extvector = 1;
  nevery = 1;
  restart_peratom = 1;
  restart_global = 1;

  seed = utils::inumeric(FLERR,arg[3],false,lmp);
  electronic_specific_heat = utils::numeric(FLERR,arg[4],false,lmp);
  electronic_density = utils::numeric(FLERR,arg[5],false,lmp);
  electronic_thermal_conductivity = utils::numeric(FLERR,arg[6],false,lmp);
  gamma_p = utils::numeric(FLERR,arg[7],false,lmp);
  gamma_s = utils::numeric(FLERR,arg[8],false,lmp);
  v_0 = utils::numeric(FLERR,arg[9],false,lmp);
  nxgrid = utils::inumeric(FLERR,arg[10],false,lmp);
  nygrid = utils::inumeric(FLERR,arg[11],false,lmp);
  nzgrid = utils::inumeric(FLERR,arg[12],false,lmp);

  int n = strlen(arg[13]) + 1;
  infile = new char[n];
  strcpy(infile,arg[13]);

  // error check

  if (seed <= 0)
    error->all(FLERR,"Invalid random number seed in fix ttm command");
  if (electronic_specific_heat <= 0.0)
    error->all(FLERR,"Fix ttm electronic_specific_heat must be > 0.0");
  if (electronic_density <= 0.0)
    error->all(FLERR,"Fix ttm electronic_density must be > 0.0");
  if (electronic_thermal_conductivity < 0.0)
    error->all(FLERR,"Fix ttm electronic_thermal_conductivity must be >= 0.0");
  if (gamma_p <= 0.0) error->all(FLERR,"Fix ttm gamma_p must be > 0.0");
  if (gamma_s < 0.0) error->all(FLERR,"Fix ttm gamma_s must be >= 0.0");
  if (v_0 < 0.0) error->all(FLERR,"Fix ttm v_0 must be >= 0.0");
  if (nxgrid <= 0 || nygrid <= 0 || nzgrid <= 0)
    error->all(FLERR,"Fix ttm number of nodes must be > 0");

  v_0_sq = v_0*v_0;

  // OFFSET to make
  // SHIFT to map atom to nearest or lower-left grid point

  shift = OFFSET + SHIFT;

  // initialize Marsaglia RNG with processor-unique seed

  random = new RanMars(lmp,seed + comm->me);

  // allocate per-type arrays for force prefactors

  gfactor1 = new double[atom->ntypes+1];
  gfactor2 = new double[atom->ntypes+1];

  // check for allowed maxium number of total grid nodes

  bigint totalgrid = (bigint) nxgrid * nygrid * nzgrid;
  if (totalgrid > MAXSMALLINT)
    error->all(FLERR,"Too many grid points in fix ttm");
  ngridtotal = totalgrid;

  // allocate per-atom flangevin and zero it

  flangevin = nullptr;
  grow_arrays(atom->nmax);

  for (int i = 0; i < atom->nmax; i++) {
    flangevin[i][0] = 0.0;
    flangevin[i][1] = 0.0;
    flangevin[i][2] = 0.0;
  }

  // set 2 callbacks

  atom->add_callback(Atom::GROW);
  atom->add_callback(Atom::RESTART);

  // determines which class deallocate_grid() is called from

  deallocate_flag = 0;
}

/* ---------------------------------------------------------------------- */

FixTTM::~FixTTM()
{
  delete [] infile;

  delete random;

  delete [] gfactor1;
  delete [] gfactor2;

  memory->destroy(flangevin);

  if (!deallocate_flag) deallocate_grid();
}

/* ---------------------------------------------------------------------- */

void FixTTM::post_constructor()
{
  // allocate global grid on each proc
  // needs to be done in post_contructor() beccause is virtual method

  allocate_grid();

  // zero net_energy_transfer
  // in case compute_vector accesses it on timestep 0

  outflag = 0;
  memset(&net_energy_transfer[0][0][0],0,ngridtotal*sizeof(double));

  // set initial electron temperatures from user input file

  read_electron_temperatures(infile);
}

/* ---------------------------------------------------------------------- */

int FixTTM::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixTTM::init()
{
  if (domain->dimension == 2)
    error->all(FLERR,"Cannot use fix ttm with 2d simulation");
  if (domain->nonperiodic != 0)
    error->all(FLERR,"Cannot use non-periodic boundares with fix ttm");
  if (domain->triclinic)
    error->all(FLERR,"Cannot use fix ttm with triclinic box");

  // set force prefactors

  for (int i = 1; i <= atom->ntypes; i++) {
    gfactor1[i] = - gamma_p / force->ftm2v;
    gfactor2[i] =
      sqrt(24.0*force->boltz*gamma_p/update->dt/force->mvv2e) / force->ftm2v;
  }

  if (utils::strmatch(update->integrate_style,"^respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixTTM::setup(int vflag)
{
  if (utils::strmatch(update->integrate_style,"^verlet")) {
    post_force_setup(vflag);
  } else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa_setup(vflag,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }
}

/* ---------------------------------------------------------------------- */

void FixTTM::post_force_setup(int /*vflag*/)
{
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // apply langevin forces that have been stored from previous run

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      f[i][0] += flangevin[i][0];
      f[i][1] += flangevin[i][1];
      f[i][2] += flangevin[i][2];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixTTM::post_force(int /*vflag*/)
{
  int ix,iy,iz;
  double xscale,yscale,zscale;
  double gamma1,gamma2;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // apply damping and thermostat to all atoms in fix group

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      xscale = (x[i][0] - domain->boxlo[0])/domain->xprd;
      yscale = (x[i][1] - domain->boxlo[1])/domain->yprd;
      zscale = (x[i][2] - domain->boxlo[2])/domain->zprd;
      ix = static_cast<int>(xscale*nxgrid + shift) - OFFSET;
      iy = static_cast<int>(yscale*nygrid + shift) - OFFSET;
      iz = static_cast<int>(zscale*nzgrid + shift) - OFFSET;
      if (ix < 0) ix += nxgrid;
      if (iy < 0) iy += nygrid;
      if (iz < 0) iz += nzgrid;
      if (ix >= nxgrid) ix -= nxgrid;
      if (iy >= nygrid) iy -= nygrid;
      if (iz >= nzgrid) iz -= nzgrid;

      if (T_electron[iz][iy][ix] < 0)
        error->all(FLERR,"Electronic temperature dropped below zero");

      double tsqrt = sqrt(T_electron[iz][iy][ix]);

      gamma1 = gfactor1[type[i]];
      double vsq = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
      if (vsq > v_0_sq) gamma1 *= (gamma_p + gamma_s)/gamma_p;
      gamma2 = gfactor2[type[i]] * tsqrt;

      flangevin[i][0] = gamma1*v[i][0] + gamma2*(random->uniform()-0.5);
      flangevin[i][1] = gamma1*v[i][1] + gamma2*(random->uniform()-0.5);
      flangevin[i][2] = gamma1*v[i][2] + gamma2*(random->uniform()-0.5);

      f[i][0] += flangevin[i][0];
      f[i][1] += flangevin[i][1];
      f[i][2] += flangevin[i][2];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixTTM::post_force_respa_setup(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == nlevels_respa-1) post_force_setup(vflag);
}

/* ---------------------------------------------------------------------- */

void FixTTM::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixTTM::end_of_step()
{
  int ix,iy,iz;
  double xscale,yscale,zscale;

  double **x = atom->x;
  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  outflag = 0;
  for (iz = 0; iz < nzgrid; iz++)
    for (iy = 0; iy < nygrid; iy++)
      for (ix = 0; ix < nxgrid; ix++)
        net_energy_transfer[iz][iy][ix] = 0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      xscale = (x[i][0] - domain->boxlo[0])/domain->xprd;
      yscale = (x[i][1] - domain->boxlo[1])/domain->yprd;
      zscale = (x[i][2] - domain->boxlo[2])/domain->zprd;
      ix = static_cast<int>(xscale*nxgrid + shift) - OFFSET;
      iy = static_cast<int>(yscale*nygrid + shift) - OFFSET;
      iz = static_cast<int>(zscale*nzgrid + shift) - OFFSET;
      if (ix < 0) ix += nxgrid;
      if (iy < 0) iy += nygrid;
      if (iz < 0) iz += nzgrid;
      if (ix >= nxgrid) ix -= nxgrid;
      if (iy >= nygrid) iy -= nygrid;
      if (iz >= nzgrid) iz -= nzgrid;
      net_energy_transfer[iz][iy][ix] +=
        (flangevin[i][0]*v[i][0] + flangevin[i][1]*v[i][1] +
         flangevin[i][2]*v[i][2]);
    }

  MPI_Allreduce(&net_energy_transfer[0][0][0],
                &net_energy_transfer_all[0][0][0],
                ngridtotal,MPI_DOUBLE,MPI_SUM,world);

  double dx = domain->xprd/nxgrid;
  double dy = domain->yprd/nygrid;
  double dz = domain->zprd/nzgrid;
  double del_vol = dx*dy*dz;

  // num_inner_timesteps = # of inner steps (thermal solves)
  // required this MD step to maintain a stable explicit solve

  int num_inner_timesteps = 1;
  double inner_dt = update->dt;

  double stability_criterion = 1.0 -
    2.0*inner_dt/(electronic_specific_heat*electronic_density) *
    (electronic_thermal_conductivity*(1.0/dx/dx + 1.0/dy/dy + 1.0/dz/dz));

  if (stability_criterion < 0.0) {
    inner_dt = 0.5*(electronic_specific_heat*electronic_density) /
      (electronic_thermal_conductivity*(1.0/dx/dx + 1.0/dy/dy + 1.0/dz/dz));
    num_inner_timesteps = static_cast<int>(update->dt/inner_dt) + 1;
    inner_dt = update->dt/double(num_inner_timesteps);
    if (num_inner_timesteps > 1000000)
      error->warning(FLERR,"Too many inner timesteps in fix ttm");
  }

  // finite difference iterations to update T_electron

  for (int istep = 0; istep < num_inner_timesteps; istep++) {
    
    for (iz = 0; iz < nzgrid; iz++)
      for (iy = 0; iy < nygrid; iy++)
        for (ix = 0; ix < nxgrid; ix++)
          T_electron_old[iz][iy][ix] =
            T_electron[iz][iy][ix];

    // compute new electron T profile

    for (iz = 0; iz < nzgrid; iz++)
      for (iy = 0; iy < nygrid; iy++)
        for (ix = 0; ix < nxgrid; ix++) {
          int right_xnode = ix + 1;
          int right_ynode = iy + 1;
          int right_znode = iz + 1;
          if (right_xnode == nxgrid) right_xnode = 0;
          if (right_ynode == nygrid) right_ynode = 0;
          if (right_znode == nzgrid) right_znode = 0;
          int left_xnode = ix - 1;
          int left_ynode = iy - 1;
          int left_znode = iz - 1;
          if (left_xnode == -1) left_xnode = nxgrid - 1;
          if (left_ynode == -1) left_ynode = nygrid - 1;
          if (left_znode == -1) left_znode = nzgrid - 1;

          T_electron[iz][iy][ix] =
            T_electron_old[iz][iy][ix] +
            inner_dt/(electronic_specific_heat*electronic_density) *
            (electronic_thermal_conductivity *

             ((T_electron_old[iz][iy][right_xnode] +
               T_electron_old[iz][iy][left_xnode] -
               2*T_electron_old[iz][iy][ix])/dx/dx +
              (T_electron_old[iz][right_ynode][ix] +
               T_electron_old[iz][left_ynode][ix] -
               2*T_electron_old[iz][iy][ix])/dy/dy +
              (T_electron_old[right_znode][iy][iz] +
               T_electron_old[left_znode][iy][ix] -
               2*T_electron_old[iz][iy][ix])/dz/dz) -
             
             (net_energy_transfer_all[iz][iy][ix])/del_vol);
        }
  }
}

/* ----------------------------------------------------------------------
   read in initial electron temperatures from a user-specified file
   only called by proc 0
------------------------------------------------------------------------- */

void FixTTM::read_electron_temperatures(const char *filename)
{
  int ***T_initial_set;
  memory->create(T_initial_set,nxgrid,nygrid,nzgrid,"ttm:T_initial_set");
  memset(&T_initial_set[0][0][0],0,ngridtotal*sizeof(int));

  std::string name = utils::get_potential_file_path(filename);
  if (name.empty())
    error->one(FLERR,"Cannot open input file: {}",
                                 filename);
  FILE *fp = fopen(name.c_str(),"r");

  // read initial electron temperature values from file

  char line[MAXLINE];
  int ix,iy,iz;
  double T_tmp;

  while (1) {
    if (fgets(line,MAXLINE,fp) == nullptr) break;
    ValueTokenizer values(line);
    if (values.has_next()) ix = values.next_int();
    if (values.has_next()) iy = values.next_int();
    if (values.has_next()) iz = values.next_int();
    if (values.has_next()) T_tmp  = values.next_double();
    else error->one(FLERR,"Incorrect format in fix ttm input file");

    // check correctness of input data

    if ((ix < 0) || (ix >= nxgrid)
        || (iy < 0) || (iy >= nygrid)
        || (iz < 0) || (iz >= nzgrid))
      error->one(FLERR,"Fix ttm invalide node index in fix ttm input");

    if (T_tmp < 0.0)
      error->one(FLERR,"Fix ttm electron temperatures must be > 0.0");

    T_electron[iz][iy][ix] = T_tmp;
    T_initial_set[iz][iy][ix] = 1;
  }

  fclose(fp);

  // check completeness of input data

  for (int iz = 0; iz < nzgrid; iz++)
    for (int iy = 0; iy < nygrid; iy++)
      for (int ix = 0; ix < nxgrid; ix++)
        if (T_initial_set[iz][iy][ix] == 0)
          error->one(FLERR,"Initial temperatures not all set in fix ttm");

  memory->destroy(T_initial_set);

  MPI_Bcast(&T_electron[0][0][0],ngridtotal,MPI_DOUBLE,0,world);
}

/* ---------------------------------------------------------------------- */

void FixTTM::reset_dt()
{
  for (int i = 1; i <= atom->ntypes; i++)
    gfactor2[i] =
      sqrt(24.0*force->boltz*gamma_p/update->dt/force->mvv2e) / force->ftm2v;
}

/* ---------------------------------------------------------------------- */

void FixTTM::grow_arrays(int ngrow)
{
  memory->grow(flangevin,ngrow,3,"ttm:flangevin");
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixTTM::write_restart(FILE *fp)
{
  double *rlist;
  memory->create(rlist,nxgrid*nygrid*nzgrid+4,"ttm:rlist");

  int n = 0;
  rlist[n++] = nxgrid;
  rlist[n++] = nygrid;
  rlist[n++] = nzgrid;
  rlist[n++] = seed;

  // store global grid values

  for (int iz = 0; iz < nzgrid; iz++)
    for (int iy = 0; iy < nygrid; iy++)
      for (int ix = 0; ix < nxgrid; ix++)
        rlist[n++] =  T_electron[iz][iy][ix];

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(rlist,sizeof(double),n,fp);
  }

  memory->destroy(rlist);
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixTTM::restart(char *buf)
{
  int n = 0;
  double *rlist = (double *) buf;

  // check that restart grid size is same as current grid size

  int nxgrid_old = static_cast<int> (rlist[n++]);
  int nygrid_old = static_cast<int> (rlist[n++]);
  int nzgrid_old = static_cast<int> (rlist[n++]);

  if (nxgrid_old != nxgrid || nygrid_old != nygrid || nzgrid_old != nzgrid)
    error->all(FLERR,"Must restart fix ttm with same grid size");

  // change RN seed from initial seed, to avoid same Langevin factors
  // just increment by 1, since for RanMars that is a new RN stream

  seed = static_cast<int> (rlist[n++]) + 1;
  delete random;
  random = new RanMars(lmp,seed+comm->me);

  // restore global grid values

  for (int iz = 0; iz < nzgrid; iz++)
    for (int iy = 0; iy < nygrid; iy++)
      for (int ix = 0; ix < nxgrid; ix++)
        T_electron[iz][iy][ix] = rlist[n++];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixTTM::pack_restart(int i, double *buf)
{
  // pack buf[0] this way because other fixes unpack it

  buf[0] = 4;
  buf[1] = flangevin[i][0];
  buf[2] = flangevin[i][1];
  buf[3] = flangevin[i][2];
  return 4;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixTTM::unpack_restart(int nlocal, int nth)
{
  double **extra = atom->extra;

  // skip to Nth set of extra values
  // unpack the Nth first values this way because other fixes pack them

  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int> (extra[nlocal][m]);
  m++;

  flangevin[nlocal][0] = extra[nlocal][m++];
  flangevin[nlocal][1] = extra[nlocal][m++];
  flangevin[nlocal][2] = extra[nlocal][m++];
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixTTM::size_restart(int /*nlocal*/)
{
  return 4;
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixTTM::maxsize_restart()
{
  return 4;
}

/* ----------------------------------------------------------------------
   return the energy of the electronic subsystem or the net_energy transfer
   between the subsystems
------------------------------------------------------------------------- */

double FixTTM::compute_vector(int n)
{
  if (outflag == 0) {
    e_energy = 0.0;
    transfer_energy = 0.0;

    int ix,iy,iz;

    double dx = domain->xprd/nxgrid;
    double dy = domain->yprd/nygrid;
    double dz = domain->zprd/nzgrid;
    double del_vol = dx*dy*dz;

    for (iz = 0; iz < nzgrid; iz++)
      for (iy = 0; iy < nygrid; iy++)
        for (ix = 0; ix < nxgrid; ix++) {
          e_energy +=
            T_electron[iz][iy][ix]*electronic_specific_heat*
            electronic_density*del_vol;
          transfer_energy +=
            net_energy_transfer_all[iz][iy][ix]*update->dt;
        }

    outflag = 1;
  }

  if (n == 0) return e_energy;
  if (n == 1) return transfer_energy;
  return 0.0;
}

/* ----------------------------------------------------------------------
   memory usage for flangevin and 3d grid
------------------------------------------------------------------------- */

double FixTTM::memory_usage()
{
  double bytes = 0.0;
  bytes += (double)atom->nmax * 3 * sizeof(double);
  bytes += (double)4*ngridtotal * sizeof(int);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate 3d grid quantities
------------------------------------------------------------------------- */

void FixTTM::allocate_grid()
{
  memory->create(T_electron_old,nzgrid,nygrid,nxgrid,"ttm:T_electron_old");
  memory->create(T_electron,nzgrid,nygrid,nxgrid,"ttm:T_electron");
  memory->create(net_energy_transfer,nzgrid,nygrid,nxgrid,
                 "ttm:net_energy_transfer");
  memory->create(net_energy_transfer_all,nzgrid,nygrid,nxgrid,
                 "ttm:net_energy_transfer_all");
}

/* ----------------------------------------------------------------------
   deallocate 3d grid quantities
------------------------------------------------------------------------- */

void FixTTM::deallocate_grid()
{
  memory->destroy(T_electron_old);
  memory->destroy(T_electron);
  memory->destroy(net_energy_transfer);
  memory->destroy(net_energy_transfer_all);
}
