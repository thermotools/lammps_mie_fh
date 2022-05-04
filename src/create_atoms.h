/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMMAND_CLASS
// clang-format off
CommandStyle(create_atoms,CreateAtoms);
// clang-format on
#else

#ifndef LMP_CREATE_ATOMS_H
#define LMP_CREATE_ATOMS_H

#include "command.h"

namespace LAMMPS_NS {

class CreateAtoms : public Command {
 public:
  CreateAtoms(class LAMMPS *);
  void command(int, char **) override;

 private:
  int ntype, style, mode, nbasis, nrandom, seed;
  int remapflag;
  int maxtry;
  int quat_user;
  int overlapflag;
  double overlap_radius;
  int subsetflag;
  bigint nsubset;
  double subsetfrac;
  int *basistype;
  double xone[3], quatone[4];

  int varflag, vvar, xvar, yvar, zvar;
  char *vstr, *xstr, *ystr, *zstr;
  char *xstr_copy, *ystr_copy, *zstr_copy;

  int ilo, ihi, jlo, jhi, klo, khi;

  int nlatt;             // number of owned lattice sites
  int nlatt_overflow;    // 1 if local nlatt exceeds a 32-bit int

  int *flag;    // flag subset of particles to insert on lattice
  int *next;

  class Region *region;
  class Molecule *onemol;
  class RanMars *ranmol;
  class RanMars *ranlatt;
  double **temp_mol_coords;

  int triclinic;
  double sublo[3], subhi[3];    // epsilon-extended proc sub-box for adding atoms

  void add_single();
  void add_random();
  void add_lattice();
  void loop_lattice(int);
  void gen_mol_coords(double *, double * = nullptr);
  void create_mol();
  int vartest(double *);    // evaluate a variable with new atom position
};

}    // namespace LAMMPS_NS

#endif
#endif
