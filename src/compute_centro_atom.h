/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(centro/atom,ComputeCentroAtom);
// clang-format on
#else

#ifndef COMPUTE_CENTRO_ATOM_H
#define COMPUTE_CENTRO_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeCentroAtom : public Compute {
 public:
  ComputeCentroAtom(class LAMMPS *, int, char **);
  ~ComputeCentroAtom() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void compute_peratom() override;
  double memory_usage() override;

 private:
  int nmax, maxneigh, nnn;
  double *distsq;
  int *nearest;
  class NeighList *list;
  double *centro;
  int axes_flag;

  void select(int, int, double *);
  void select2(int, int, double *, int *);
};

}    // namespace LAMMPS_NS

#endif
#endif
