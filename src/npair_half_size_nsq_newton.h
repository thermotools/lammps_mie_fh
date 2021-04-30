/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://lammps.sandia.gov/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef NPAIR_CLASS
// clang-format off
NPairStyle(half/size/nsq/newton,
           NPairHalfSizeNsqNewton,
           NP_HALF | NP_SIZE | NP_NSQ | NP_NEWTON | NP_ORTHO | NP_TRI);
// clang-format on
#else

#ifndef LMP_NPAIR_HALF_SIZE_NSQ_NEWTON_H
#define LMP_NPAIR_HALF_SIZE_NSQ_NEWTON_H

#include "npair.h"

namespace LAMMPS_NS {

class NPairHalfSizeNsqNewton : public NPair {
 public:
  NPairHalfSizeNsqNewton(class LAMMPS *);
  ~NPairHalfSizeNsqNewton() {}
  void build(class NeighList *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Neighbor list overflow, boost neigh_modify one

UNDOCUMENTED

*/
