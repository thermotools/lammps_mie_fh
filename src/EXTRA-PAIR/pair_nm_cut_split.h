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

#ifdef PAIR_CLASS
PairStyle(nm/cut/split,PairNMCutSplit)
#else

#ifndef LMP_PAIR_NM_CUT_SPLIT_H
#define LMP_PAIR_NM_CUT_SPLIT_H

#include "pair_nm_cut.h"
namespace LAMMPS_NS {

class PairNMCutSplit : public PairNMCut {
    public :
    PairNMCutSplit(class LAMMPS *);
    double single(int, int, int, int, double, double, double, double &);
    virtual void compute(int, int);
    };
}
#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

*/
