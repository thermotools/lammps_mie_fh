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

#ifndef LMP_TABULAR_FUNCTION_H
#define LMP_TABULAR_FUNCTION_H

namespace LAMMPS_NS {
class TabularFunction {
 public:
  TabularFunction();
  virtual ~TabularFunction();

  void set_values(int, double, double, double *);

 private:
  int size;
  double xmin, xmax, xmaxsq, rdx, vmax;
  double *xs, *ys, *ys1, *ys2, *ys3, *ys4, *ys5, *ys6;

  void set_xrange(double x1, double x2);
  void reset_size(int);
  void initialize();

 public:
  // used by pair style bop
  void value(double x, double &y, int ny, double &y1, int ny1)
  {
    double ps = (x - xmin) * rdx + 1.0;
    int ks = ps;
    if (ks > size - 1) ks = size - 1;
    ps = ps - ks;
    if (ps > 1.0) ps = 1.0;
    if (ny)
      y = ((ys3[ks - 1] * ps + ys2[ks - 1]) * ps + ys1[ks - 1]) * ps +
          ys[ks - 1];
    if (ny1) y1 = (ys6[ks - 1] * ps + ys5[ks - 1]) * ps + ys4[ks - 1];
  }
  // used by pair style polymorphic
  void value2(double x, double &y, int ny, double &y1, int ny1)
  {
    double ps = (x - xmin) * rdx;
    int ks = ps + 0.5;
    if (ks > size-1) ks = size-1;
    if (ks < 0 ) ks = 0;
    ps = ps - ks;
    if (ny) y = ((ys3[ks]*ps + ys2[ks])*ps + ys1[ks])*ps + ys[ks];
    if (ny1) y1 = (ys6[ks]*ps + ys5[ks])*ps + ys4[ks];
  }

  double get_xmin() const { return xmin; }
  double get_xmax() const { return xmax; }
  double get_xmaxsq() const { return xmaxsq; }
  double get_rdx() const { return rdx; }
  double get_vmax() { return vmax; }
};
}    // namespace LAMMPS_NS

#endif
