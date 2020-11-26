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

#include "nstencil_full_multi2_3d.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "nbin.h"
#include "memory.h"
#include "atom.h"
#include <math.h>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NStencilFullMulti23d::NStencilFullMulti23d(LAMMPS *lmp) : NStencil(lmp) {}

/* ---------------------------------------------------------------------- */

void NStencilFullMulti23d::set_stencil_properties()
{
  int n = atom->ntypes;
  int i, j;
  
  // like -> like => use half stencil
  for (i = 1; i <= n; i++) {
    stencil_half[i][i] = 0;
    stencil_skip[i][i] = 0;
    stencil_bin_type[i][i] = i;
    stencil_cutsq[i][i] = cutneighsq[i][i];
  }

  // smaller -> larger => use existing newtoff stencil in larger bin
  // larger -> smaller => use multi-like stencil for small-large in smaller bin
  // If types are same cutoff, use existing like-like stencil.

  for (i = 1; i <= n; i++) {
    for (j = 1; j <= n; j++) {
      if(i == j) continue;

      stencil_half[i][j] = 0;
      stencil_skip[i][j] = 0;
      
      if(cutneighsq[i][i] <= cutneighsq[j][j]){
        stencil_cutsq[i][j] = cutneighsq[j][j];   
        stencil_bin_type[i][j] = j;
      } else {
        stencil_cutsq[i][j] = cutneighsq[i][j];   
        stencil_bin_type[i][j] = j;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   create stencils based on bin geometry and cutoff
------------------------------------------------------------------------- */

void NStencilFullMulti23d::create()
{
  int itype, jtype, bin_type, i, j, k, ns;
  int n = atom->ntypes;
  double cutsq;
  
  
  for (itype = 1; itype <= n; itype++) {
    for (jtype = 1; jtype <= n; jtype++) {
      if (stencil_skip[itype][jtype]) continue;
      
      ns = 0;
      
      sx = stencil_sx_multi2[itype][jtype];
      sy = stencil_sy_multi2[itype][jtype];
      sz = stencil_sz_multi2[itype][jtype];
      
      mbinx = stencil_mbinx_multi2[itype][jtype];
      mbiny = stencil_mbiny_multi2[itype][jtype];
      mbinz = stencil_mbinz_multi2[itype][jtype];  
      
      bin_type = stencil_bin_type[itype][jtype];
      
      cutsq = stencil_cutsq[itype][jtype];
      
      for (k = -sz; k <= sz; k++)
        for (j = -sy; j <= sy; j++)
          for (i = -sx; i <= sx; i++)
	        if (bin_distance_multi2(i,j,k,bin_type) < cutsq)
	          stencil_multi2[itype][jtype][ns++] = 
                      k*mbiny*mbinx + j*mbinx + i;
      
      nstencil_multi2[itype][jtype] = ns;
    }
  }
}
