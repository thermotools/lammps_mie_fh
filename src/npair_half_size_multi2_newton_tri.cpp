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

#include "npair_half_size_multi2_newton_tri.h"
#include "neigh_list.h"
#include "atom.h"
#include "atom_vec.h"
#include "my_page.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NPairHalfSizeMulti2NewtonTri::NPairHalfSizeMulti2NewtonTri(LAMMPS *lmp) : NPair(lmp) {}

/* ----------------------------------------------------------------------
   size particles
   binned neighbor list construction with Newton's 3rd law for triclinic
   multi2-type stencil is itype-jtype dependent   
   each owned atom i checks its own bin and other bins in triclinic stencil
   every pair stored exactly once by some processor
------------------------------------------------------------------------- */

void NPairHalfSizeMulti2NewtonTri::build(NeighList *list)
{
  int i,j,k,n,itype,jtype,ibin,jbin,ns,js;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  double radi,radsum,cutdistsq;
  int *neighptr,*s;

  double **x = atom->x;
  double *radius = atom->radius;
  int *type = atom->type;
  int *mask = atom->mask;
  tagint *molecule = atom->molecule;
  int nlocal = atom->nlocal;
  if (includegroup) nlocal = atom->nfirst;

  int history = list->history;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  MyPage<int> *ipage = list->ipage;

  int mask_history = 3 << SBBITS;

  int inum = 0;
  ipage->reset();

  for (i = 0; i < nlocal; i++) {
    n = 0;
    neighptr = ipage->vget();

    itype = type[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];

    ibin = atom2bin_multi2[itype][i];

    // loop through stencils for all types
    for (jtype = 1; jtype <= atom->ntypes; jtype++) {

      // if same type use own bin
      if(itype == jtype) jbin = ibin;
	  else jbin = coord2bin(x[i], jtype);


      // loop over all atoms in bins in stencil
      // stencil is empty if i larger than j
      // stencil is half if i same size as j
      // stencil is full if i smaller than j
      // if half: pairs for atoms j "below" i are excluded
      // below = lower z or (equal z and lower y) or (equal zy and lower x)
      //         (equal zyx and j <= i)
      // latter excludes self-self interaction but allows superposed atoms

	  s = stencil_multi2[itype][jtype];
	  ns = nstencil_multi2[itype][jtype];
      
	  for (k = 0; k < ns; k++) {
	    js = binhead_multi2[jtype][jbin + s[k]];
	    for (j = js; j >= 0; j = bins_multi2[jtype][j]) {
                  
          // if same size (e.g. same type), use half stencil            
          if(cutneighsq[itype][itype] == cutneighsq[jtype][jtype]){
            if (x[j][2] < ztmp) continue;
            if (x[j][2] == ztmp) {
              if (x[j][1] < ytmp) continue;
              if (x[j][1] == ytmp) {
                if (x[j][0] < xtmp) continue;
                if (x[j][0] == xtmp && j <= i) continue;
              }
            }                
          }  
          
          if (exclude && exclusion(i,j,itype,jtype,mask,molecule)) continue;
        
	      delx = xtmp - x[j][0];
	      dely = ytmp - x[j][1];
	      delz = ztmp - x[j][2];
	      rsq = delx*delx + dely*dely + delz*delz;
	      radsum = radi + radius[j];
	      cutdistsq = (radsum+skin) * (radsum+skin);
      
	      if (rsq <= cutdistsq) {
	        if (history && rsq < radsum*radsum) 
	  	    neighptr[n++] = j ^ mask_history;
	        else
	  	    neighptr[n++] = j;
	      }
	    }
	  }
    }

    ilist[inum++] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    ipage->vgot(n);
    if (ipage->status())
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
  }

  list->inum = inum;
}
