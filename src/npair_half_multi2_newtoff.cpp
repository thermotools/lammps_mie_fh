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

#include "npair_half_multi2_newtoff.h"
#include "neigh_list.h"
#include "atom.h"
#include "atom_vec.h"
#include "molecule.h"
#include "domain.h"
#include "my_page.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NPairHalfMulti2Newtoff::NPairHalfMulti2Newtoff(LAMMPS *lmp) : NPair(lmp) {}

/* ----------------------------------------------------------------------
   REWRITE
   binned neighbor list construction with partial Newton's 3rd law
   each owned atom i checks own bin and other bins in stencil
   multi-type stencil is itype dependent and is distance checked
   pair stored once if i,j are both owned and i < j
   pair stored by me if j is ghost (also stored by proc owning j)
------------------------------------------------------------------------- */

void NPairHalfMulti2Newtoff::build(NeighList *list)
{
  int i,j,k,n,itype,jtype,ktype,ibin,kbin,which,ns,imol,iatom,moltemplate;
  tagint tagprev;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *neighptr,*s;
  int js;

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  tagint *tag = atom->tag;
  tagint *molecule = atom->molecule;
  tagint **special = atom->special;
  int **nspecial = atom->nspecial;
  int nlocal = atom->nlocal;
  if (includegroup) nlocal = atom->nfirst;

  int *molindex = atom->molindex;
  int *molatom = atom->molatom;
  Molecule **onemols = atom->avec->onemols;
  if (molecular == Atom::TEMPLATE) moltemplate = 1;
  else moltemplate = 0;

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  MyPage<int> *ipage = list->ipage;

  int inum = 0;
  ipage->reset();

  for (i = 0; i < nlocal; i++) {
    n = 0;
    neighptr = ipage->vget();

    itype = type[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    if (moltemplate) {
      imol = molindex[i];
      iatom = molatom[i];
      tagprev = tag[i] - iatom - 1;
    }

    // loop over all atoms in other bins in stencil including self
    // only store pair if i < j
    // skip if i,j neighbor cutoff is less than bin distance
    // stores own/own pairs only once
    // stores own/ghost pairs on both procs

    ibin = nb->atom2bin_type[itype][i];
    for (ktype = 1; ktype <= atom->ntypes; ktype++) {
      if (itype == ktype) {
	    kbin = ibin;
      } else {
	    // Locate i in ktype bin
	    kbin = coord2bin(x[i], ktype);
      }
      
      s = stencil_multi2[itype][ktype];
      ns = nstencil_multi2[itype][ktype];
      for (k = 0; k < ns; k++) {
	    js = binhead_multi2[ktype][kbin + s[k]];
	    for (j = js; j >=0; j = nb->bins_multi2[ktype][j]) {
	      if (j <= i) continue;
          
	      jtype = type[j];
 
	      if (exclude && exclusion(i,j,itype,jtype,mask,molecule)) continue;
          
	      delx = xtmp - x[j][0];
	      dely = ytmp - x[j][1];
	      delz = ztmp - x[j][2];
          rsq = delx*delx + dely*dely + delz*delz;
          
          if (rsq <= cutneighsq[itype][jtype]) {
            if (molecular != Atom::ATOMIC) {
              if (!moltemplate)
                which = find_special(special[i],nspecial[i],tag[j]);
              else if (imol >= 0)
                which = find_special(onemols[imol]->special[iatom],
                                     onemols[imol]->nspecial[iatom],
                                     tag[j]-tagprev);
              else which = 0;
              if (which == 0) neighptr[n++] = j;
              else if (domain->minimum_image_check(delx,dely,delz))
                neighptr[n++] = j;
              else if (which > 0) neighptr[n++] = j ^ (which << SBBITS);
            } else neighptr[n++] = j;
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
