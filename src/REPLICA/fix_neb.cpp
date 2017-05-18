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

#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "fix_neb.h"
#include "universe.h"
#include "update.h"
#include "atom.h"
#include "domain.h"
#include "comm.h"
#include "modify.h"
#include "compute.h"
#include "group.h"
#include "memory.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{SINGLE_PROC_DIRECT,SINGLE_PROC_MAP,MULTI_PROC};
/* ---------------------------------------------------------------------- */

FixNEB::FixNEB(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), id_pe(NULL), pe(NULL), xprev(NULL), xnext(NULL), fnext(NULL), 
  tangent(NULL),springF(NULL), xsend(NULL), xrecv(NULL),fsend(NULL),frecv(NULL), tagsend(NULL), tagrecv(NULL),   xsendall(NULL), xrecvall(NULL),fsendall(NULL), frecvall(NULL), tagsendall(NULL), tagrecvall(NULL),   counts(NULL), displacements(NULL),nlenall(NULL)
{


  StandardNEB=false;
  NEBLongRange=true;
  PerpSpring=false;
  FreeEndIni=false;
  FreeEndFinal=false;
  FreeEndFinalWithRespToEIni =false;
  FinalAndInterWithRespToEIni = false;
  kspringPerp=0;
  if (narg < 4) error->all(FLERR,"Illegal fix neb command, argument missing");

  kspring = force->numeric(FLERR,arg[3]);
  if (kspring <= 0.0) error->all(FLERR,"Illegal fix neb command. The spring force was not provided properly");

  int iarg =4;
  while (iarg < narg){
    if (strcmp (arg[iarg],"idealpos")==0) 
      {NEBLongRange = true;
	iarg+=1;}
    else if (strcmp (arg[iarg],"nearestneigh")==0)
      {NEBLongRange = false;
	StandardNEB = true;
	iarg+=1;}
    else if (strcmp (arg[iarg],"PerpSpring")==0) 
      {PerpSpring=true;
	kspringPerp = force->numeric(FLERR,arg[iarg+1]);
	if (kspringPerp < 0.0) error->all(FLERR,"Illegal fix neb command. The perpendicular spring force was not provided properly");
      iarg+=2;
      }
    else if (strcmp (arg[iarg],"freeend")==0) 
      {
	if (strcmp (arg[iarg+1],"ini")==0) 
	  FreeEndIni=true;
	else if (strcmp (arg[iarg+1],"final")==0) 
	  FreeEndFinal=true;
	else if (strcmp (arg[iarg+1],"final")==0) 
	  FreeEndFinal=true;
	else if (strcmp (arg[iarg+1],"finalWithRespToIni")==0) 
	  FreeEndFinalWithRespToEIni=true;
	else if (strcmp (arg[iarg+1],"finalAndInterWithRespToIni")==0) 
	  {FinalAndInterWithRespToEIni=true;
	    FreeEndFinalWithRespToEIni=true;}
	iarg+=2;}
    else {error->all(FLERR,"Illegal fix neb command. Unknown keyword");}
  }

  // nreplica = number of partitions
  // ireplica = which world I am in universe
  // nprocs_universe = # of procs in all replicase
  // procprev,procnext = root proc in adjacent replicas


  me = comm->me;
  nprocs = comm->nprocs;

  nprocs_universe = universe->nprocs;
  nreplica = universe->nworlds;
  ireplica = universe->iworld;
  
  if (ireplica > 0) procprev = universe->root_proc[ireplica-1];
  else procprev = -1;
  if (ireplica < nreplica-1) procnext = universe->root_proc[ireplica+1];
  else procnext = -1;
  uworld = universe->uworld;
  int *iroots = new int[nreplica];
  MPI_Group uworldgroup,rootgroup;
  if (NEBLongRange ){
    for (int iIm =0; iIm < nreplica;iIm++)
      {
	iroots[iIm]=universe->root_proc[iIm];
      }
    MPI_Comm_group(uworld, &uworldgroup);
    MPI_Group_incl(uworldgroup, nreplica, iroots, &rootgroup);
    
    //    MPI_Comm_create_group(uworld, rootgroup, 0, &rootworld);
    MPI_Comm_create(uworld, rootgroup, &rootworld);
  }
  // create a new compute pe style
  // id = fix-ID + pe, compute group = all

  int n = strlen(id) + 4;
  id_pe = new char[n];
  strcpy(id_pe,id);
  strcat(id_pe,"_pe");

  char **newarg = new char*[3];
  newarg[0] = id_pe;
  newarg[1] = (char *) "all";
  newarg[2] = (char *) "pe";
  modify->add_compute(3,newarg);
  delete [] newarg;

  // initialize local storage

  maxlocal = 0;
  ntotal = 0;

  xprev = xnext = tangent = springF= NULL;
  fnext = NULL;
  xsend = xrecv = NULL;
  fsend = frecv = NULL;
  tagsend = tagrecv = NULL;
  xsendall = xrecvall = NULL;
  fsendall = frecvall = NULL;
  tagsendall = tagrecvall = NULL;
  counts = displacements = NULL;


  if (NEBLongRange)
    {nlenall=NULL;}
}

/* ---------------------------------------------------------------------- */

FixNEB::~FixNEB()
{
  modify->delete_compute(id_pe);
  delete [] id_pe;

  memory->destroy(xprev);
  memory->destroy(xnext);
  memory->destroy(tangent);
  memory->destroy(springF);
  memory->destroy(xsend);
  memory->destroy(xrecv);
  memory->destroy(fsend);
  memory->destroy(frecv);
  memory->destroy(tagsend);
  memory->destroy(tagrecv);

  memory->destroy(xsendall);
  memory->destroy(xrecvall);
  memory->destroy(fsendall);
  memory->destroy(frecvall);
  memory->destroy(tagsendall);
  memory->destroy(tagrecvall);

  memory->destroy(counts);
  memory->destroy(displacements);

  if (NEBLongRange) 
    memory->destroy(nlenall);

}

/* ---------------------------------------------------------------------- */

int FixNEB::setmask()
{
  int mask = 0;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixNEB::init()
{
  int icompute = modify->find_compute(id_pe);
  if (icompute < 0)
    error->all(FLERR,"Potential energy ID for fix neb does not exist");
  pe = modify->compute[icompute];

  // turn off climbing mode, NEB command turns it on after init()

  rclimber = -1;


  // nebatoms = # of atoms in fix group = atoms with inter-replica forces

  bigint count = group->count(igroup);
  if (count > MAXSMALLINT) error->all(FLERR,"Too many active NEB atoms");
  nebatoms = count;

  // comm mode for inter-replica exchange of coords

  if (nreplica == nprocs_universe &&
      nebatoms == atom->natoms && atom->sortfreq == 0) 
    cmode = SINGLE_PROC_DIRECT;
  else if (nreplica == nprocs_universe) cmode = SINGLE_PROC_MAP;
  else cmode = MULTI_PROC;

  // ntotal = total # of atoms in system, NEB atoms or not

  if (atom->natoms > MAXSMALLINT) error->all(FLERR,"Too many atoms for NEB");
  ntotal = atom->natoms;

  if (atom->nlocal > maxlocal) reallocate();

  if (MULTI_PROC && counts == NULL) {
    memory->create(xsendall,ntotal,3,"neb:xsendall");
    memory->create(xrecvall,ntotal,3,"neb:xrecvall");
    memory->create(fsendall,ntotal,3,"neb:fsendall");
    memory->create(frecvall,ntotal,3,"neb:frecvall");
    memory->create(tagsendall,ntotal,"neb:tagsendall");
    memory->create(tagrecvall,ntotal,"neb:tagrecvall");
    memory->create(counts,nprocs,"neb:counts");
    memory->create(displacements,nprocs,"neb:displacements");
  }





}

/* ---------------------------------------------------------------------- */

void FixNEB::min_setup(int vflag)
{
  min_post_force(vflag);

  // trigger potential energy computation on next timestep

  pe->addstep(update->ntimestep+1);
}

/* ---------------------------------------------------------------------- */

void FixNEB::min_post_force(int vflag)
{
  double vprev,vnext,vmax,vmin;
  double delxp,delyp,delzp,delxn,delyn,delzn;
  double delta1[3],delta2[3];
  MPI_Status status;
  MPI_Request request;
  double vIni =0.0;
  // veng = PE of this replica
  // vprev,vnext = PEs of adjacent replicas
  // only proc 0 in each replica communicates

  vprev=vnext=veng = pe->compute_scalar();

  if (ireplica < nreplica-1 && me ==0)
    MPI_Send(&veng,1,MPI_DOUBLE,procnext,0,uworld);
  if (ireplica > 0 && me ==0) 
    MPI_Recv(&vprev,1,MPI_DOUBLE,procprev,0,uworld,&status);

  if (ireplica > 0 && me == 0)
    MPI_Send(&veng,1,MPI_DOUBLE,procprev,0,uworld);
  if (ireplica < nreplica-1 && me == 0)
    MPI_Recv(&vnext,1,MPI_DOUBLE,procnext,0,uworld,&status);

  if (cmode == MULTI_PROC) {
    MPI_Bcast(&vprev,1,MPI_DOUBLE,0,world);
    MPI_Bcast(&vnext,1,MPI_DOUBLE,0,world);
  }

  if (FreeEndFinal)
    {
      if (update->ntimestep==0)
	{EFinalIni = veng;}
    }

  
  if (ireplica==0)
    vIni=veng;
  
  if (FreeEndFinalWithRespToEIni )    {
    if ( me ==0){
      int procFirst;
      procFirst=universe->root_proc[0];
      MPI_Bcast(&vIni,1,MPI_DOUBLE,procFirst,uworld);	  //MPI_Recv(&vIni,1,MPI_DOUBLE,procFirst,0,uworld,&status);
    }
    if (cmode == MULTI_PROC) {
      MPI_Bcast(&vIni,1,MPI_DOUBLE,0,world);
    }
  }
    if (FreeEndIni && ireplica==0 )
    {
    if (me == 0 )
      if (update->ntimestep==0) 
	{
	  EIniIni = veng;
	  if (cmode == MULTI_PROC) 
	    MPI_Bcast(&EIniIni,1,MPI_DOUBLE,0,world);
	}
    }

  // communicate atoms to/from adjacent replicas to fill xprev,xnext
  
  inter_replica_comm();

  // trigger potential energy computation on next timestep

  pe->addstep(update->ntimestep+1);

  /*  ALL THIS SHOULD BE DONE IN inter_replica_comm()
  // xprev,xnext = atom coords of adjacent replicas
  // assume order of atoms in all replicas is the same
  // check that number of atoms hasn't changed

  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double dot = 0.0;
  double prefactor;

  if (nlocal != nebatoms) error->one(FLERR,"Atom count changed in fix neb");

  if (ireplica > 0)
    MPI_Irecv(xprev[0],3*nlocal,MPI_DOUBLE,procprev,0,uworld,&request);
  if (ireplica < nreplica-1)
    MPI_Send(x[0],3*nlocal,MPI_DOUBLE,procnext,0,uworld);
  if (ireplica > 0) MPI_Wait(&request,&status);




  if (ireplica < nreplica-1)
    MPI_Irecv(xnext[0],3*nlocal,MPI_DOUBLE,procnext,0,uworld,&request);
  if (ireplica > 0)
    MPI_Send(x[0],3*nlocal,MPI_DOUBLE,procprev,0,uworld);
  if (ireplica < nreplica-1) MPI_Wait(&request,&status);
  */





  double **x = atom->x;
  int *mask = atom->mask;
  double dot = 0.0;
  double prefactor;

  double **f = atom->f;
  int nlocal = atom->nlocal;

  /* SHOULD BE DONE IN  inter_replica_comm()
  if (ireplica < nreplica-1)
    MPI_Irecv(fnext[0],3*nlocal,MPI_DOUBLE,procnext,0,uworld,&request);
  if (ireplica > 0)
    MPI_Send(f[0],3*nlocal,MPI_DOUBLE,procprev,0,uworld);
  if (ireplica < nreplica-1) MPI_Wait(&request,&status);
  */


  //calculating separation between images
  plen = 0.0;
  nlen = 0.0;
  double tlen = 0.0;
  double  gradnextlen = 0.0;
  double dotFreeEndIniOld=0.0;
  double dotFreeEndFinalOld=0.0;

  dotgrad = 0.0;

  gradlen = 0.0;


  dotpath = 0.0;
  dottangrad = 0.0;

  


  if (ireplica ==nreplica-1){

    
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {

	delxp = x[i][0] - xprev[i][0];
	delyp = x[i][1] - xprev[i][1];
	delzp = x[i][2] - xprev[i][2];
	domain->minimum_image(delxp,delyp,delzp);
	plen += delxp*delxp + delyp*delyp + delzp*delzp;
	dottangrad += delxp* f[i][0]+ delyp*f[i][1]+delzp*f[i][2];
	gradlen += f[i][0]*f[i][0] + f[i][1]*f[i][1] + f[i][2]*f[i][2];
	    
	if (FreeEndFinal||FreeEndFinalWithRespToEIni){
	  tangent[i][0]=delxp;
	  tangent[i][1]=delyp;
	  tangent[i][2]=delzp;
	  tlen += tangent[i][0]*tangent[i][0] + tangent[i][1]*tangent[i][1] +
	    tangent[i][2]*tangent[i][2];
	  dot += f[i][0]*tangent[i][0] + f[i][1]*tangent[i][1] + f[i][2]*tangent[i][2];
	}
      }
    }

  
  else if (ireplica == 0){
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	delxn = xnext[i][0] - x[i][0];
	delyn = xnext[i][1] - x[i][1];
	delzn = xnext[i][2] - x[i][2];
	domain->minimum_image(delxn,delyn,delzn);
	nlen += delxn*delxn + delyn*delyn + delzn*delzn;
	gradnextlen += fnext[i][0]*fnext[i][0] + fnext[i][1]*fnext[i][1] +fnext[i][2] * fnext[i][2];
	dotgrad += f[i][0]*fnext[i][0] + f[i][1]*fnext[i][1] +
	  f[i][2]*fnext[i][2];
	dottangrad += delxn* f[i][0]+ delyn*f[i][1]+delzn*f[i][2];
	gradlen += f[i][0]*f[i][0] + f[i][1]*f[i][1] + f[i][2]*f[i][2];
	if (FreeEndIni)
	  {
	    tangent[i][0]=delxn;
	    tangent[i][1]=delyn;
	    tangent[i][2]=delzn;
	    tlen += tangent[i][0]*tangent[i][0] + tangent[i][1]*tangent[i][1] +
	      tangent[i][2]*tangent[i][2];
	    dot += f[i][0]*tangent[i][0] + f[i][1]*tangent[i][1] + f[i][2]*tangent[i][2];
	  }
      }
  }
  else //not the first or last replica
    {
      vmax = MAX(fabs(vnext-veng),fabs(vprev-veng));
      vmin = MIN(fabs(vnext-veng),fabs(vprev-veng));


      for (int i = 0; i < nlocal; i++)
	if (mask[i] & groupbit) {
	  delxp = x[i][0] - xprev[i][0];
	  delyp = x[i][1] - xprev[i][1];
	  delzp = x[i][2] - xprev[i][2];
	  domain->minimum_image(delxp,delyp,delzp);
	  plen += delxp*delxp + delyp*delyp + delzp*delzp;

	  delxn = xnext[i][0] - x[i][0];
	  delyn = xnext[i][1] - x[i][1];
	  delzn = xnext[i][2] - x[i][2];
	  domain->minimum_image(delxn,delyn,delzn);	  domain->minimum_image(delxn,delyn,delzn);

	  if (vnext > veng && veng > vprev) {	    
	    tangent[i][0]=delxn;
	    tangent[i][1]=delyn;
	    tangent[i][2]=delzn;
	  }
	  else if (vnext < veng && veng < vprev) {
	    tangent[i][0]=delxp;
	    tangent[i][1]=delyp;
	    tangent[i][2]=delzp;
	  }
	  else {
	    if (vnext > vprev) {
	      tangent[i][0] = vmax*delxn + vmin*delxp;
	      tangent[i][1] = vmax*delyn + vmin*delyp;
	      tangent[i][2] = vmax*delzn + vmin*delzp;
	    } else {
	      tangent[i][0] = vmin*delxn + vmax*delxp;
	      tangent[i][1] = vmin*delyn + vmax*delyp;
	      tangent[i][2] = vmin*delzn + vmax*delzp;
	    }

	  }

	  nlen += delxn*delxn + delyn*delyn + delzn*delzn;
	  tlen += tangent[i][0]*tangent[i][0] + tangent[i][1]*tangent[i][1] +
	    tangent[i][2]*tangent[i][2];
	  gradlen += f[i][0]*f[i][0] + f[i][1]*f[i][1] + f[i][2]*f[i][2];
	  dotpath += delxp*delxn + delyp*delyn + delzp*delzn;
	  dottangrad += tangent[i][0]* f[i][0]+ tangent[i][1]*f[i][1]+tangent[i][2]*f[i][2];



	  gradnextlen += fnext[i][0]*fnext[i][0] + fnext[i][1]*fnext[i][1] +fnext[i][2] * fnext[i][2];
	  dotgrad += f[i][0]*fnext[i][0] + f[i][1]*fnext[i][1] +
	    f[i][2]*fnext[i][2];

	  
	  springF[i][0]=kspringPerp*(delxn-delxp);
	  springF[i][1]=kspringPerp*(delyn-delyp);
	  springF[i][2]=kspringPerp*(delzn-delzp);
	 
	}
    }

  double lenall;
  MPI_Allreduce(&nlen,&lenall,1,MPI_DOUBLE,MPI_SUM,world);
  nlen = sqrt(lenall);

  MPI_Allreduce(&plen,&lenall,1,MPI_DOUBLE,MPI_SUM,world);
  plen = sqrt(lenall);

  MPI_Allreduce(&tlen,&lenall,1,MPI_DOUBLE,MPI_SUM,world);
  tlen = sqrt(lenall);

  MPI_Allreduce(&gradlen,&lenall,1,MPI_DOUBLE,MPI_SUM,world);
  gradlen = sqrt(lenall);

  MPI_Allreduce(&gradnextlen,&lenall,1,MPI_DOUBLE,MPI_SUM,world);
  gradnextlen = sqrt(lenall);


  

  double dotpathall;
  double dottangradall;
  double dotgradall;

  MPI_Allreduce(&dotpath,&dotpathall,1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&dottangrad,&dottangradall,1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&dotgrad,&dotgradall,1,MPI_DOUBLE,MPI_SUM,world);


  dotpath=dotpathall;
  dottangrad=dottangradall;
  dotgrad=dotgradall;

  // normalize tangent vector

  if (tlen > 0.0) {
    double tleninv = 1.0/tlen;
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        tangent[i][0] *= tleninv;
        tangent[i][1] *= tleninv;
        tangent[i][2] *= tleninv;
      }
  }


  // first or last replica has no change to forces, just return

  if(ireplica >0 && ireplica<nreplica-1) {
    dottangrad = dottangrad/(tlen*gradlen);

  }
  if(ireplica==0)     dottangrad = dottangrad/(nlen*gradlen);
  if(ireplica==nreplica-1)     dottangrad = dottangrad/(plen*gradlen);
  if(ireplica < nreplica-1)
    {  
      dotgrad = dotgrad /(gradlen*gradnextlen);
    }
    


  if(FreeEndIni&&ireplica == 0)
    {
      if (tlen > 0.0) {
	double dotall;
	MPI_Allreduce(&dot,&dotall,1,MPI_DOUBLE,MPI_SUM,world);
	dot=dotall;

	double tleninv = 1.0/tlen;
	dot *= tleninv;
	if (dot<0)
	  prefactor = -dot - (veng-EIniIni);
	else  prefactor = -dot + (veng-EIniIni);
	for (int i = 0; i < nlocal; i++)
	  if (mask[i] & groupbit) {
	    f[i][0] += prefactor *tangent[i][0];
	    f[i][1] += prefactor *tangent[i][1]; 
	    f[i][2] += prefactor *tangent[i][2];
	  }

      }
    

    }


 
    


  if(FreeEndFinal&&ireplica == nreplica -1)
    {if (tlen > 0.0) {
	double dotall;
	MPI_Allreduce(&dot,&dotall,1,MPI_DOUBLE,MPI_SUM,world);
	dot=dotall;
	double tleninv = 1.0/tlen;
	dot *= tleninv;
	if (dot<0)
	  prefactor = -dot - (veng-EFinalIni);
	else  prefactor = -dot + (veng-EFinalIni);
	for (int i = 0; i < nlocal; i++)
	  if (mask[i] & groupbit) {
	    f[i][0] += prefactor *tangent[i][0];
	    f[i][1] += prefactor *tangent[i][1]; 
	    f[i][2] += prefactor *tangent[i][2];
	  }

      }
    }

  if(FreeEndFinalWithRespToEIni&&ireplica == nreplica -1)
    {if (tlen > 0.0) {
	  double dotall;
	  MPI_Allreduce(&dot,&dotall,1,MPI_DOUBLE,MPI_SUM,world);
	  dot=dotall;
	  double tleninv = 1.0/tlen;
	dot *= tleninv;
	if (dot<0)
	  prefactor = -dot - (veng-vIni);
	else  prefactor = -dot + (veng-vIni);
	for (int i = 0; i < nlocal; i++)
	  if (mask[i] & groupbit) {
	    f[i][0] += prefactor *tangent[i][0];
	    f[i][1] += prefactor *tangent[i][1]; 
	    f[i][2] += prefactor *tangent[i][2];
	  }

      }
    }
   
  double lentot = 0;

  double meanDist,idealPos,lenuntilIm,lenuntilClimber;
  lenuntilClimber=0;
  if(NEBLongRange)
    {
      if (cmode == SINGLE_PROC_DIRECT or cmode == SINGLE_PROC_MAP) 
	{MPI_Allgather(&nlen,1,MPI_DOUBLE,&nlenall[0],1,MPI_DOUBLE,uworld);}
      else
	{
	  /*	  int procRootiIm;
	  double nlentmp;

	  for (int iIm = 0; i < nreplica; i++)
	    {
	      procRootiIm=universe->root_proc[iIm];
	      if (ireplica == iIm && me ==0)
		{ nlentmp=nlen;
		  MPI_Bcast(&nlentmp,1,MPI_DOUBLE,procRootiIm,uworld);
		}
	      else
		{
		  MPI_Bcast(&nlentmp,1,MPI_DOUBLE,procRootiIm,uworld);
		}
	      nlenall[iIm]=nlen;
	    }
      */
	  if (me == 0)
	    MPI_Allgather(&nlen,1,MPI_DOUBLE,&nlenall[0],1,MPI_DOUBLE,rootworld);

	  MPI_Bcast(nlenall,nreplica,MPI_DOUBLE,0,world);
	    
	}



      lenuntilIm = 0;
      for (int i = 0; i < ireplica; i++)
	lenuntilIm += nlenall[i];
      

      for (int i = 0; i < nreplica; i++)
	lentot += nlenall[i];
	
      meanDist = lentot/(nreplica -1);

      if (rclimber>0)
	{
	  for (int i = 0; i < rclimber; i++)
	    lenuntilClimber += nlenall[i];
	  
	  double meanDistBeforeClimber = lenuntilClimber/rclimber;
	  double meanDistAfterClimber = (lentot-lenuntilClimber)/(nreplica-rclimber-1);
      
	  
	  if (ireplica<rclimber)
	    idealPos = ireplica * meanDistBeforeClimber;
	  else
	    idealPos = lenuntilClimber+ (ireplica-rclimber)*meanDistAfterClimber;
	}
      else{
	
	
	idealPos = ireplica * meanDist;
      }


    }


  if (ireplica == 0 || ireplica == nreplica-1)     return ;

  double AngularContr;
  double pi;
  double thetapath;
  pi = 4.0*atan(1.0);
dotpath = dotpath/(plen*nlen);

AngularContr = 0.5 *(1+cos(pi * dotpath));
    
  



  double dotSpringTangent;
  dotSpringTangent=0;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit)
      {
	dot += f[i][0]*tangent[i][0] + f[i][1]*tangent[i][1] + f[i][2]*tangent[i][2];
	dotSpringTangent+=springF[i][0]*tangent[i][0]+springF[i][1]*tangent[i][1]+springF[i][2]*tangent[i][2];}
  }


  double dotSpringTangentall;
  MPI_Allreduce(&dotSpringTangent,&dotSpringTangentall,1,MPI_DOUBLE,MPI_SUM,world);
  dotSpringTangent=dotSpringTangentall;
  double dotall;
  MPI_Allreduce(&dot,&dotall,1,MPI_DOUBLE,MPI_SUM,world);
  dot=dotall;


 
  //  double prefactor, prefSpring;
  double ToDisp;
  if (ireplica == rclimber) {
    prefactor = -2.0*dot;
 
  }
  else {
    
    if(NEBLongRange)
      {prefactor = -dot - kspring*(lenuntilIm-idealPos)/(2*meanDist);}
    else if (StandardNEB)
      {prefactor = -dot + kspring*(nlen-plen);}

    if (FinalAndInterWithRespToEIni&& veng<vIni)
      {
	
	for (int i = 0; i < nlocal; i++)
	  if (mask[i] & groupbit) {
	    f[i][0] = 0;
	    f[i][1] = 0;
	    f[i][2] = 0;
	  }
	prefactor =  kspring*(nlen-plen);
	AngularContr=0;
      }
  }
    
  

  

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {

      f[i][0] += prefactor*tangent[i][0] +AngularContr*(springF[i][0] -dotSpringTangent*tangent[i][0]);
      f[i][1] += prefactor*tangent[i][1]+ AngularContr*(springF[i][1] - dotSpringTangent*tangent[i][1]);
      f[i][2] += prefactor*tangent[i][2]+ AngularContr*(springF[i][2] - dotSpringTangent*tangent[i][2]);



    }

}

/* ----------------------------------------------------------------------
   send/recv NEB atoms to/from adjacent replicas
   received atoms matching my local atoms are stored in xprev,xnext
   replicas 0 and N-1 send but do not receive any atoms
------------------------------------------------------------------------- */


void FixNEB::inter_replica_comm()
{
  int i,m;
  MPI_Request request;
  MPI_Request requests[2];
  MPI_Status statuses[2];

  // reallocate memory if necessary

  if (atom->nlocal > maxlocal) reallocate();

  double **x = atom->x;
  double **f = atom->f;
  tagint *tag = atom->tag;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // -----------------------------------------------------
  // 3 cases: two for single proc per replica
  //          one for multiple procs per replica
  // -----------------------------------------------------

  // single proc per replica
  // all atoms are NEB atoms and no atom sorting
  // direct comm of x -> xprev and x -> xnext

  if (cmode == SINGLE_PROC_DIRECT) {
    if (ireplica > 0)
      MPI_Irecv(xprev[0],3*nlocal,MPI_DOUBLE,procprev,0,uworld,&request);
    if (ireplica < nreplica-1)
      MPI_Send(x[0],3*nlocal,MPI_DOUBLE,procnext,0,uworld);
    if (ireplica > 0) MPI_Wait(&request,MPI_STATUS_IGNORE);
    
    if (ireplica < nreplica-1)
      MPI_Irecv(xnext[0],3*nlocal,MPI_DOUBLE,procnext,0,uworld,&request);
    if (ireplica > 0)
      MPI_Send(x[0],3*nlocal,MPI_DOUBLE,procprev,0,uworld);
    if (ireplica < nreplica-1) MPI_Wait(&request,MPI_STATUS_IGNORE);

    if (ireplica < nreplica-1)
      MPI_Irecv(fnext[0],3*nlocal,MPI_DOUBLE,procnext,0,uworld,&request);
    if (ireplica > 0)
      MPI_Send(f[0],3*nlocal,MPI_DOUBLE,procprev,0,uworld);
    if (ireplica < nreplica-1) MPI_Wait(&request,MPI_STATUS_IGNORE);


    return;
  }

  // single proc per replica
  // but only some atoms are NEB atoms or atom sorting is enabled
  // send atom IDs and coords of only NEB atoms to prev/next proc
  // recv procs use atom->map() to match received coords to owned atoms

  if (cmode == SINGLE_PROC_MAP) {
    m = 0;
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        tagsend[m] = tag[i];
        xsend[m][0] = x[i][0];
        xsend[m][1] = x[i][1];
        xsend[m][2] = x[i][2];
        fsend[m][0] = f[i][0];
        fsend[m][1] = f[i][1];
        fsend[m][2] = f[i][2];
        m++;
      }

    if (ireplica > 0) {
      MPI_Irecv(xrecv[0],3*nebatoms,MPI_DOUBLE,procprev,0,uworld,&requests[0]);
      MPI_Irecv(tagrecv,nebatoms,MPI_LMP_TAGINT,procprev,0,uworld,&requests[1]);
    }
    if (ireplica < nreplica-1) {
      MPI_Send(xsend[0],3*nebatoms,MPI_DOUBLE,procnext,0,uworld);
      MPI_Send(tagsend,nebatoms,MPI_LMP_TAGINT,procnext,0,uworld);
    }

    if (ireplica > 0) {
      MPI_Waitall(2,requests,statuses);
      for (i = 0; i < nebatoms; i++) {
        m = atom->map(tagrecv[i]);
        xprev[m][0] = xrecv[i][0];
        xprev[m][1] = xrecv[i][1];
        xprev[m][2] = xrecv[i][2];
      }
    }
      
    if (ireplica < nreplica-1) {
      MPI_Irecv(xrecv[0],3*nebatoms,MPI_DOUBLE,procnext,0,uworld,&requests[0]);
      MPI_Irecv(frecv[0],3*nebatoms,MPI_DOUBLE,procnext,0,uworld,&requests[0]);
      MPI_Irecv(tagrecv,nebatoms,MPI_LMP_TAGINT,procnext,0,uworld,&requests[1]);
    }
    if (ireplica > 0) {
      MPI_Send(xsend[0],3*nebatoms,MPI_DOUBLE,procprev,0,uworld);
      MPI_Send(fsend[0],3*nebatoms,MPI_DOUBLE,procprev,0,uworld);
      MPI_Send(tagsend,nebatoms,MPI_LMP_TAGINT,procprev,0,uworld);
    }

    if (ireplica < nreplica-1) {
      MPI_Waitall(2,requests,statuses);
      for (i = 0; i < nebatoms; i++) {
        m = atom->map(tagrecv[i]);
        xnext[m][0] = xrecv[i][0];
        xnext[m][1] = xrecv[i][1];
        xnext[m][2] = xrecv[i][2];
        fnext[m][0] = frecv[i][0];
        fnext[m][1] = frecv[i][1];
        fnext[m][2] = frecv[i][2];
      }
    }

    return;
  }

  // multiple procs per replica
  // MPI_Gather all coords and atom IDs to root proc of each replica
  // send to root of adjacent replicas
  // bcast within each replica
  // each proc extracts info for atoms it owns via atom->map()

  m = 0;
  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      tagsend[m] = tag[i];
      xsend[m][0] = x[i][0];
      xsend[m][1] = x[i][1];
      xsend[m][2] = x[i][2];
      fsend[m][0] = f[i][0];
      fsend[m][1] = f[i][1];
      fsend[m][2] = f[i][2];
      m++;
    }

  MPI_Gather(&m,1,MPI_INT,counts,1,MPI_INT,0,world);
  displacements[0] = 0;
  for (i = 0; i < nprocs-1; i++)
    displacements[i+1] = displacements[i] + counts[i];
  MPI_Gatherv(tagsend,m,MPI_LMP_TAGINT,
              tagsendall,counts,displacements,MPI_LMP_TAGINT,0,world);
  for (i = 0; i < nprocs; i++) counts[i] *= 3;
  for (i = 0; i < nprocs-1; i++)
    displacements[i+1] = displacements[i] + counts[i];
  if (xsend){
    MPI_Gatherv(xsend[0],3*m,MPI_DOUBLE,
                xsendall[0],counts,displacements,MPI_DOUBLE,0,world);
    MPI_Gatherv(fsend[0],3*m,MPI_DOUBLE,
                fsendall[0],counts,displacements,MPI_DOUBLE,0,world);
  }
  else {
    MPI_Gatherv(NULL,3*m,MPI_DOUBLE,
                xsendall[0],counts,displacements,MPI_DOUBLE,0,world);
    MPI_Gatherv(NULL,3*m,MPI_DOUBLE,
                fsendall[0],counts,displacements,MPI_DOUBLE,0,world);
  }

  if (ireplica > 0 && me == 0) {
    MPI_Irecv(xrecvall[0],3*nebatoms,MPI_DOUBLE,procprev,0,uworld,&requests[0]);
    MPI_Irecv(tagrecvall,nebatoms,MPI_LMP_TAGINT,procprev,0,uworld,
              &requests[1]);
  }
  if (ireplica < nreplica-1 && me == 0) {
    MPI_Send(xsendall[0],3*nebatoms,MPI_DOUBLE,procnext,0,uworld);
    MPI_Send(tagsendall,nebatoms,MPI_LMP_TAGINT,procnext,0,uworld);
  }

  if (ireplica > 0) {
    if (me == 0) MPI_Waitall(2,requests,statuses);

    MPI_Bcast(tagrecvall,nebatoms,MPI_INT,0,world);
    MPI_Bcast(xrecvall[0],3*nebatoms,MPI_DOUBLE,0,world);

    for (i = 0; i < nebatoms; i++) {
      m = atom->map(tagrecvall[i]);
      if (m < 0 || m >= nlocal) continue;
      xprev[m][0] = xrecvall[i][0];
      xprev[m][1] = xrecvall[i][1];
      xprev[m][2] = xrecvall[i][2];
    }
  }

  if (ireplica < nreplica-1 && me == 0) {
    MPI_Irecv(xrecvall[0],3*nebatoms,MPI_DOUBLE,procnext,0,uworld,&requests[0]);
    MPI_Irecv(frecvall[0],3*nebatoms,MPI_DOUBLE,procnext,0,uworld,&requests[0]);
    MPI_Irecv(tagrecvall,nebatoms,MPI_LMP_TAGINT,procnext,0,uworld,
              &requests[1]);
  }
  if (ireplica > 0 && me == 0) {
    MPI_Send(xsendall[0],3*nebatoms,MPI_DOUBLE,procprev,0,uworld);
    MPI_Send(fsendall[0],3*nebatoms,MPI_DOUBLE,procprev,0,uworld);
    MPI_Send(tagsendall,nebatoms,MPI_LMP_TAGINT,procprev,0,uworld);
  }

  if (ireplica < nreplica-1) {
    if (me == 0) MPI_Waitall(2,requests,statuses);

    MPI_Bcast(tagrecvall,nebatoms,MPI_INT,0,world);
    MPI_Bcast(xrecvall[0],3*nebatoms,MPI_DOUBLE,0,world);
    MPI_Bcast(frecvall[0],3*nebatoms,MPI_DOUBLE,0,world);

    for (i = 0; i < nebatoms; i++) {
      m = atom->map(tagrecvall[i]);
      if (m < 0 || m >= nlocal) continue;
      xnext[m][0] = xrecvall[i][0];
      xnext[m][1] = xrecvall[i][1];
      xnext[m][2] = xrecvall[i][2];
      fnext[m][0] = frecvall[i][0];
      fnext[m][1] = frecvall[i][1];
      fnext[m][2] = frecvall[i][2];
    }
  }
}


/* ----------------------------------------------------------------------
   reallocate xprev,xnext,tangent arrays if necessary
   reallocate communication arrays if necessary
------------------------------------------------------------------------- */

void FixNEB::reallocate()
{

  memory->destroy(xprev);
  memory->destroy(xnext);
  memory->destroy(tangent);
  memory->destroy(fnext);
  memory->destroy(springF);

  if (NEBLongRange)
        {memory->destroy(nlenall);
	  memory->create(nlenall,nreplica,"neb:nlenall");
	}
  


  if (cmode != SINGLE_PROC_DIRECT) {
    memory->destroy(xsend);
    memory->destroy(fsend);
    memory->destroy(xrecv);
    memory->destroy(frecv);
    memory->destroy(tagsend);
    memory->destroy(tagrecv);
  }

  maxlocal = atom->nmax;

  memory->create(xprev,maxlocal,3,"neb:xprev");
  memory->create(xnext,maxlocal,3,"neb:xnext");
  memory->create(tangent,maxlocal,3,"neb:tangent");
  memory->create(fnext,maxlocal,3,"neb:fnext");
  memory->create(springF,maxlocal,3,"neb:springF");


  if (cmode != SINGLE_PROC_DIRECT) {
    memory->create(xsend,maxlocal,3,"neb:xsend");
    memory->create(fsend,maxlocal,3,"neb:fsend");
    memory->create(xrecv,maxlocal,3,"neb:xrecv");
    memory->create(frecv,maxlocal,3,"neb:frecv");
    memory->create(tagsend,maxlocal,"neb:tagsend");
    memory->create(tagrecv,maxlocal,"neb:tagrecv");
  }
}
