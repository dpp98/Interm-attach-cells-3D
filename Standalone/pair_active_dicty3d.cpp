/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "pair_active_dicty3d.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "domain.h"

#include <cmath>
#include <cstring>
#include <vector>
#include <algorithm>
#include <random>
#include <iostream>
#include <ctime>
#include <chrono>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairActiveDicty3d::PairActiveDicty3d(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 1;
  respa_enable = 0;
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairActiveDicty3d::~PairActiveDicty3d()
{
  if (allocated) {
    memory->destroy(cutsq);
    memory->destroy(setflag);
    memory->destroy(cut_near);
    memory->destroy(cut_far);
    memory->destroy(f0);
    memory->destroy(ta);
    memory->destroy(var_ta);
  }
}

/* ---------------------------------------------------------------------- */

void PairActiveDicty3d::compute(int eflag, int vflag)
{
    
    std::srand(std::time(nullptr));
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
   
    std::random_device rd;
    std::mt19937 g(rd());

    int i, j, ii, jj, inum, jnum, itype, jtype, sid, nsid;
    double xtmp, ytmp, ztmp, delx, dely, delz, evdwl, fpair;
    double nxtmp, nytmp, nztmp, ndist, rsq, r2inv, r6inv, r7inv;
    double vxtmp, vytmp, vztmp, nvxtmp, nvytmp, nvztmp;
    double fmag, fmagid, tmp, l0, id_x, id_y, id_z, fmax; 

    std::uniform_real_distribution<double> uni_dis(0.0, 1.0); 
    int *ilist, *jlist, *numneigh, **firstneigh;
    int atsteps, neigh_id;
    int ch_n;

    evdwl = 0.0;
    ev_init(eflag, vflag);
    
    double **x = atom->x;
    double **f = atom->f;
    double **v = atom->v;
    int *type = atom->type;
    int *id = atom->tag;
    int nlocal = atom->nlocal;

    double *tas = atom->q;           // Number of attachment steps
    double *Nid = atom->radius;      // Global ID of neighbor
    //double *bias = atom->rmass;
    double p, pcnew, psnew;
    double pc; 
    int lnid, num_neigh, num_lneigh, num_rneigh, num_bneigh, num_sneigh;
    std::vector<int> neigh;
    std::vector<double> pdf, cdf;
    double totprob;
    
    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
    
    // loop over neighbors of my atoms
        
    for (ii = 0; ii < inum; ii++) {
        
        i = ilist[ii];
        xtmp = x[i][0];
        ytmp = x[i][1];
        ztmp = x[i][2];
	
        vxtmp = v[i][0];
        vytmp = v[i][1];
        vztmp = v[i][2];
 
        itype = type[i];
        jlist = firstneigh[i];
        jnum = numneigh[i];
    
        if (itype != 1) continue;
        
        atsteps = int(tas[i]);
        neigh_id = int(Nid[i]);
        
        neigh.clear(); pdf.clear(); cdf.clear();
        num_rneigh = 0; num_lneigh = 0; num_sneigh = 0; num_bneigh = 0; num_neigh = 0;
        totprob = 0.0;
        
        // check distance with neigh
        // if dist with neigh is > cutoff then atsteps = 0
        if (neigh_id > 0) {
            
            lnid = atom->map(neigh_id);
                
            nxtmp = x[lnid][0];
            nytmp = x[lnid][1];
       
            // correct for periodicity
            if (nxtmp-xtmp > 0.5*(domain->boxhi[0]-domain->boxlo[0])) nxtmp = nxtmp - (domain->boxhi[0] - domain->boxlo[0]);
            if (nxtmp-xtmp < 0.5*(domain->boxlo[0]-domain->boxhi[0])) nxtmp = nxtmp +  (domain->boxhi[0] - domain->boxlo[0]);
            if (nytmp-ytmp > 0.5*(domain->boxhi[1]-domain->boxlo[1])) nytmp = nytmp - (domain->boxhi[1] - domain->boxlo[1]);
            if (nytmp-ytmp < 0.5*(domain->boxlo[1]-domain->boxhi[1])) nytmp = nytmp +  (domain->boxhi[1] - domain->boxlo[1]);

                      
            nztmp = x[lnid][2];
            jtype = type[lnid];
        
            delx = nxtmp - xtmp;
            dely = nytmp - ytmp;
            delz = nztmp - ztmp;
        
            ndist = std::sqrt(delx*delx + dely*dely + delz*delz);
        
            // break the attachment if the distance exceeds maximum limit
            if (ndist > cut_far[itype][jtype]) {
                atsteps = 0;
            }
            
        }
        
        if (neigh_id <= 0) atsteps = 0;
        
        if (atsteps <= 0) { // choose new neighbor
            
            // store the neighs in a vector
            for (jj = 0; jj < jnum; jj++) {
                
                j = jlist[jj];
                j &= NEIGHMASK;
                
                nxtmp = x[j][0];
                nytmp = x[j][1];
                nztmp = x[j][2];
                jtype = type[j];
        
                delx = nxtmp - xtmp;
                dely = nytmp - ytmp;
                delz = nztmp - ztmp;
                
                ndist = std::sqrt(delx*delx + dely*dely + delz*delz);
                                
                if (ndist > cut_near[itype][jtype] && ndist < cut_far[itype][jtype]) {   // insert && jtype == 1 is only for measuring surface tension
                     
                    neigh.push_back(id[j]);
                    pdf.push_back(0.0);
                    cdf.push_back(0.0);
                    
                    if (delx > 0.0) num_rneigh++;
                    else num_lneigh++;
        	
        	    if (jtype == 2) num_sneigh++;
        	    else num_bneigh++;
        	  
                }
            }
	    
            std::shuffle(neigh.begin(), neigh.end(), g); 
            num_neigh = static_cast<int>(neigh.size());
	    
	    for (int k = 0; k < num_neigh; k++) {
               
               lnid = atom->map(neigh[k]);
               
               nxtmp = x[lnid][0];
               nytmp = x[lnid][1];
               
	       // correct for periodicity
	       if (nxtmp-xtmp > 0.5*(domain->boxhi[0]-domain->boxlo[0])) nxtmp = nxtmp - (domain->boxhi[0] - domain->boxlo[0]);
               if (nxtmp-xtmp < 0.5*(domain->boxlo[0]-domain->boxhi[0])) nxtmp = nxtmp +  (domain->boxhi[0] - domain->boxlo[0]);
               if (nytmp-ytmp > 0.5*(domain->boxhi[1]-domain->boxlo[1])) nytmp = nytmp - (domain->boxhi[1] - domain->boxlo[1]);
               if (nytmp-ytmp < 0.5*(domain->boxlo[1]-domain->boxhi[1])) nytmp = nytmp +  (domain->boxhi[1] - domain->boxlo[1]);

               delx = nxtmp - xtmp;
               jtype = type[lnid];
 
            	if (delx >= 0.0 && jtype==2)
            	    pdf[k] = pc;
            	else if (delx < 0.0 && jtype==2)
            	    pdf[k] = 1.0-pc;
            	else if (delx >= 0.0 && jtype==1)
            	    pdf[k] = pc;
            	else 
            	    pdf[k] = 1.0-pc;	
            }
	    
            // normalize the pdf
            totprob = 0.0;
            for (int k = 0; k < num_neigh; k++)
                totprob += pdf[k];
            for (int k = 0; k < num_neigh; k++) 
                pdf[k] = pdf[k] / totprob;
	     
            cdf[0] = pdf[0];
            for (int k = 1; k < num_neigh; k++)
                cdf[k] = cdf[k-1]+pdf[k];
            
            // error check
            if (cdf[num_neigh-1]>1.000000000001)
                std::cout << "Error stop!" << std::endl;
            
            p = uni_dis(g);
            for (int k = 0; k < num_neigh; k++) {
                if (p < cdf[k]) {
                    neigh_id = neigh[k];
                    break; 
                }
            }
            
            std::normal_distribution<double> ndis{ta[itype][jtype], var_ta[itype][jtype]*ta[itype][jtype]};
            
            tas[i] = ndis(g)/timestep;

	    // If you want cell-cell attachments to be longer by a factor of 'n', then uncomment 
	    // the below 3 lines and replace n by the number you want 
	    //lnid = atom->map(neigh_id);
            //jtype = type[lnid];
	    //if (jtype == 1) tas[i] = 100000.0*tas[i];	

            if (tas[i] < 0) tas[i] = 1.0;
            
            Nid[i] = neigh_id;

        }
        
        // apply force due to pseudopod traction
        if (neigh_id > 0) {
        
            lnid = atom->map(neigh_id);

            nxtmp = x[lnid][0];
            nytmp = x[lnid][1];
            nztmp = x[lnid][2];
            jtype = type[lnid];
            
            // correct for periodicity
	    if (nxtmp-xtmp > 0.5*(domain->boxhi[0]-domain->boxlo[0])) nxtmp = nxtmp - (domain->boxhi[0] - domain->boxlo[0]);
            if (nxtmp-xtmp < 0.5*(domain->boxlo[0]-domain->boxhi[0])) nxtmp = nxtmp +  (domain->boxhi[0] - domain->boxlo[0]);
            if (nytmp-ytmp > 0.5*(domain->boxhi[1]-domain->boxlo[1])) nytmp = nytmp - (domain->boxhi[1] - domain->boxlo[1]);
            if (nytmp-ytmp < 0.5*(domain->boxlo[1]-domain->boxhi[1])) nytmp = nytmp +  (domain->boxhi[1] - domain->boxlo[1]);

            delx = nxtmp - xtmp;
            dely = nytmp - ytmp;
            delz = nztmp - ztmp;
                    
            ndist = std::sqrt(delx*delx + dely*dely + delz*delz);

            if (ndist < cut_near[itype][jtype])
                fmag = 0.0;
            else
            fmag = f0[itype][jtype] * (ndist - cut_near[itype][jtype]);

            if (jtype == 2) fmag *= wf;
            
            f[i][0] += (delx/ndist) * fmag;
            f[i][1] += (dely/ndist) * fmag;
            f[i][2] += (delz/ndist) * fmag;
                         
            f[lnid][0] -= (delx/ndist) * fmag;
            f[lnid][1] -= (dely/ndist) * fmag;
            f[lnid][2] -= (delz/ndist) * fmag;
            
        }
               
        ///////////////////////////////
        // reduce the attachment steps
        tas[i] -= 1.0;
        
    }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairActiveDicty3d::allocate()
{
    
  allocated = 1;
  int np1 = atom->ntypes + 1;

  memory->create(setflag, np1, np1, "pair:setflag");
  for (int i = 1; i < np1; i++)
    for (int j = i; j < np1; j++) setflag[i][j] = 0;

  memory->create(cutsq, np1, np1, "pair:cutsq");
  memory->create(cut_far, np1, np1, "pair:cut_far");
  memory->create(cut_near, np1, np1, "pair:cut::cut_near");
  memory->create(f0, np1, np1, "pair:f0");
  memory->create(offset, np1, np1, "pair:offset");
  memory->create(ta,np1,np1,"pair::ta");
  memory->create(var_ta,np1,np1,"pair::var_ta");
  
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairActiveDicty3d::settings(int narg, char **arg)
{
    
  if (narg != 1) error->all(FLERR, "Pair style dicty must have exactly three");
  double cut_far_one = utils::numeric(FLERR, arg[0], false, lmp);
  // reset per-type pair cutoffs that have been explicitly set previously

  if (allocated) {
    for (int i = 1; i <= atom->ntypes; i++)
      for (int j = i; j <= atom->ntypes; j++) {
          if (setflag[i][j]) cut_far[i][j] = cut_far_one; 
      }
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairActiveDicty3d::coeff(int narg, char **arg)
{
        
  if (narg != 11) error->all(FLERR, "Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  double f0_one = utils::numeric(FLERR, arg[2], false, lmp);

  double cut_near_one = utils::numeric(FLERR, arg[3], false, lmp);
  double cut_far_one = utils::numeric(FLERR, arg[4], false, lmp);
  double ta_one = utils::numeric(FLERR, arg[5], false, lmp);
  double var_ta_one = utils::numeric(FLERR, arg[6], false, lmp);
  timestep = utils::numeric(FLERR, arg[7], false, lmp);
  ps = utils::numeric(FLERR, arg[8], false, lmp);
  wf = utils::numeric(FLERR, arg[9], false, lmp);
  pc =  utils::numeric(FLERR, arg[10], false, lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = jlo; j <= jhi; j++) {
      f0[i][j] = f0_one;
      cut_far[i][j] = cut_far_one;
      cut_near[i][j] = cut_near_one;
      ta[i][j] = ta_one;
      var_ta[i][j] = var_ta_one;
      setflag[i][j] = 1;
      count++;
    }
  }
 
  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairActiveDicty3d::init_style()
{
  //if (atom->tag_enable == 0)
  //  error->all(FLERR,"Pair style Tersoff requires atom IDs");
  //if (force->newton_pair == 0)
  //  error->all(FLERR,"Pair style Tersoff requires newton pair on");

  // need a full neighbor list

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}


/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairActiveDicty3d::init_one(int i, int j)
{
        
  if (setflag[i][j] == 0) error->all(FLERR, "All pair coeffs are not set");

  if (offset_flag) {
    offset[i][j] = f0[i][j] * cut_far[i][j];
  } else
    offset[i][j] = 0.0;

  f0[j][i] = f0[i][j];
  offset[j][i] = offset[i][j];
  ta[j][i] = ta[i][j];
  var_ta[j][i] = var_ta[i][j];
  cut_far[j][i] = cut_far[i][j];
  cut_near[j][i] = cut_near[i][j];

  return cut_far[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */


void PairActiveDicty3d::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i, j;
  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j], sizeof(int), 1, fp);
      if (setflag[i][j]) {
        fwrite(&f0[i][j], sizeof(double), 1, fp);
        fwrite(&cut_near[i][j], sizeof(double), 1, fp);
        fwrite(&cut_far[i][j], sizeof(double), 1, fp);
        fwrite(&ta[i][j], sizeof(double), 1, fp);
        fwrite(&var_ta[i][j], sizeof(double), 1, fp);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairActiveDicty3d::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i, j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR, &setflag[i][j], sizeof(int), 1, fp, nullptr, error);
      MPI_Bcast(&setflag[i][j], 1, MPI_INT, 0, world);
      if (setflag[i][j]) {
        if (me == 0) {
          
          utils::sfread(FLERR, &f0[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &cut_near[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &cut_far[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &ta[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &var_ta[i][j], sizeof(double), 1, fp, nullptr, error);
        }
    
        MPI_Bcast(&f0[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&cut_near[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&cut_far[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&ta[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&var_ta[i][j], 1, MPI_DOUBLE, 0, world);
        
      }
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairActiveDicty3d::write_restart_settings(FILE *fp)
{
  fwrite(&timestep, sizeof(double), 1, fp);
  fwrite(&offset_flag, sizeof(int), 1, fp);
  fwrite(&mix_flag, sizeof(int), 1, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairActiveDicty3d::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR, &timestep, sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &offset_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &mix_flag, sizeof(int), 1, fp, nullptr, error);
  }
  MPI_Bcast(&timestep, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&offset_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&mix_flag, 1, MPI_INT, 0, world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairActiveDicty3d::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp, "%d %g\n", i, f0[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairActiveDicty3d::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp, "%d %d %g %g %g %g %g\n", i, j, ta[i][j], var_ta[i][j], f0[i][j], cut_near[i][j], cut_far[i][j]);
}

/* ---------------------------------------------------------------------- */

//double PairActiveDicty::single(int /*i*/, int /*j*/, int itype, int jtype, double rsq,
//                             double /*factor_coul*/, double factor_lj, double &fforce)
//{
//  double r, dr, aexp, bexp;

//  r = sqrt(rsq);
  
//  fforce = f0[itype][jtype];
//  return f0[itype][jtype];
//}


/* ---------------------------------------------------------------------- */

void *PairActiveDicty3d::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str, "f0") == 0) return (void *) f0;
  return nullptr;
}

