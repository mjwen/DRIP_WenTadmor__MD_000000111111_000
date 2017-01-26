/*
*
* CDDL HEADER START
*
* The contents of this file are subject to the terms of the Common Development
* and Distribution License Version 1.0 (the "License").
*
* You can obtain a copy of the license at
* http://www.opensource.org/licenses/CDDL-1.0.  See the License for the
* specific language governing permissions and limitations under the License.
*
* When distributing Covered Code, include this CDDL HEADER in each file and
* include the License file in a prominent location with the name LICENSE.CDDL.
* If applicable, add the following below this CDDL HEADER, with the fields
* enclosed by brackets "[]" replaced with your own identifying information:
*
* Portions Copyright (c) [yyyy] [name of copyright owner]. All rights reserved.
*
* CDDL HEADER END
*
*
* Copyright (c) 2017, Regents of the University of Minnesota.  All rights reserved.
*
* Contributors:
*    Ryan S. Elliott
*    Ellad B. Tadmor
*    Valeriu Smirichinski
*    Mingjian Wen
*/

/*******************************************************************************
*
*  Registry-dependent interlayer potential of Kolmogorov and Crespi Model Driver
*
*  Language: C
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "KIM_API_C.h"
#include "KIM_API_status.h"

/******************************************************************************
 * Below are the definitions and values of all Model parameters
 *******************************************************************************/
#define DIM 3       /* dimensionality of space */
#define SPEC1 1     /* internal species code */
#define SPEC2 2     /* internal species code */


/* Define prototypes for Model Driver init */
int model_driver_init(void* km, char* paramfile_names, int* nmstrlen, int* numparamfiles);

/* Define prototypes for Model (Driver) reinit, compute, and destroy */
/* defined as static to avoid namespace clashes with other Models    */
static int reinit(void* km);
static int destroy(void* km);
static int compute(void* km);

/* local function prototypes */
static void calc_phi_repulsive(const double C0, const double C2, const double C4,
    const double C, const double delta, const double lambda, const double z0,
    const double r, const double rij[DIM], const double normali[DIM],
    const double normalj[DIM], double *const phi);

static void calc_phi_attractive(const double A, const double z0, const double r,
    double *const phi);

static void get_normal(const double Rvec1[DIM], const double Rvec2[DIM],
    const double Rvec3[DIM], double normal[DIM]);

static double get_rho_sq(const double r, const double rvec[DIM], const double n[DIM]);

static double get_td(const double C0, const double C2, const double C4,
    const double rho_sq, const double delta);

static void cross(const double x[DIM], const double y[DIM], double rslt[DIM]);

/* Define model_buffer structure */
struct model_buffer {
  int NBC;
  int IterOrLoca;
  int energy_ind;
  int forces_ind;
  int particleEnergy_ind;
  int process_dEdr_ind;
  int process_d2Edr2_ind;
  int model_index_shift;
  int numberOfParticles_ind;
  int numberOfSpecies_ind;
  int particleSpecies_ind;
  int particleStatus_ind;
  int coordinates_ind;
  int boxSideLengths_ind;
  int get_neigh_ind;
  int cutoff_ind;

  double* cutoff;
  double* cutsq;
  double* C0;
  double* C2;
  double* C4;
  double* C;
  double* delta;
  double* lambda;
  double* A;
  double* z0;
};


/* the repulsive part due to collapse of electron clouds */
static void calc_phi_repulsive(const double C0, const double C2, const double C4,
    const double C, const double delta, const double lambda, const double z0,
    const double r, const double rij[DIM],const double ni[DIM], const double nj[DIM],
    double *const phi)
{
  double rhoij_sq, rhoji_sq;
  double tdij, tdji;

  rhoij_sq = get_rho_sq(r, rij, ni);
  rhoji_sq = get_rho_sq(r, rij, nj);

  tdij = get_td(C0, C2, C4, rhoij_sq, delta);
  tdji = get_td(C0, C2, C4, rhoji_sq, delta);

  *phi = exp(-lambda*(r-z0)) * (C+tdij+tdji);
}

/* the r^{-6} attractive part */
static void calc_phi_attractive(const double A, const double z0, const double r,
    double *const phi)
{
  double roz0_sq = r*r/(z0*z0);
  *phi = -A/(roz0_sq*roz0_sq*roz0_sq);
}

/* compute normal */
/* the direction of the normal does not matter */
static void get_normal(const double Rvec1[DIM], const double Rvec2[DIM],
    const double Rvec3[DIM], double normal[DIM])
{
  /* local variables */
  double mag_normal;
  int i;

  cross(Rvec1, Rvec2, normal);

  mag_normal=0.0;
  for (i=0; i<DIM; i++) {
    mag_normal += normal[i]*normal[i];
  }
  mag_normal = sqrt(mag_normal);

  for (i=0; i<DIM; i++) {
    normal[i] /= mag_normal;
  }

}

/* compute square of rhoij */
static double get_rho_sq(const double r, const double rvec[DIM], const double n[DIM])
{
  double rhosq;
  double n_dot_r;
  int i;

  n_dot_r = 0.0;
  for (i=0; i<DIM; i++) {
    n_dot_r += n[i]*rvec[i];
  }
  rhosq = r*r - n_dot_r*n_dot_r;

  return rhosq;
}

/* transverse decay function f(rho)*/
static double get_td(const double C0, const double C2, const double C4,
    const double rho_sq, const double delta)
{
  double rod_sq;
  double rslt;

  rod_sq = rho_sq/(delta*delta);
  rslt = exp(-rod_sq) * (C0 + (rod_sq*(C2 + rod_sq*C4)));

  return rslt;
}


/* cross product, x cross y */
static void cross(const double x[DIM], const double y[DIM], double rslt[DIM])
{
  rslt[0] = x[1]*y[2] - x[2]*y[1];
  rslt[1] = x[2]*y[0] - x[0]*y[2];
  rslt[2] = x[0]*y[1] - x[1]*y[0];
}




/* compute function */
static int compute(void* km)
{
  /* local variables */
  intptr_t* pkim = *((intptr_t**) km);
  double R;
  double Rsqij;
  double phi_repul;
  double phi_attract;
/* double dphi_repul;
  double dphi_attract;
  double dEidr_attrat;
  double dEidr_repul;
  double dEidr_attract;
*/
  double Rij[DIM];

/*  double *pRij = &(Rij[0]);
*/
  int ier;
  int i;
  int j;
  int k;
  int ii;
  int jj;
  int kdim;
  int currentAtom;
  int* neighListOfCurrentAtom;
  struct model_buffer* buffer;
  int comp_energy;
  int comp_force;
  int comp_particleEnergy;
  int comp_process_dEdr;
  int comp_process_d2Edr2;
  int NBC;
  int IterOrLoca;
  int model_index_shift;
  int zero = 0;
  int one = 1;
  int request;

  int* nAtoms;
  int* nSpecies;
  int* particleSpecies;
  int* particleStatus;
  double cutoff;
  double cutsq;
  double* Rij_list;
  double* coords;
  double* energy;
  double* force;
  double* particleEnergy;
  double* boxSideLengths;
  int numOfAtomNeigh;
  int iSpecies;
  int jSpecies;
  int inter_idx;
  typedef int (*get_neigh_ptr)(void *,int *,int *,int *, int *, int **, double **);
  get_neigh_ptr get_neigh = NULL;

  /* params */
  double C0;
  double C2;
  double C4;
  double C;
  double delta;
  double lambda;
  double A;
  double z0;

  /* layers and normal */
  int* in_layer; /* atoms in which layer , -1: not classified in any layer */
  int* layer; /* layer contains which atoms */
  double (*normal)[DIM];
  int currentLayer;
  int nin; /* number of atoms in current layer */
  int nremain; /* number of atoms not included in any layer */
  int nlayers;
  double tmp_normal[DIM];
  double neigh1_Rsq, neigh2_Rsq, neigh3_Rsq;
  double neigh1_Rvec[DIM], neigh2_Rvec[DIM], neigh3_Rvec[DIM];
  int ilayer;
  int jlayer;
  double cutsq_layer = (0.72*3.44)*(0.72*3.44);





  /* get buffer from KIM object */
  buffer = (struct model_buffer*) KIM_API_get_model_buffer(pkim, &ier);
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_model_buffer", ier);
    return ier;
  }

  /* unpack info from the buffer */
  NBC = buffer->NBC;
  IterOrLoca = buffer->IterOrLoca;
  model_index_shift = buffer->model_index_shift;

  /* check to see if we have been asked to compute the forces, particleEnergy, and dEdr */
  KIM_API_getm_compute_by_index(pkim, &ier, 5*3,
      buffer->energy_ind,         &comp_energy,         1,
      buffer->forces_ind,         &comp_force,          1,
      buffer->particleEnergy_ind, &comp_particleEnergy, 1,
      buffer->process_dEdr_ind,   &comp_process_dEdr,   1,
      buffer->process_d2Edr2_ind, &comp_process_d2Edr2, 1);
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_compute_by_index", ier);
    return ier;
  }

  KIM_API_getm_data_by_index(pkim, &ier, 9*3,
      buffer->numberOfParticles_ind,           &nAtoms,         1,
      buffer->numberOfSpecies_ind,             &nSpecies,       1,
      buffer->particleSpecies_ind,             &particleSpecies,1,
      buffer->particleStatus_ind,              &particleStatus, 1,
      buffer->coordinates_ind,                 &coords,         1,
      buffer->boxSideLengths_ind,              &boxSideLengths, (NBC==2),
      buffer->energy_ind,                      &energy,         comp_energy,
      buffer->forces_ind,                      &force,          comp_force,
      buffer->particleEnergy_ind,              &particleEnergy, comp_particleEnergy);
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_data_by_index", ier);
    return ier;
  }
  if (NBC!=3) {
    get_neigh = (get_neigh_ptr) KIM_API_get_method_by_index(pkim, buffer->get_neigh_ind, &ier);
    if (KIM_STATUS_OK > ier) {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_method_by_index", ier);
      return ier;
    }
  }

  /* initialize potential energies, forces, and virial term */
  if (comp_particleEnergy) {
    for (i = 0; i < *nAtoms; ++i) {
      particleEnergy[i] = 0.0;
    }
  }

  if (comp_energy) {
    *energy = 0.0;
  }

  if (comp_force) {
    for (i = 0; i < *nAtoms; ++i) {
      for (kdim = 0; kdim < DIM; ++kdim) {
        force[i*DIM + kdim] = 0.0;
      }
    }
  }

  /* Initialize neighbor handling for CLUSTER NBC */
  if (3 == NBC) /* CLUSTER */ {
    neighListOfCurrentAtom = (int *) malloc((*nAtoms)*sizeof(int));
  }

  /* Initialize neighbor handling for Iterator mode */
  if (1 == IterOrLoca) {
    ier = (*get_neigh)(&pkim, &zero, &zero, &currentAtom, &numOfAtomNeigh,
        &neighListOfCurrentAtom, &Rij_list);
    /* check for successful initialization */
    if (KIM_STATUS_NEIGH_ITER_INIT_OK != ier) {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_neigh", ier);
      ier = KIM_STATUS_FAIL;
      return ier;
    }
  }


  /*********************************************/
  /* create layers and find normal at each atom */

  /* allocate memory */
  layer = (int*) malloc((*nAtoms)*sizeof(int));
  in_layer = (int*) malloc((*nAtoms)*sizeof(int));
  normal = (double(*)[DIM]) malloc((*nAtoms)*DIM*sizeof(double));

  /* init atoms layer to -1 (not in any layer)*/
  for (k=0; k<*nAtoms; k++) {
    in_layer[k] = -1;
  }


  nlayers = 1;
  nremain = *nAtoms;

  while(1) {  /* to create all layers */

    /* init all atoms in current layer have atom number -1 */
    for (k=0; k<*nAtoms; k++) {
      layer[k] = -1;
    }

    /* find an atom not incldued in any layer and start with it */
    currentLayer = nlayers - 1;
    for (k=0; k<*nAtoms; k++) {
      if (in_layer[k] == -1) {
        in_layer[k] = currentLayer;
        layer[0] = k;
        break;
      }
    }

    nin = 1;
    ii = 0;
    while(1) { /* to find all atoms in currentLayer */

      i = layer[ii];

      /* Set up neighbor list of atom i all NBC methods */

      if (1 == IterOrLoca) { /* ITERATOR mode */
        /*NOTE we should remove iterator mode from descriptor file */
        ier = KIM_STATUS_FAIL;
        KIM_API_report_error(__LINE__, __FILE__, "Iterator mode does not supported", ier);
        return ier;
      } else {
        if (3 == NBC) {/* CLUSTER NBC method */
          numOfAtomNeigh = *nAtoms - 1;
          for (kdim = 0; kdim < *nAtoms; ++kdim) {
            if (kdim < i)
              neighListOfCurrentAtom[kdim] = kdim - model_index_shift;
            if (kdim > i)
              neighListOfCurrentAtom[kdim - 1] = kdim - model_index_shift;
          }
          ier = KIM_STATUS_OK;
        } else {
          request = i - model_index_shift;
          ier = (*get_neigh)(&pkim, &one, &request, &currentAtom, &numOfAtomNeigh,
              &neighListOfCurrentAtom, &Rij_list);
          if (KIM_STATUS_OK != ier) {
            KIM_API_report_error(__LINE__, __FILE__, "get_neigh", ier);
            ier = KIM_STATUS_FAIL;
            return ier;
          }
        }
      }


/* NOTE delete */
/*      printf("i=%d,numneigh=%d\n",i,numOfAtomNeigh);
      if(i==702) {
        printf("0 2 %f %f %f 0 0 0\n", coords[i*DIM+0],  coords[i*DIM+1], coords[i*DIM+2]);
      }
*/


      /* init neigh1 to be the 1st nearest neigh, neigh3 the 3rd nearest */
      neigh1_Rsq = 1.0e10;
      neigh2_Rsq = 2.0e10;
      neigh3_Rsq = 3.0e10;
      for (k=0; k<DIM; k++) {
        neigh1_Rvec[k] = 0.0;
        neigh2_Rvec[k] = 0.0;
        neigh3_Rvec[k] = 0.0;
      }

      /* loop over the neighbors of atom i */
      for (jj = 0; jj < numOfAtomNeigh; ++ jj) {

        j = neighListOfCurrentAtom[jj] + model_index_shift; /* get neighbor ID */

        /* compute relative position vector and squared distance */
        Rsqij = 0.0;
        for (kdim = 0; kdim < DIM; ++kdim) {
          if (0 != NBC) { /* all methods except NEIGH_RVEC */
            Rij[kdim] = coords[j*DIM + kdim] - coords[i*DIM + kdim];
          } else {         /* NEIGH_RVEC method */
            Rij[kdim] = Rij_list[jj*DIM + kdim];
          }

          /* apply periodic boundary conditions if required */
          if (2 == NBC) {
            if (fabs(Rij[kdim]) > 0.5*boxSideLengths[kdim]) {
              Rij[kdim] -= (Rij[kdim]/fabs(Rij[kdim]))*boxSideLengths[kdim];
            }
          }

          /* compute squared distance */
          Rsqij += Rij[kdim]*Rij[kdim];
        }

        /* belongs to current layer or not*/
        if (Rsqij < cutsq_layer) {
          if (in_layer[j] != -1) {/* already been included in some layer */
            if(in_layer[j] != in_layer[i]){ /* i and j in different layer */
                                            /* if the choice of cutsq_layer is
                                               appropriate, this should not happen. */
              ier = KIM_STATUS_FAIL;
              KIM_API_report_error(__LINE__, __FILE__, "atom already in some layer", ier);
              return ier;
            }
          }
          else
          {
            nin += 1;
            layer[nin-1] = j;
            in_layer[j] = currentLayer;
          }
        }

        /* it would be more appropriate to place the following piece outside the
           if(Rsqij<cutsq_layer) statement. But it is fine as long as cutsq_layer is
           large enough to include 3 neighors. We do this to increase the speed. */
        /* find the 3 nearest neigh */
        if (Rsqij < neigh1_Rsq) {
          for (k=0; k<DIM; k++) {
            neigh3_Rvec[k] = neigh2_Rvec[k];
            neigh2_Rvec[k] = neigh1_Rvec[k];
            neigh1_Rvec[k] = Rij[k];
          }
          neigh3_Rsq = neigh2_Rsq;
          neigh2_Rsq = neigh1_Rsq;
          neigh1_Rsq = Rsqij;
        } else if (Rsqij < neigh2_Rsq) {
          for (k=0; k<DIM; k++) {
            neigh3_Rvec[k] = neigh2_Rvec[k];
            neigh2_Rvec[k] = Rij[k];
          }
          neigh3_Rsq = neigh2_Rsq;
          neigh2_Rsq = Rsqij;
        } else if (Rsqij < neigh3_Rsq) {
          for (k=0; k<DIM; k++) {
            neigh3_Rvec[k] = Rij[k];
          }
          neigh3_Rsq = Rsqij;
        }


/* NOTE delete */
/*
      if(i==702) {
        printf("%d 1 %f %f %f 0 0 0\n", kk++,coords[j*DIM+0],  coords[j*DIM+1], coords[j*DIM+2]);
      }
*/
      }

      /* note padding atoms () have 0 neigh, and do not contribute energy */
      if (numOfAtomNeigh != 0) {
        /* check whether we have enough neighs to compute normal*/
        if (neigh1_Rsq >= 1.0e10 || neigh2_Rsq >= 1.0e10 ||  neigh3_Rsq >= 1.0e10) {
          ier = KIM_STATUS_FAIL;
          KIM_API_report_error(__LINE__, __FILE__, "no enough neigh to construct normal", ier);
          ier = KIM_STATUS_FAIL;
          return ier;
        }

        /* compute normal */
        get_normal(neigh1_Rvec, neigh2_Rvec, neigh3_Rvec, tmp_normal);
        for (k=0; k<DIM; k++) {
          normal[i][k] = tmp_normal[k];
        }
      }

      /* get to the next atom in current layer */
      ii++;

      if (ii == nin) break;

    } /* end while finding one layer */

      nremain -= nin;
      if (nremain == 0) break;

      nlayers += 1;

  } /* end while finding all layers */




  /*********************************************/
  /* Compute enery and forces */

  /* loop over particles and compute enregy and forces */
  i = -1;
  while( 1 )
  {

    /* Set up neighbor list for next atom for all NBC methods */
    if (1 == IterOrLoca) /* ITERATOR mode */
    {
      ier = (*get_neigh)(&pkim, &zero, &one, &currentAtom, &numOfAtomNeigh,
          &neighListOfCurrentAtom, &Rij_list);
      if (KIM_STATUS_NEIGH_ITER_PAST_END == ier) /* the end of the list, terminate loop */
      {
        break;
      }
      if (KIM_STATUS_OK > ier) /* some sort of problem, exit */
      {
        KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_neigh", ier);
        return ier;
      }

      i = currentAtom + model_index_shift;
    }
    else
    {
      i++;

/*printf("natoms=%d, cutsq=%f, i=%d, particleStat=%d\n",*nAtoms, cutsq, i,particleStatus[i]);
*/
      if (*nAtoms <= i) /* incremented past end of list, terminate loop */
      {
        break;
      }

      if (!particleStatus[i]) continue; /* not contributing atom */

      if (3 == NBC) /* CLUSTER NBC method */
      {
        numOfAtomNeigh = *nAtoms - 1;
        for (kdim = 0; kdim < *nAtoms; ++kdim)
        {
          if (kdim < i)
            neighListOfCurrentAtom[kdim] = kdim - model_index_shift;
          if (kdim > i)
            neighListOfCurrentAtom[kdim - 1] = kdim - model_index_shift;
        }
        ier = KIM_STATUS_OK;
      }
      else
      {
        request = i - model_index_shift;
        ier = (*get_neigh)(&pkim, &one, &request, &currentAtom, &numOfAtomNeigh,
            &neighListOfCurrentAtom, &Rij_list);
        if (KIM_STATUS_OK != ier) /* some sort of problem, exit */
        {
          KIM_API_report_error(__LINE__, __FILE__, "get_neigh", ier);
          ier = KIM_STATUS_FAIL;
          return ier;
        }
      }
    }
    iSpecies = particleSpecies[i];
    ilayer = in_layer[i];


    /* loop over the neighbors of atom i */
    for (jj = 0; jj < numOfAtomNeigh; ++ jj) {

      j = neighListOfCurrentAtom[jj] + model_index_shift; /* get neighbor ID */
      jSpecies = particleSpecies[j];
      jlayer = in_layer[j];

      /* atoms in the same layer do not interact */
      if (ilayer == jlayer) continue;

      /* get corresponding parameters */
      if (iSpecies == SPEC1 && jSpecies == SPEC1) {
        inter_idx = 0;
      } else if (iSpecies == SPEC2 && jSpecies == SPEC2) {
        inter_idx = 2;
      } else {
        inter_idx = 1;
      }

      /* get two body parameters */
      cutoff = (buffer->cutoff)[inter_idx];
      cutsq  = (buffer->cutsq)[inter_idx];
      C0     = (buffer->C0)[inter_idx];
      C2     = (buffer->C2)[inter_idx];
      C4     = (buffer->C4)[inter_idx];
      C      = (buffer->C)[inter_idx];
      delta  = (buffer->delta)[inter_idx];
      lambda = (buffer->lambda)[inter_idx];
      A      = (buffer->A)[inter_idx];
      z0     = (buffer->z0)[inter_idx];



      /* compute relative position vector and squared distance */
      Rsqij = 0.0;
      for (kdim = 0; kdim < DIM; ++kdim)
      {
        if (0 != NBC) /* all methods except NEIGH_RVEC */
        {
          Rij[kdim] = coords[j*DIM + kdim] - coords[i*DIM + kdim];
        }
        else          /* NEIGH_RVEC method */
        {
          Rij[kdim] = Rij_list[jj*DIM + kdim];
        }

        /* apply periodic boundary conditions if required */
        if (2 == NBC)
        {
          if (fabs(Rij[kdim]) > 0.5*boxSideLengths[kdim])
          {
            Rij[kdim] -= (Rij[kdim]/fabs(Rij[kdim]))*boxSideLengths[kdim];
          }
        }

        /* compute squared distance */
        Rsqij += Rij[kdim]*Rij[kdim];
      }
      R = sqrt(Rsqij);

      /* compute energy and force */
      if (Rsqij > cutsq) continue; /* particles are not interacting  */

/*      if (comp_force || comp_process_dEdr)
      {
*/
        /* compute pair potential and its derivative */
/*        calc_phi_dphi_two(*A, *B, *p, *q, *cutoff, *sigma,
            R1, &phi_two, &dphi_two);

        dEidr_two = 0.5*dphi_two;
      }
      else
 */
      {
        /* compute just pair potential */
        calc_phi_repulsive(C0, C2, C4, C, delta, lambda, z0, R, Rij, normal[i],
            normal[j], &phi_repul);
          calc_phi_attractive(A, z0, R, &phi_attract);
      }

      /* contribution to energy */
      if (comp_particleEnergy)
      {
        particleEnergy[i] += 0.5*(phi_repul + phi_attract);
      }
      if (comp_energy)
      {
        *energy += 0.5*(phi_repul + phi_attract);
      }

      /* contribution to process_dEdr */
/*      if (comp_process_dEdr)
      {
        ier = KIM_API_process_dEdr(km, &dEidr_two, &R1, &pRij, &i, &j);
      }
*/

      /* contribution to forces */
/*      if (comp_force) {
        for (kdim = 0; kdim < DIM; ++kdim) {
          force[i*DIM + kdim] += (dEidr_repul + dEidr_attract)*Rij[kdim]/R;
          force[j*DIM + kdim] -= (dEidr_repul + dEidr_attract)*Rij[kdim]/R;
        }
      }
*/

    } /* loop on jj */


  }    /* infinite while loop terminated by break statements above */

  /* Free temporary storage */
  if (3 == NBC)
  {
    free(neighListOfCurrentAtom);
  }


  /* free local memory */
  free(layer);
  free(in_layer);
  free(normal);


  /* everything is great */
  ier = KIM_STATUS_OK;

  return ier;
}

/* Initialization function */
int model_driver_init(void *km, char* paramfile_names, int* nmstrlen, int* numparamfiles)
{
  /* KIM variables */
  intptr_t* pkim = *((intptr_t**) km);
  char* paramfile1name;

  /* Local variables */
  FILE* fid;
  double* model_cutoff;
  double* model_cutsq;
  double* model_C0;
  double* model_C2;
  double* model_C4;
  double* model_C;
  double* model_delta;
  double* model_lambda;
  double* model_A;
  double* model_z0;

  struct model_buffer* buffer;
  const char* NBCstr;
  double* cutoff;
  int nSpecies;
  int nInteract;
  int ier;
  int i;

  /* set paramfile1name */
  if (*numparamfiles != 1) {
    ier = KIM_STATUS_FAIL;
    KIM_API_report_error(__LINE__, __FILE__, "Incorrect number of parameter files.", ier);
    return ier;
  }
  paramfile1name = paramfile_names;

  /* store pointer to functions in KIM object */
  KIM_API_setm_method(pkim, &ier, 3*4,
      "compute", 1, &compute, 1,
      "reinit",  1, &reinit,  1,
      "destroy", 1, &destroy, 1);
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_setm_method", ier);
    return ier;
  }

  /* Read in model parameters from parameter file */
  fid = fopen(paramfile1name, "r");
  if (fid == NULL) {
    ier = KIM_STATUS_FAIL;
    KIM_API_report_error(__LINE__, __FILE__, "Unable to open parameter file", ier);
    return ier;
  }

  /* read number of species */
  ier = fscanf(fid, "%d\n", &nSpecies);
  if (ier != 1) {
    ier = KIM_STATUS_FAIL;
    KIM_API_report_error(__LINE__, __FILE__, "error reading first line of parameter file", ier);
    return ier;
  }
  nInteract = (nSpecies + 1)*nSpecies/2;

  /* allocate memory for parameters */
  model_C0     = (double*) malloc(nInteract*sizeof(double));
  model_C2     = (double*) malloc(nInteract*sizeof(double));
  model_C4     = (double*) malloc(nInteract*sizeof(double));
  model_C      = (double*) malloc(nInteract*sizeof(double));
  model_delta  = (double*) malloc(nInteract*sizeof(double));
  model_lambda = (double*) malloc(nInteract*sizeof(double));
  model_A      = (double*) malloc(nInteract*sizeof(double));
  model_z0     = (double*) malloc(nInteract*sizeof(double));
  model_cutoff = (double*) malloc(nInteract*sizeof(double));
  model_cutsq  = (double*) malloc(nInteract*sizeof(double));
  if (model_C0==NULL
      || model_C2==NULL
      || model_C4==NULL
      || model_C==NULL
      || model_delta==NULL
      || model_lambda==NULL
      || model_A==NULL
      || model_z0==NULL
      || model_cutoff==NULL
      || model_cutsq==NULL) {
    ier = KIM_STATUS_FAIL;
    KIM_API_report_error(__LINE__, __FILE__, "malloc", ier);
    return ier;
  }

  /* read parameters */
  for (i=0; i< nInteract; ++i) {
    ier = fscanf(fid, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
        &model_C0[i],
        &model_C2[i],
        &model_C4[i],
        &model_C[i],
        &model_delta[i],
        &model_lambda[i],
        &model_A[i],
        &model_z0[i],
        &model_cutoff[i]);
    /* check that we read the right number of parameters */
    if (9 != ier) {
      ier = KIM_STATUS_FAIL;
      KIM_API_report_error(__LINE__, __FILE__, "corrupted parameter file", ier);
      return ier;
    }
  }

  /* close param file */
  fclose(fid);

  /* convert units */
  for (i=0; i< nInteract; ++i) {
    model_delta[i] *= KIM_API_convert_to_act_unit(pkim, "A", "eV", "e", "K", "fs",
        1.0, 0.0, 0.0, 0.0, 0.0, &ier);
    if (KIM_STATUS_OK > ier) {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_convert_to_act_unit", ier);
      return ier;
    }

    model_lambda[i] *= KIM_API_convert_to_act_unit(pkim, "A", "eV", "e", "K", "fs",
        -1.0, 0.0, 0.0, 0.0, 0.0, &ier);
    if (KIM_STATUS_OK > ier) {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_convert_to_act_unit", ier);
      return ier;
    }

    model_z0[i] *= KIM_API_convert_to_act_unit(pkim, "A", "eV", "e", "K", "fs",
        1.0, 0.0, 0.0, 0.0, 0.0, &ier);
    if (KIM_STATUS_OK > ier) {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_convert_to_act_unit", ier);
      return ier;
    }
  }

  /* store model cutoff in KIM object */
  cutoff = (double*) KIM_API_get_data(pkim, "cutoff", &ier);
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_data", ier);
    return ier;
  }
  *cutoff = 0.0;
  for (i=0; i< nInteract; ++i) {
    model_cutsq[i] = model_cutoff[i]*model_cutoff[i];
    if ( model_cutoff[i] > *cutoff)
      *cutoff = model_cutoff[i];
  }

  /* store parameters in KIM object */
  KIM_API_setm_data(pkim, &ier, 9*4,
      "PARAM_FREE_cutoff",    nInteract,  model_cutoff,  1,
      "PARAM_FREE_C0",        nInteract,  model_C0,      1,
      "PARAM_FREE_C2",        nInteract,  model_C2,      1,
      "PARAM_FREE_C4",        nInteract,  model_C4,      1,
      "PARAM_FREE_C",         nInteract,  model_C,       1,
      "PARAM_FREE_delta",     nInteract,  model_delta,   1,
      "PARAM_FREE_lambda",    nInteract,  model_lambda,  1,
      "PARAM_FREE_A",         nInteract,  model_A,       1,
      "PARAM_FREE_z0",        nInteract,  model_z0,      1);
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_setm_data", ier);
    return ier;
  }

  /* allocate buffer */
  buffer = (struct model_buffer*) malloc(sizeof(struct model_buffer));
  if (NULL==buffer) {
    ier = KIM_STATUS_FAIL;
    KIM_API_report_error(__LINE__, __FILE__, "malloc", ier);
    return ier;
  }

  /* setup buffer */
  /* Determine neighbor list boundary condition (NBC) */
  ier = KIM_API_get_NBC_method(pkim, &NBCstr);
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_NBC_method", ier);
    return ier;
  }
  if (!strcmp("NEIGH_RVEC_F",NBCstr)) {
    buffer->NBC = 0;
  }
  else if (!strcmp("NEIGH_PURE_F",NBCstr)) {
    buffer->NBC = 1;
  }
  else if (!strcmp("MI_OPBC_F",NBCstr)) {
    buffer->NBC = 2;
  }
  else if (!strcmp("CLUSTER",NBCstr)) {
    buffer->NBC = 3;
  }
  else {
    ier = KIM_STATUS_FAIL;
    KIM_API_report_error(__LINE__, __FILE__, "Unknown NBC method", ier);
    return ier;
  }

  /* determine neighbor list handling mode */
  if (buffer->NBC != 3) {
    /*****************************
     * IterOrLoca = 1 -- Iterator
     *            = 2 -- Locator
     *****************************/
    buffer->IterOrLoca = KIM_API_get_neigh_mode(pkim, &ier);
    if (KIM_STATUS_OK > ier) {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_neigh_mode", ier);
      return ier;
    }
    if ((buffer->IterOrLoca != 1) && (buffer->IterOrLoca != 2)) {
      printf("* ERROR: Unsupported IterOrLoca mode = %i\n", buffer->IterOrLoca);
      return ier;
    }
  }
  else {
    buffer->IterOrLoca = 2;   /* for CLUSTER NBC */
  }

  buffer->model_index_shift = KIM_API_get_model_index_shift(pkim);

  KIM_API_getm_index(pkim, &ier, 13*3,
      "cutoff",                      &(buffer->cutoff_ind),                      1,
      "numberOfParticles",           &(buffer->numberOfParticles_ind),           1,
      "numberOfSpecies",             &(buffer->numberOfSpecies_ind),             1,
      "particleSpecies",             &(buffer->particleSpecies_ind),             1,
      "particleStatus",              &(buffer->particleStatus_ind),              1,
      "coordinates",                 &(buffer->coordinates_ind),                 1,
      "get_neigh",                   &(buffer->get_neigh_ind),                   1,
      "boxSideLengths",              &(buffer->boxSideLengths_ind),              1,
      "energy",                      &(buffer->energy_ind),                      1,
      "forces",                      &(buffer->forces_ind),                      1,
      "particleEnergy",              &(buffer->particleEnergy_ind),              1,
      "process_dEdr",                &(buffer->process_dEdr_ind),                1,
      "process_d2Edr2",              &(buffer->process_d2Edr2_ind),              1);
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_index", ier);
    return ier;
  }

  /* store parameters in buffer */
  buffer->cutoff = model_cutoff;
  buffer->cutsq  = model_cutsq;
  buffer->C0     = model_C0;
  buffer->C2     = model_C2;
  buffer->C4     = model_C4;
  buffer->C      = model_C;
  buffer->delta  = model_delta;
  buffer->lambda = model_lambda;
  buffer->A      = model_A;
  buffer->z0     = model_z0;
  /* end setup buffer */

  /* store in model buffer */
  KIM_API_set_model_buffer(pkim, (void*) buffer, &ier);
  if (KIM_STATUS_OK > ier)
  {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_set_model_buffer", ier);
    return ier;
  }

  ier = KIM_STATUS_OK;
  return ier;
}

/* Reinitialization function */
static int reinit(void *km)
{
  /* Local variables */
  intptr_t* pkim = *((intptr_t**) km);
  double *cutoff;
  struct model_buffer* buffer;
  int* nSpecies;
  int nInteract;
  int ier;
  int i;

  /* get buffer from KIM object */
  buffer = (struct model_buffer*) KIM_API_get_model_buffer(pkim, &ier);
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_model_buffer", ier);
    return ier;
  }

  nSpecies = KIM_API_get_data_by_index(pkim, buffer->numberOfSpecies_ind, &ier);
  nInteract = (*nSpecies + 1)*(*nSpecies)/2;

  /* update cutoff in KIM API and also cutsq */
  cutoff = KIM_API_get_data(pkim, "cutoff", &ier);
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_data", ier);
    return ier;
  }
  *cutoff = 0.0;
  for (i=0; i< nInteract; ++i) {
    buffer->cutsq[i] = buffer->cutoff[i]*buffer->cutoff[i];
    if ( buffer->cutoff[i] > *cutoff) {
      *cutoff = buffer->cutoff[i];
    }
  }

  ier = KIM_STATUS_OK;
  return ier;
}

/* destroy function */
static int destroy(void *km)
{
  /* Local variables */
  intptr_t* pkim = *((intptr_t**) km);
  struct model_buffer* buffer;
  int ier;

  /* get model buffer from KIM object */
  buffer = (struct model_buffer*) KIM_API_get_model_buffer(pkim, &ier);
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_model_buffer", ier);
    return ier;
  }

  /* free parameters */
  free(buffer->cutoff);
  free(buffer->cutsq);
  free(buffer->C0);
  free(buffer->C2);
  free(buffer->C4);
  free(buffer->C);
  free(buffer->delta);
  free(buffer->lambda);
  free(buffer->A);
  free(buffer->z0);

  /* destroy the buffer */
  free(buffer);

  ier = KIM_STATUS_OK;
  return ier;
}
