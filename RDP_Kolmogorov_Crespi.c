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
* typedef and function prototype
*******************************************************************************/

#define HALF 0.5
#define DIM 3       /* dimensionality of space */
#define SPEC1 1     /* internal species code */
#define SPEC2 2     /* internal species code */

typedef double vectorDIM[DIM];

typedef int (*get_neigh_ptr)(void *,int *,int *,int *, int *, int **, double **);

typedef struct {
  int model_index_shift;

  /* compute flag */
  int comp_energy;
  int comp_forces;
  int comp_particleEnergy;
  int comp_process_dEdr;

  /* model input output */
  int* nAtoms;
  int* nSpecies;
  int* particleSpecies;
  int* particleStatus;
  double* energy;
  double* particleEnergy;
  vectorDIM* coords;
  vectorDIM* forces;
  get_neigh_ptr get_neigh;

  /* layers */
  int* in_layer; /* atoms in which layer. -1: not classified in any layer */
  int (*nearest3neigh)[3];

  /* derivatives of ni (normal at atom i) w.r.t ri and its 3 neighbors */
  int nbi1;
  int nbi2;
  int nbi3;
  vectorDIM ni;
  vectorDIM dni_dri[DIM];
  vectorDIM dni_drnb1[DIM];
  vectorDIM dni_drnb2[DIM];
  vectorDIM dni_drnb3[DIM];

  double* cutoff;
  double* C0;
  double* C2;
  double* C4;
  double* C;
  double* delta;
  double* lambda;
  double* A;
  double* z0;
  double* B;
  double* eta;
} model_buffer;

/* function prototypes that all Model Drivers need */
int model_driver_init(void* km, char* paramfile_names, int* nmstrlen, int* numparamfiles);
static int reinit(void* km);
static int destroy(void* km);
static int compute(void* km);

/* local function prototypes */
static int init_compute(void* km);
static int create_layers(void* km);
static int calc_energy_forces(void* km);
static int free_data(void* km);

static double calc_attractive(model_buffer *const buffer, const int i, const int j,
    double *const rij, const double r);

static double calc_repulsive(model_buffer *const buffer, const int i, const int j,
    double *const rij, const double r);

static void normal(model_buffer *const buffer, const int i, int *const nbi1,
    int *const nbi2, int *const nbi3, double normal[DIM], double dn_dri[DIM][DIM],
    double dn_drk[DIM][DIM], double dn_drl[DIM][DIM], double dn_drm[DIM][DIM]);

static double td(double C0, double C2, double C4, double delta,
    const double* rvec, double r, const double* n,  double *const rho_sq,
    double *const dtd);

static void get_drhosqij(const double rij[DIM], const double ni[DIM],
    double dni_dri[DIM][DIM], double dni_drn1[DIM][DIM],
    double dni_drn2[DIM][DIM], double dni_drn3[DIM][DIM],
    double drhosq_dri[DIM], double drhosq_drj[DIM],
    double drhosq_drn1[DIM], double drhosq_drn2[DIM],
    double drhosq_drn3[DIM]);

static double  dihedral(model_buffer *const buffer, double rhosq, const int i, const int j,
    double *const d_drhosq, double d_dri[DIM], double d_drj[DIM],
    double d_drk1[DIM], double d_drk2[DIM], double d_drk3[DIM],
    double d_drl1[DIM], double d_drl2[DIM], double d_drl3[DIM]);

static double deriv_cos_omega(double rk[DIM], double ri[DIM], double rj[DIM],
  double rl[DIM], double dcos_drk[DIM], double dcos_dri[DIM],
  double dcos_drj[DIM], double dcos_drl[DIM]);


static double tap(model_buffer* buffer, int i, int j, double r, double *const dtap);

static double tap_rho(double rhosq, double cut_rhosq, double *const drhosq);

/* helper */
static double rsq_rij(model_buffer* const buffer, const int i, const int j,
    double *const rij);

static int param_index(model_buffer* const buffer, const int i, const int j);

double dot(const double x[DIM], const double y[DIM]);

static void deriv_cross(double rk[DIM], double rl[DIM], double rm[DIM], double cross[DIM],
    double dcross_drk[DIM][DIM], double dcross_drl[DIM][DIM], double dcross_drm[DIM][DIM]);

/*static void cross(const double x[DIM], const double y[DIM], double rslt[DIM]);
*/
static void mat_dot_vec(double X[DIM][DIM], const double y[DIM], double z[DIM]);


/******************************************************************************
* function definitions
******************************************************************************/

static int compute(void* km)
{
  int ier;

  init_compute(km);
  create_layers(km);
  calc_energy_forces(km);
  free_data(km);

  /* everything is great */
  ier = KIM_STATUS_OK;

  return ier;
}

/* store KIM API pointers in buffer */
static int init_compute(void* km)
{
  /* local variables */
  intptr_t* pkim = *((intptr_t**) km);

  model_buffer* buffer;
  int i;
  int ier;


  /* get buffer from KIM object */
  buffer = (model_buffer*) KIM_API_get_model_buffer(pkim, &ier);
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_model_buffer", ier);
    return ier;
  }

  /* check to see if we have been asked to compute the forces, particleEnergy, and dEdr */
  KIM_API_getm_compute(pkim, &ier, 4*3,
      "energy",          &buffer->comp_energy,          1,
      "forces",          &buffer->comp_forces,          1,
      "particleEnergy",  &buffer->comp_particleEnergy,  1,
      "process_dEdr",    &buffer->comp_process_dEdr,    1);
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_compute", ier);
    return ier;
  }

  KIM_API_getm_data(pkim, &ier, 8*3,
      "numberOfParticles",  &buffer->nAtoms,           1,
      "numberOfSpecies",    &buffer->nSpecies,         1,
      "particleSpecies",    &buffer->particleSpecies,  1,
      "particleStatus",     &buffer->particleStatus,   1,
      "coordinates",        &buffer->coords,           1,
      "energy",             &buffer->energy,           buffer->comp_energy,
      "forces",             &buffer->forces,           buffer->comp_forces,
      "particleEnergy",     &buffer->particleEnergy,   buffer->comp_particleEnergy);
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_data", ier);
    return ier;
  }

  /* neigh access funtion */
  buffer->get_neigh = (get_neigh_ptr) KIM_API_get_method(pkim, "get_neigh", &ier);
  if (KIM_STATUS_OK > ier) {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_method", ier);
      return ier;
  }

  /* initialize energy and forces */
  if (buffer->comp_particleEnergy) {
    for (i = 0; i < *buffer->nAtoms; ++i) {
      buffer->particleEnergy[i] = 0.0;
    }
  }
  if (buffer->comp_energy) {
    *buffer->energy = 0.0;
  }
  if (buffer->comp_forces) {
    for (i = 0; i < *buffer->nAtoms; ++i) {
      buffer->forces[i][0] = 0.0;
      buffer->forces[i][1] = 0.0;
      buffer->forces[i][2] = 0.0;
    }
  }

  ier = KIM_STATUS_OK;
  return ier;
}


/* create layers and find the nearest 3 neighbors of each atom */
static int create_layers(void* km)
{
  /* local variables */
  intptr_t* pkim = *((intptr_t**) km);

  model_buffer* buffer;
  int nAtoms;
  int model_index_shift;

  /* neighbor list */
  get_neigh_ptr get_neigh;
  int currentAtom;
  int request;
  int num_neigh;
  int* ilist;
  double* Rij_list;
  int one = 1;

  /* layers and nearest3neighs */
  int* in_layer; /* atoms in which layer , -1: not classified in any layer */
  int* layer; /* layer contains which atoms */
  int (*nearest3neigh)[3];
  int currentLayer;
  int nin; /* number of atoms in current layer */
  int nremain; /* number of atoms not included in any layer */
  int nlayers;
  int nb1, nb2, nb3;
  double nb1_rsq, nb2_rsq, nb3_rsq;
  double cutsq_layer = (0.72*3.44)*(0.72*3.44);

  vectorDIM rij;
  double rsq;
  int i,j,k,ii,jj;
  int ier;
  int more_than_one_layer;


  /* unpack data from buffer */
  buffer = (model_buffer*) KIM_API_get_model_buffer(pkim, &ier);
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_model_buffer", ier);
    return ier;
  }
  nAtoms = *buffer->nAtoms;
  model_index_shift = buffer->model_index_shift;
  get_neigh = buffer->get_neigh;

  /* no atoms, which can happen when doing comain decompisition */
  if (nAtoms <= 0) return KIM_STATUS_OK;

  /* allocate memory for layer infor*/
  layer = (int*) malloc(nAtoms*sizeof(int));
  in_layer = (int*) malloc(nAtoms*sizeof(int));
  nearest3neigh = (int(*)[3]) malloc(nAtoms*sizeof(double[3]));

  /* init atoms layer to -1 (not in any layer) */
  for (k=0; k<nAtoms; k++) {
    in_layer[k] = -1;
  }

  nlayers = 1;
  nremain = nAtoms;

  while(1) {  /* to create all layers */

    /* init all atoms in current layer have atom number -1 */
    /* We do not always have nAtoms (except the very first loop), we allocate layer
       to the size of nAtom just because we do not want to reallocate it each time. */
    for (k=0; k<nAtoms; k++) {
      layer[k] = -1;
    }

    /* find an atom not incldued in any layer and start with it */
    currentLayer = nlayers - 1;
    for (k=0; k<nAtoms; k++) {
      if (in_layer[k] == -1) {
        in_layer[k] = currentLayer;
        layer[0] = k; /* first atom in current layer */
        break;
      }
    }

    nin = 1;
    ii = 0;
    while(1) { /* to find all atoms in currentLayer */

      i = layer[ii];

      /* get neighbors of atom i */
      request = i - model_index_shift;
      ier = (*get_neigh)(&pkim, &one, &request, &currentAtom, &num_neigh, &ilist, &Rij_list);
      if (KIM_STATUS_OK != ier) {
        KIM_API_report_error(__LINE__, __FILE__, "get_neigh", ier);
        ier = KIM_STATUS_FAIL;
        return ier;
      }

      /* init nb1 to be the 1st nearest neigh, nb3 the 3rd nearest */
      nb1 = -1;
      nb2 = -1;
      nb3 = -1;
      nb1_rsq = 1.1e10;
      nb2_rsq = 2.0e10;
      nb3_rsq = 3.0e10;

      /* loop over the neighbors of atom i */
      for (jj=0; jj<num_neigh; ++jj) {

        j = ilist[jj] + model_index_shift; /* get neighbor ID */

        /* compute relative position vector and squared distance */
        rsq = rsq_rij(buffer, i, j, rij);


        /* belongs to current layer */
        if (rsq < cutsq_layer) {
          if (in_layer[j] == -1) { /* has not been included in some layer */
            nin += 1;
            layer[nin-1] = j;
            in_layer[j] = currentLayer;
          } else {
            if(in_layer[j] != in_layer[i]){ /* i and j in different layer. If the choice
                                               of cutsq_layer is appropriate, this should not happen. */
              ier = KIM_STATUS_FAIL;
              KIM_API_report_error(__LINE__, __FILE__, "atom already in some layer", ier);
              return ier;
            }
          }
        }

        /* find the 3 nearest neigh */
        if (rsq < nb1_rsq) {
          nb3 = nb2;
          nb2 = nb1;
          nb1 = j;
          nb3_rsq = nb2_rsq;
          nb2_rsq = nb1_rsq;
          nb1_rsq = rsq;
        } else if (rsq < nb2_rsq) {
          nb3 = nb2;
          nb2 = j;
          nb3_rsq = nb2_rsq;
          nb2_rsq = rsq;
        } else if (rsq < nb3_rsq) {
          nb3 = j;
          nb3_rsq = rsq;
        }

      } /* loop on jj */


      /* check whether we have enough neighs to compute normal*/
      if (nb1_rsq >= 1.0e10 || nb2_rsq >= 1.0e10 ||  nb3_rsq >= 1.0e10) {
        ier = KIM_STATUS_FAIL;
        KIM_API_report_error(__LINE__, __FILE__, "no enough neigh to construct normal", ier);
        ier = KIM_STATUS_FAIL;
        return ier;
      }

      /* store nearest 3 neigh. This is for latter use to comptue the derivatives
         of normal; Naively, we can compute the derivatives here, but each atom has
         3 nerest neighs, and the derivative with respect to each neigh has 9
         component, resulting in storing an array of size N*(1+3)*9. Seems not good.*/
      nearest3neigh[i][0] = nb1;
      nearest3neigh[i][1] = nb2;
      nearest3neigh[i][2] = nb3;

      /* get to the next atom in current layer */
      ii++;
      if (ii == nin) break;

    } /* end while finding one layer */

      nremain -= nin;
      if (nremain == 0) break;
      nlayers += 1;
  } /* end while finding all layers */

  /* whether all atoms in the same layer? */
  more_than_one_layer = 0;
  currentLayer = 0;
  for (i=0; i<nAtoms; i++) {
    if (in_layer[i] != currentLayer) { /* find a atom not in current layer */
      more_than_one_layer = 1;
      break;
    }
  }
  if (! more_than_one_layer) {
    ier = KIM_STATUS_FAIL;
    KIM_API_report_error(__LINE__, __FILE__, "Only one layer detected. The layer "
    "sepration seems to small. You can either use larger layer separation or "
    "tune `cutsq_layer' in the Model.", ier);
    ier = KIM_STATUS_FAIL;
    return ier;
  }


  /* store in_layer nearest3neigh in bubber */
  buffer->in_layer = in_layer;
  buffer->nearest3neigh = nearest3neigh;

  /* dealloate */
  free(layer);

  return KIM_STATUS_OK;
}


/* contribution to energy and forces */
static int calc_energy_forces(void* km)
{
  intptr_t* pkim = *((intptr_t**) km);

  model_buffer* buffer;
  int model_index_shift;
  int nAtoms;
  double* energy;
  double* particleEnergy;
  int* particleStatus;
  int* in_layer;

  /* neighbor list */
  get_neigh_ptr get_neigh;
  int currentAtom;
  int request;
  int num_neigh;
  int* ilist;
  double* Rij_list;
  int one = 1;

  double rij[DIM];
  double rsq;
  double r;
  double phi_attr, phi_repul;
  int inter_idx;
  double cutoff;
  int ilayer, jlayer;
  int i, j, k, jj;
  int ier;

  /* get buffer from kim object */
  buffer = (model_buffer*) KIM_API_get_model_buffer(pkim, &ier);
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_model_buffer", ier);
    return ier;
  }
  /* unpack data from buffer */
  model_index_shift = buffer->model_index_shift;
  nAtoms = *buffer->nAtoms;
  particleStatus = buffer->particleStatus;
  particleEnergy = buffer->particleEnergy;
  energy = buffer->energy;
  in_layer = buffer->in_layer;
  get_neigh = buffer->get_neigh;

  /* loop over all particles */
  for (i=0; i<nAtoms; i++) {

    /* non-contriburing atom */
    if(!particleStatus[i]) continue;

    ilayer = in_layer[i];

    /* get neighbors of atom i */
    request = i - model_index_shift;
    ier = (*get_neigh)(&pkim, &one, &request, &currentAtom, &num_neigh, &ilist, &Rij_list);
    if (KIM_STATUS_OK != ier) {
      KIM_API_report_error(__LINE__, __FILE__, "get_neigh", ier);
      ier = KIM_STATUS_FAIL;
      return ier;
    }

    /* get 3 nearest neighs of atom i and compute the derivative of ni w.r.t. them */
    normal(buffer, i, &buffer->nbi1, &buffer->nbi2, &buffer->nbi3, buffer->ni,
        buffer->dni_dri, buffer->dni_drnb1, buffer->dni_drnb2, buffer->dni_drnb3);

    /* loop over the neighbors of atom i */
    for (jj=0; jj<num_neigh; ++jj) {

      /* get neighbor ID */
      j = ilist[jj] + model_index_shift;

      jlayer = in_layer[j];

      /* atoms in the same layer do not interact */
      if (ilayer == jlayer) continue;

      /* get params */
      inter_idx = param_index(buffer, i, j);
      cutoff = buffer->cutoff[inter_idx];

      /* relative position vector and squared distance */
      rsq = rsq_rij(buffer, i, j, rij);
      r = sqrt(rsq);

      /* compute energy and forces */
      phi_repul = 0.0;
      phi_attr = 0.0;
      if (rsq < cutoff*cutoff) {
        phi_repul = calc_repulsive(buffer, i, j, rij, r);
        phi_attr = calc_attractive(buffer, i, j, rij, r);

        if (buffer->comp_energy) {
          /* HALF because of full neighborlist */
          *energy += HALF * (phi_repul + phi_attr);
        }
        if (buffer->comp_particleEnergy) {
          for (k=0; k<nAtoms; k++) {
            particleEnergy[k] +=  HALF * (phi_repul + phi_attr);
          }
        }

      }

    }  /* loop on jj */
  }  /* loop on i */

  ier = KIM_STATUS_OK;
  return ier;
}


/* the r^{-6} attractive part */
static double calc_attractive(model_buffer *const buffer, const int i, const int j,
    double *const rij, const double r)
{
  vectorDIM* forces;

  /* params */
  int inter_idx;
  double A;
  double z0;

  double tp;
  double dtp;
  double r6;
  double dr6;

  double roz0_sq;
  double phi;
  double fpair;
  int k;

  forces = buffer->forces;

  /* get corresponding parameters */
  inter_idx = param_index(buffer, i, j);
  z0 = buffer->z0[inter_idx];
  A = buffer->A[inter_idx];

  roz0_sq = r*r/(z0*z0);
  tp = tap(buffer, i, j, r, &dtp);
  r6 = A/(roz0_sq*roz0_sq*roz0_sq);
  dr6 = -6*r6/r;

  phi = - r6*tp;

  /* forces */
  if (buffer->comp_forces) {
    fpair = - HALF * (r6*dtp + dr6*tp);
    for (k=0; k<DIM; k++) {
      forces[i][k] += rij[k]*fpair/r;
      forces[j][k] -= rij[k]*fpair/r;
    }
  }

  return phi;
}


/* the repulsive part */
static double calc_repulsive(model_buffer *const buffer, const int i, const int j,
    double *const rij, const double r)
{
  vectorDIM* forces;

  /* params */
  int inter_idx;
  double C0;
  double C2;
  double C4;
  double C;
  double delta;
  double lambda;
  double z0;

  /* nearest 3 neigh of i */
  int nbi1;
  int nbi2;
  int nbi3;
  double* ni;
  vectorDIM* dni_dri;
  vectorDIM* dni_drnb1;
  vectorDIM* dni_drnb2;
  vectorDIM* dni_drnb3;

  /* nearest 3 neigh of j*/
  int nbj1;
  int nbj2;
  int nbj3;
  vectorDIM nj;
  vectorDIM dnj_drj[DIM];
  vectorDIM dnj_drnb1[DIM];
  vectorDIM dnj_drnb2[DIM];
  vectorDIM dnj_drnb3[DIM];

  double V1, dV1;
  double V2;
  double tp, dtp;
  double tdij, tdji;
  double rhosqij, rhosqji;
  double dtdij, dtdji;
  double gij, gji;
  double dgij_drhosq, dgji_drhosq;
  vectorDIM dgij_dri;
  vectorDIM dgij_drj;
  vectorDIM dgij_drk1;
  vectorDIM dgij_drk2;
  vectorDIM dgij_drk3;
  vectorDIM dgij_drl1;
  vectorDIM dgij_drl2;
  vectorDIM dgij_drl3;
  vectorDIM dgji_dri;
  vectorDIM dgji_drj;
  vectorDIM dgji_drk1;
  vectorDIM dgji_drk2;
  vectorDIM dgji_drk3;
  vectorDIM dgji_drl1;
  vectorDIM dgji_drl2;
  vectorDIM dgji_drl3;

  vectorDIM drhosqij_dri;
  vectorDIM drhosqij_drj;
  vectorDIM drhosqij_drnb1;
  vectorDIM drhosqij_drnb2;
  vectorDIM drhosqij_drnb3;
  vectorDIM drhosqji_dri;
  vectorDIM drhosqji_drj;
  vectorDIM drhosqji_drnb1;
  vectorDIM drhosqji_drnb2;
  vectorDIM drhosqji_drnb3;
  vectorDIM rji;

  double tmp;
  double phi;
  int k;


  /* unpack data from buffer */
  forces = buffer->forces;

  nbi1 = buffer->nbi1;
  nbi2 = buffer->nbi2;
  nbi3 = buffer->nbi3;
  ni = buffer->ni;
  dni_dri = buffer->dni_dri;
  dni_drnb1 = buffer->dni_drnb1;
  dni_drnb2 = buffer->dni_drnb2;
  dni_drnb3 = buffer->dni_drnb3;

  /* get corresponding parameters */
  inter_idx = param_index(buffer, i, j);
  C0 = buffer->C0[inter_idx];
  C2 = buffer->C2[inter_idx];
  C4 = buffer->C4[inter_idx];
  C = buffer->C[inter_idx];
  delta = buffer->delta[inter_idx];
  lambda = buffer->lambda[inter_idx];
  z0 = buffer->z0[inter_idx];

  /* rji */
  rji[0] = -rij[0];
  rji[1] = -rij[1];
  rji[2] = -rij[2];

  /* get 3 nearest neighs of atom i and compute the derivative of ni w.r.t. them */
  normal(buffer, j, &nbj1, &nbj2, &nbj3, nj, dnj_drj, dnj_drnb1, dnj_drnb2, dnj_drnb3);

  /* derivative of rhosq w.r.t coords of atoms i, j, and the nearests 3 neighs of i */
  get_drhosqij(rij, ni, dni_dri, dni_drnb1, dni_drnb2, dni_drnb3, drhosqij_dri,
      drhosqij_drj, drhosqij_drnb1, drhosqij_drnb2, drhosqij_drnb3);
  get_drhosqij(rji, nj, dnj_drj, dnj_drnb1, dnj_drnb2, dnj_drnb3, drhosqji_drj,
      drhosqji_dri, drhosqji_drnb1, drhosqji_drnb2, drhosqji_drnb3);

  /* transverse decay function f(rho) and its derivative w.r.t. rhosq */
  tdij = td(C0, C2, C4, delta, rij, r, ni, &rhosqij, &dtdij);
  tdji = td(C0, C2, C4, delta, rji, r, nj, &rhosqji, &dtdji);

  /* dihedral angle function and its derivateives */
  gij = dihedral(buffer, rhosqij, i, j, &dgij_drhosq, dgij_dri, dgij_drj,
      dgij_drk1, dgij_drk2, dgij_drk3, dgij_drl1, dgij_drl2, dgij_drl3);
  gji = dihedral(buffer, rhosqji, j, i, &dgji_drhosq, dgji_drj, dgji_dri,
      dgji_drl1, dgji_drl2, dgji_drl3, dgji_drk1, dgji_drk2, dgji_drk3);

  V2 = C + tdij + tdji + gij + gji;


  /* tap part */
/* NOTE this could be modified to pass cutoff, not passing buffer and i, j*/
  tp = tap(buffer, i, j, r, &dtp);

  /* exponential part */
  V1 = exp(-lambda*(r-z0));
  dV1 = -V1*lambda;

  phi = tp*V1*V2;

  /* forces */
  if (buffer->comp_forces) {
    for (k=0; k<DIM; k++) {

      /* forces due to derivatives of tap and V1 */
      tmp =  HALF* (dtp*V1 + tp*dV1) *V2 *rij[k]/r;
      forces[i][k] += tmp;
      forces[j][k] -= tmp;

      /* the following incldue the transverse decay part tdij or tdji, and the dihedral
         part gij or gji */
      /* derivative of V2 contribute to atoms i, j */
      forces[i][k] -= HALF*tp*V1* ( (dtdij+dgij_drhosq)*drhosqij_dri[k]
          + (dtdji+dgji_drhosq)*drhosqji_dri[k] + dgij_dri[k] + dgji_dri[k] );
      forces[j][k] -= HALF*tp*V1* ( (dtdij+dgij_drhosq)*drhosqij_drj[k]
          + (dtdji+dgji_drhosq)*drhosqji_drj[k] + dgij_drj[k] + dgji_drj[k] );


      /* derivative of V2 contribute to neighs of atom i */
      forces[nbi1][k] -= HALF*tp*V1* ( (dtdij+dgij_drhosq)*drhosqij_drnb1[k]
          + dgij_drk1[k] + dgji_drk1[k] );
      forces[nbi2][k] -= HALF*tp*V1* ( (dtdij+dgij_drhosq)*drhosqij_drnb2[k]
          + dgij_drk2[k] + dgji_drk2[k] );
      forces[nbi3][k] -= HALF*tp*V1* ( (dtdij+dgij_drhosq)*drhosqij_drnb3[k]
          + dgij_drk3[k] + dgji_drk3[k] );

      /* derivative of V2 contribute to neighs of atom j */
      forces[nbj1][k] -= HALF*tp*V1* ( (dtdji+dgji_drhosq)*drhosqji_drnb1[k]
          + dgij_drl1[k] + dgji_drl1[k]);
      forces[nbj2][k] -= HALF*tp*V1* ( (dtdji+dgji_drhosq)*drhosqji_drnb2[k]
          + dgij_drl2[k] + dgji_drl2[k]);
      forces[nbj3][k] -= HALF*tp*V1* ( (dtdji+dgji_drhosq)*drhosqji_drnb3[k]
          + dgij_drl3[k] + dgji_drl3[k]);

    }
  }

  return phi;
}


/* free data allocated in compute */
static int free_data(void* km)
{
  intptr_t* pkim = *((intptr_t**) km);
  model_buffer* buffer;
  int ier;

  /* get buffer from kim object */
  buffer = (model_buffer*) KIM_API_get_model_buffer(pkim, &ier);
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_model_buffer", ier);
    return ier;
  }

  free(buffer->in_layer);
  free(buffer->nearest3neigh);

  ier = KIM_STATUS_OK;
  return ier;
}


/* compute normal and its derivatives w.r.t atom ri, and its 3 nearest neighs k1, k2, k3*/
static void normal(model_buffer *const buffer, const int i, int *const k1,
    int* const k2, int* const k3, double normal[DIM], double dn_dri[DIM][DIM],
    double dn_drk1[DIM][DIM], double dn_drk2[DIM][DIM], double dn_drk3[DIM][DIM])
{
  /* local variables */
  int (*nearest3neigh)[3];
  vectorDIM* coords;
  int j,k;

  /* unpack data from buffer */
  nearest3neigh = buffer->nearest3neigh;
  coords = buffer->coords;

  *k1 = nearest3neigh[i][0];
  *k2 = nearest3neigh[i][1];
  *k3 = nearest3neigh[i][2];

  /* normal does not depend on i, setting to zero */
  for (j=0; j<DIM; j++) {
    for (k=0; k<DIM; k++) {
      dn_dri[j][k] = 0.0;
    }
  }

  /* get normal and derives of normal w.r.t to neighbors of  */
  deriv_cross(coords[*k1], coords[*k2], coords[*k3], normal, dn_drk1, dn_drk2, dn_drk3);
}


static void get_drhosqij(const double rij[DIM], const double ni[DIM],
    double dni_dri[DIM][DIM], double dni_drn1[DIM][DIM],
    double dni_drn2[DIM][DIM], double dni_drn3[DIM][DIM],
    double drhosq_dri[DIM], double drhosq_drj[DIM],
    double drhosq_drn1[DIM], double drhosq_drn2[DIM],
    double drhosq_drn3[DIM])
{

  int k;
  double ni_dot_rij = 0;
  double dni_dri_dot_rij[DIM];
  double dni_drn1_dot_rij[DIM];
  double dni_drn2_dot_rij[DIM];
  double dni_drn3_dot_rij[DIM];

  ni_dot_rij = dot(ni, rij);
  mat_dot_vec(dni_dri, rij, dni_dri_dot_rij);
  mat_dot_vec(dni_drn1, rij, dni_drn1_dot_rij);
  mat_dot_vec(dni_drn2, rij, dni_drn2_dot_rij);
  mat_dot_vec(dni_drn3, rij, dni_drn3_dot_rij);

  for (k=0; k<DIM; k++) {
    drhosq_dri[k] = -2*rij[k] -2*ni_dot_rij * (-ni[k] + dni_dri_dot_rij[k]);
    drhosq_drj[k] = 2*rij[k] - 2*ni_dot_rij*ni[k];
    drhosq_drn1[k] = -2*ni_dot_rij * dni_drn1_dot_rij[k];
    drhosq_drn2[k] = -2*ni_dot_rij * dni_drn2_dot_rij[k];
    drhosq_drn3[k] = -2*ni_dot_rij * dni_drn3_dot_rij[k];
  }

}



/* derivartive of transverse decay function f(rho) w.r.t rho */
static double td(double C0, double C2, double C4, double delta,
    const double* rvec, double r, const double* n, double *const rho_sq,
    double *const dtd)
{
  double n_dot_r;
  double del_sq;
  double rod_sq;
  double td;

  n_dot_r = dot(n, rvec);
  *rho_sq = r*r - n_dot_r*n_dot_r;

  if(*rho_sq<0) {   /* in case n is [0, 0, 1] and rho_sq is negative due to numerical error */
    *rho_sq = 0;
  }

  del_sq = delta*delta;
  rod_sq = (*rho_sq)/del_sq;
  td = exp(-rod_sq) * (C0 + rod_sq*(C2 + rod_sq*C4));
  *dtd = -td/del_sq + exp(-rod_sq) * (C2 + 2*C4*rod_sq)/del_sq;

  return td;
}


/* derivartive of dihedral angle func gij w.r.t rho, and atom positions */
static double dihedral(model_buffer *const buffer, double rhosq, const int i, const int j,
    double *const d_drhosq, double d_dri[DIM], double d_drj[DIM],
    double d_drk1[DIM], double d_drk2[DIM], double d_drk3[DIM],
    double d_drl1[DIM], double d_drl2[DIM], double d_drl3[DIM])
{
  /* params */
  int inter_idx;
  double B;
  double eta;

  int (*nearest3neigh)[3];
  int k[3];
  int l[3];
  vectorDIM* coords;

  /* local vars */
  double dihe;
  double D0;
  double D2;
  double d_drhosq_tap;
  double epart1;
  double epart2;
  double epart3;
  double cos_kl[3][3];          /* cos_omega_k1ijl1, cos_omega_k1ijl2 ...  */
  double d_dcos_kl[3][3];       /* deriv of dihedral w.r.t to cos_omega_kijl */
  double dcos_kl[3][3][4][DIM]; /* 4 indicates k, i, j, l, e.g. dcoskl[0][1][0] means
                                dcos_omega_k1ijl2 / drk */
  int m, n, dim;


  /* NOTE this should be a params and read in from file */
  double cut_rhosq = (1.42*1.1)*(1.42*1.1);

  /* if larger than cutoff of rho, return 0 */
  if (rhosq >= cut_rhosq) {
    *d_drhosq = 0;
    for (dim=0; dim<DIM; dim++) {
      d_dri[dim] = 0;
      d_drj[dim] = 0;
      d_drk1[dim] = 0;
      d_drk2[dim] = 0;
      d_drk3[dim] = 0;
      d_drl1[dim] = 0;
      d_drl2[dim] = 0;
      d_drl3[dim] = 0;
    }
    return 0;
  }


  /* get corresponding parameters */
  inter_idx = param_index(buffer, i, j);
  B = buffer->B[inter_idx];
  eta = buffer->eta[inter_idx];

  /* unpack data from buffer */
  nearest3neigh = buffer->nearest3neigh;
  coords = buffer->coords;


  /* 3 neighs of atoms i and j */
  for (m=0; m<3; m++) {
    k[m] = nearest3neigh[i][m];
    l[m] = nearest3neigh[j][m];
  }

  /* cos_omega_kijl and the derivatives w.r.t coords */
  for (m=0; m<3; m++) {
    for (n=0; n<3; n++) {
      cos_kl[m][n] = deriv_cos_omega(coords[k[m]], coords[i], coords[j], coords[l[n]],
          dcos_kl[m][n][0], dcos_kl[m][n][1], dcos_kl[m][n][2], dcos_kl[m][n][3]);
    }
  }

  epart1 = exp(eta * cos_kl[0][0]*cos_kl[0][1]*cos_kl[0][2]);
  epart2 = exp(eta * cos_kl[1][0]*cos_kl[1][1]*cos_kl[1][2]);
  epart3 = exp(eta * cos_kl[2][0]*cos_kl[2][1]*cos_kl[2][2]);
  D2 = epart1 + epart2 + epart3;

  /* cutoff function */
  D0 = tap_rho(rhosq, cut_rhosq, &d_drhosq_tap);
  D0 *= B;

  /* dihedral energy */
  dihe = D0*D2;

  /* deriv of dihedral w.r.t rhosq */
  *d_drhosq = B*d_drhosq_tap*D2;


  /* deriv of dihedral w.r.t cos_omega_kijl */
  d_dcos_kl[0][0] = D0* epart1 *eta *cos_kl[0][1]*cos_kl[0][2];
  d_dcos_kl[0][1] = D0* epart1 *eta *cos_kl[0][0]*cos_kl[0][2];
  d_dcos_kl[0][2] = D0* epart1 *eta *cos_kl[0][0]*cos_kl[0][1];
  d_dcos_kl[1][0] = D0* epart2 *eta *cos_kl[1][1]*cos_kl[1][2];
  d_dcos_kl[1][1] = D0* epart2 *eta *cos_kl[1][0]*cos_kl[1][2];
  d_dcos_kl[1][2] = D0* epart2 *eta *cos_kl[1][0]*cos_kl[1][1];
  d_dcos_kl[2][0] = D0* epart3 *eta *cos_kl[2][1]*cos_kl[2][2];
  d_dcos_kl[2][1] = D0* epart3 *eta *cos_kl[2][0]*cos_kl[2][2];
  d_dcos_kl[2][2] = D0* epart3 *eta *cos_kl[2][0]*cos_kl[2][1];

  /* initialization to be zero and later add values */
  for (dim=0; dim<DIM; dim++) {
    d_drk1[dim] = 0.;
    d_drk2[dim] = 0.;
    d_drk3[dim] = 0.;
    d_dri[dim] = 0.;
    d_drj[dim] = 0.;
    d_drl1[dim] = 0.;
    d_drl2[dim] = 0.;
    d_drl3[dim] = 0.;
  }

  for (m=0; m<3; m++) {
    for (dim=0; dim<3; dim++) {
      d_drk1[dim] += d_dcos_kl[0][m]*dcos_kl[0][m][0][dim];
      d_drk2[dim] += d_dcos_kl[1][m]*dcos_kl[1][m][0][dim];
      d_drk3[dim] += d_dcos_kl[2][m]*dcos_kl[2][m][0][dim];
      d_drl1[dim] += d_dcos_kl[m][0]*dcos_kl[m][0][3][dim];
      d_drl2[dim] += d_dcos_kl[m][1]*dcos_kl[m][1][3][dim];
      d_drl3[dim] += d_dcos_kl[m][2]*dcos_kl[m][2][3][dim];
    }
    for (n=0; n<3; n++) {
      for (dim=0; dim<3; dim++) {
        d_dri[dim] += d_dcos_kl[m][n]*dcos_kl[m][n][1][dim];
        d_drj[dim] += d_dcos_kl[m][n]*dcos_kl[m][n][2][dim];
      }
    }
  }

  return dihe;
}



/* compute cos(omega_kijl) and the derivateives */
static double deriv_cos_omega(double rk[DIM], double ri[DIM], double rj[DIM],
    double rl[DIM], double dcos_drk[DIM], double dcos_dri[DIM], double dcos_drj[DIM],
    double dcos_drl[DIM])

{

  double ejik[DIM];
  double eijl[DIM];
  double tmp1[DIM];
  double tmp2[DIM];
  double dejik_dri[DIM][DIM];
  double dejik_drj[DIM][DIM];
  double dejik_drk[DIM][DIM];
  double deijl_dri[DIM][DIM];
  double deijl_drj[DIM][DIM];
  double deijl_drl[DIM][DIM];

  double cos_omega;

  int m,n;


  /* ejik and derivatives (Note the dejik_dri ... returned are actually the transpose) */
  deriv_cross(ri, rj, rk, ejik, dejik_dri, dejik_drj, dejik_drk);
  /* flip sign because deriv_cross computes rij cross rik, here we need rji cross rik*/
  for (m=0; m<DIM; m++) {
    ejik[m] = - ejik[m];
    for (n=0; n<DIM; n++) {
      dejik_dri[m][n] = - dejik_dri[m][n];
      dejik_drj[m][n] = - dejik_drj[m][n];
      dejik_drk[m][n] = - dejik_drk[m][n];
    }
  }

  /* eijl and derivatives */
  deriv_cross(rj, ri, rl, eijl, deijl_drj, deijl_dri, deijl_drl);
  /* flip sign */
  for (m=0; m<DIM; m++) {
    eijl[m] = - eijl[m];
    for (n=0; n<DIM; n++) {
      deijl_drj[m][n] = - deijl_drj[m][n];
      deijl_dri[m][n] = - deijl_dri[m][n];
      deijl_drl[m][n] = - deijl_drl[m][n];
    }
  }

  /* dcos_drk */
  mat_dot_vec(dejik_drk, eijl, dcos_drk);
  /* dcos_dri */
  mat_dot_vec(dejik_dri, eijl, tmp1);
  mat_dot_vec(deijl_dri, ejik, tmp2);
  for (m=0; m<DIM; m++) {
    dcos_dri[m] = tmp1[m] + tmp2[m];
  }
  /* dcos_drj */
  mat_dot_vec(dejik_drj, eijl, tmp1);
  mat_dot_vec(deijl_drj, ejik, tmp2);
  for (m=0; m<DIM; m++) {
    dcos_drj[m] = tmp1[m] + tmp2[m];
  }
  /* dcos drl */
  mat_dot_vec(deijl_drl, ejik, dcos_drl);

  /* cos_oemga_kijl */
  cos_omega = dot(ejik, eijl);

  return cos_omega;
}


/* tap */
/*
static double tap(model_buffer* buffer, int i, int j, double r, double *const dtap)
{
  int inter_idx;
  double cutoff;
  double roc;
  double roc_sq;
  double t;

  inter_idx = param_index(buffer, i, j);
  cutoff = buffer->cutoff[inter_idx];

  roc = r/cutoff;
  roc_sq = roc*roc;
  t = roc_sq*roc_sq* (-35.0 + 84.0*roc + roc_sq* (-70.0 + 20.0*roc)) + 1;
  *dtap = roc_sq*roc/cutoff* (-140.0 + 420.0*roc + roc_sq* (-420.0 + 140.0*roc));

  return t;
}
*/

/* tap */
static double tap(model_buffer* buffer, int i, int j, double r, double *const dtap)
{
  int inter_idx;
  double cutoff;
  double roc;
  double roc_sq;
  double t;
  double r_min = 0;

  inter_idx = param_index(buffer, i, j);
  cutoff = buffer->cutoff[inter_idx];

  if (r <= r_min) {
    t = 1;
    *dtap = 0;
  } else {
    roc = (r-r_min)/(cutoff-r_min);
    roc_sq = roc*roc;
    t = roc_sq*roc_sq* (-35.0 + 84.0*roc + roc_sq* (-70.0 + 20.0*roc)) + 1;
    *dtap = roc_sq*roc/(cutoff-r_min)* (-140.0 + 420.0*roc + roc_sq* (-420.0 + 140.0*roc));
  }

  return t;
}




/* tap rho */
static double tap_rho(double rhosq, double cut_rhosq, double *const drhosq)
{
  double roc_sq;
  double roc;
  double t;

  roc_sq = rhosq/cut_rhosq;
  roc = sqrt(roc_sq);
  t = roc_sq*roc_sq* (-35.0 + 84.0*roc + roc_sq* (-70.0 + 20.0*roc)) + 1;
  *drhosq = roc_sq/cut_rhosq* (-70.0 + 210.0*roc + roc_sq* (-210.0 + 70.0*roc));

  return t;
}






/* helper functions */

/* parameters index */
static int param_index(model_buffer* const buffer, const int i, const int j)
{
  int iSpecies;
  int jSpecies;
  int index;

  iSpecies = buffer->particleSpecies[i];
  jSpecies = buffer->particleSpecies[j];
  if (iSpecies == SPEC1 && jSpecies == SPEC1) {
    index = 0;
  } else {
    printf("KIM error, param_index\n");
    exit(0);
  }

  return index;
}


/* relative position vector and squared distance */
static double rsq_rij(model_buffer* const buffer, const int i, const int j,
    double *const rij)
{
  double rsq;
  vectorDIM* coords = buffer->coords;

  rij[0] = coords[j][0] - coords[i][0];
  rij[1] = coords[j][1] - coords[i][1];
  rij[2] = coords[j][2] - coords[i][2];
  rsq = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];

  return rsq;
}


/* dot product of two vector */
double dot(const double x[DIM], const double y[DIM])
{
  return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}

/* matrix dot product a vector, return a vector */
static void mat_dot_vec(double X[DIM][DIM], const double y[DIM], double z[DIM])
{
  int k;
  for (k=0; k<DIM; k++) {
    z[k] = X[k][0]*y[0] + X[k][1]*y[1] + X[k][2]*y[2];
  }
}


/* Compute the normalized cross product of two vector rkl, rkm, and the derivates
   w.r.t rk, rl, rm.
  NOTE, the dcross_drk, dcross_drl, and dcross_drm is actually the transpose of the
  actual one.
*/
static void deriv_cross(double rk[DIM], double rl[DIM], double rm[DIM], double cross[DIM],
  double dcross_drk[DIM][DIM], double dcross_drl[DIM][DIM], double dcross_drm[DIM][DIM])
{

  double x[DIM];
  double y[DIM];
  double p[DIM];
  double q;
  double q_cubic;
  double d_invq_d_x0;
  double d_invq_d_x1;
  double d_invq_d_x2;
  double d_invq_d_y0;
  double d_invq_d_y1;
  double d_invq_d_y2;

  int i,j;


  /* get x = rkl and y = rkm */
  for (i=0; i<DIM; i++) {
    x[i] = rl[i] - rk[i];
    y[i] = rm[i] - rk[i];
  }

  /* cross product */
  p[0] = x[1]*y[2] - x[2]*y[1];
  p[1] = x[2]*y[0] - x[0]*y[2];
  p[2] = x[0]*y[1] - x[1]*y[0];

  q = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);

  /* normalized cross */
  cross[0] = p[0]/q;
  cross[1] = p[1]/q;
  cross[2] = p[2]/q;

  /* compute derivatives */
  /* derivative of inverse q (i.e. 1/q) w.r.t x and y*/
  q_cubic = q*q*q;
  d_invq_d_x0 = (           + p[1]*y[2] - p[2]*y[1])/q_cubic;
  d_invq_d_x1 = (-p[0]*y[2]             + p[2]*y[0])/q_cubic;
  d_invq_d_x2 = ( p[0]*y[1] - p[1]*y[0]            )/q_cubic;
  d_invq_d_y0 = (           - p[1]*x[2] + p[2]*x[1])/q_cubic;
  d_invq_d_y1 = ( p[0]*x[2]             - p[2]*x[0])/q_cubic;
  d_invq_d_y2 = (-p[0]*x[1] + p[1]*x[0]            )/q_cubic;

  /* dcross/drl transposed */
  dcross_drl[0][0] =           p[0]*d_invq_d_x0;
  dcross_drl[0][1] = -y[2]/q + p[1]*d_invq_d_x0;
  dcross_drl[0][2] =  y[1]/q + p[2]*d_invq_d_x0;

  dcross_drl[1][0] =  y[2]/q + p[0]*d_invq_d_x1;
  dcross_drl[1][1] =           p[1]*d_invq_d_x1;
  dcross_drl[1][2] = -y[0]/q + p[2]*d_invq_d_x1;

  dcross_drl[2][0] = -y[1]/q + p[0]*d_invq_d_x2;
  dcross_drl[2][1] =  y[0]/q + p[1]*d_invq_d_x2;
  dcross_drl[2][2] =           p[2]*d_invq_d_x2;

  /* dcross/drm transposed */
  dcross_drm[0][0] =           p[0]*d_invq_d_y0;
  dcross_drm[0][1] =  x[2]/q + p[1]*d_invq_d_y0;
  dcross_drm[0][2] = -x[1]/q + p[2]*d_invq_d_y0;

  dcross_drm[1][0] = -x[2]/q + p[0]*d_invq_d_y1;
  dcross_drm[1][1] =           p[1]*d_invq_d_y1;
  dcross_drm[1][2] =  x[0]/q + p[2]*d_invq_d_y1;

  dcross_drm[2][0] =  x[1]/q + p[0]*d_invq_d_y2;
  dcross_drm[2][1] = -x[0]/q + p[1]*d_invq_d_y2;
  dcross_drm[2][2] =           p[2]*d_invq_d_y2;

  /* dcross/drk transposed */
  for (i=0; i<DIM; i++) {
    for (j=0; j<DIM; j++) {
      dcross_drk[i][j] = - (dcross_drl[i][j] + dcross_drm[i][j]);
    }
  }

}


/* cross product, z = x cross y */
/*static void cross(const double x[DIM], const double y[DIM], double z[DIM])
{
  z[0] = x[1]*y[2] - x[2]*y[1];
  z[1] = x[2]*y[0] - x[0]*y[2];
  z[2] = x[0]*y[1] - x[1]*y[0];
}
*/


/* Initialization function */
int model_driver_init(void *km, char* paramfile_names, int* nmstrlen, int* numparamfiles)
{
  /* KIM variables */
  intptr_t* pkim = *((intptr_t**) km);
  char* paramfile1name;

  /* Local variables */
  FILE* fid;
  double* model_cutoff;
  double* model_C0;
  double* model_C2;
  double* model_C4;
  double* model_C;
  double* model_delta;
  double* model_lambda;
  double* model_A;
  double* model_z0;
  double* model_B;
  double* model_eta;

  model_buffer* buffer;
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
  model_B      = (double*) malloc(nInteract*sizeof(double));
  model_eta    = (double*) malloc(nInteract*sizeof(double));
  model_cutoff = (double*) malloc(nInteract*sizeof(double));
  if (model_C0==NULL
      || model_C2==NULL
      || model_C4==NULL
      || model_C==NULL
      || model_delta==NULL
      || model_lambda==NULL
      || model_A==NULL
      || model_z0==NULL
      || model_B==NULL
      || model_eta==NULL
      || model_cutoff==NULL) {
    ier = KIM_STATUS_FAIL;
    KIM_API_report_error(__LINE__, __FILE__, "malloc", ier);
    return ier;
  }

  /* read parameters */
  for (i=0; i< nInteract; ++i) {
    ier = fscanf(fid, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
        &model_C0[i],
        &model_C2[i],
        &model_C4[i],
        &model_C[i],
        &model_delta[i],
        &model_lambda[i],
        &model_A[i],
        &model_z0[i],
        &model_B[i],
        &model_eta[i],
        &model_cutoff[i]);
    /* check that we read the right number of parameters */
    if (11 != ier) {
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
    if ( model_cutoff[i] > *cutoff)
      *cutoff = model_cutoff[i];
  }

  /* store parameters in KIM object */
  KIM_API_setm_data(pkim, &ier, 11*4,
      "PARAM_FREE_cutoff",    nInteract,  model_cutoff,  1,
      "PARAM_FREE_C0",        nInteract,  model_C0,      1,
      "PARAM_FREE_C2",        nInteract,  model_C2,      1,
      "PARAM_FREE_C4",        nInteract,  model_C4,      1,
      "PARAM_FREE_C",         nInteract,  model_C,       1,
      "PARAM_FREE_delta",     nInteract,  model_delta,   1,
      "PARAM_FREE_lambda",    nInteract,  model_lambda,  1,
      "PARAM_FREE_A",         nInteract,  model_A,       1,
      "PARAM_FREE_z0",        nInteract,  model_z0,      1,
      "PARAM_FREE_B",         nInteract,  model_B,       1,
      "PARAM_FREE_eta",       nInteract,  model_eta,     1);
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_setm_data", ier);
    return ier;
  }

  /* allocate buffer */
  buffer = (model_buffer*) malloc(sizeof(model_buffer));
  if (NULL==buffer) {
    ier = KIM_STATUS_FAIL;
    KIM_API_report_error(__LINE__, __FILE__, "malloc", ier);
    return ier;
  }

  /* store parameters in buffer */
  buffer->model_index_shift = KIM_API_get_model_index_shift(pkim);
  buffer->cutoff = model_cutoff;
  buffer->C0     = model_C0;
  buffer->C2     = model_C2;
  buffer->C4     = model_C4;
  buffer->C      = model_C;
  buffer->delta  = model_delta;
  buffer->lambda = model_lambda;
  buffer->A      = model_A;
  buffer->z0     = model_z0;
  buffer->B     = model_B;
  buffer->eta     = model_eta;
  /* end setup buffer */

  /* store in model buffer */
  KIM_API_set_model_buffer(pkim, (void*) buffer, &ier);
  if (KIM_STATUS_OK > ier) {
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
  model_buffer* buffer;
  int* nSpecies;
  int nInteract;
  int ier;
  int i;

  /* get buffer from KIM object */
  buffer = (model_buffer*) KIM_API_get_model_buffer(pkim, &ier);
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_model_buffer", ier);
    return ier;
  }

  nSpecies = KIM_API_get_data(pkim, "numberOfSpecies", &ier);
  nInteract = (*nSpecies + 1)*(*nSpecies)/2;

  /* update cutoff in KIM API */
  cutoff = (double*) KIM_API_get_data(pkim, "cutoff", &ier);
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_data", ier);
    return ier;
  }
  *cutoff = 0.0;
  for (i=0; i< nInteract; ++i) {
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
  model_buffer* buffer;
  int ier;

  /* get model buffer from KIM object */
  buffer = (model_buffer*) KIM_API_get_model_buffer(pkim, &ier);
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_model_buffer", ier);
    return ier;
  }

  /* free parameters */
  free(buffer->cutoff);
  free(buffer->C0);
  free(buffer->C2);
  free(buffer->C4);
  free(buffer->C);
  free(buffer->delta);
  free(buffer->lambda);
  free(buffer->A);
  free(buffer->z0);
  free(buffer->B);
  free(buffer->eta);

  /* destroy the buffer */
  free(buffer);

  ier = KIM_STATUS_OK;
  return ier;
}
