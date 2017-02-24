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
    const double* rvec, double r, const double* n, double *const dtd);

static void get_drhosqij(const double rij[DIM], const double ni[DIM],
    double dni_dri[DIM][DIM], double dni_drn1[DIM][DIM],
    double dni_drn2[DIM][DIM], double dni_drn3[DIM][DIM],
    double drhosq_dri[DIM], double drhosq_drj[DIM],
    double drhosq_drn1[DIM], double drhosq_drn2[DIM],
    double drhosq_drn3[DIM]);

/* helper */
static double rsq_rij(model_buffer* const buffer, const int i, const int j,
    double *const rij);

static int param_index(model_buffer* const buffer, const int i, const int j);

double dot(const double x[DIM], const double y[DIM]);

static void cross(const double x[DIM], const double y[DIM], double rslt[DIM]);

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


  /* unpack data from buffer */
  buffer = (model_buffer*) KIM_API_get_model_buffer(pkim, &ier);
  if (KIM_STATUS_OK > ier) {
    KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_model_buffer", ier);
    return ier;
  }
  nAtoms = *buffer->nAtoms;
  model_index_shift = buffer->model_index_shift;
  get_neigh = buffer->get_neigh;

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
  phi = - A/(roz0_sq*roz0_sq*roz0_sq);

  /* forces */
  if (buffer->comp_forces) {
    fpair = - HALF * 6*phi/r;
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

  double V1;
  double V2;
  double tdij, tdji;
  double dtdij, dtdji;
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

  /* derivative of rho w.r.t coords of atoms i, j, and the nearests 3 neighs of i */
  get_drhosqij(rij, ni, dni_dri, dni_drnb1, dni_drnb2, dni_drnb3, drhosqij_dri,
      drhosqij_drj, drhosqij_drnb1, drhosqij_drnb2, drhosqij_drnb3);
  get_drhosqij(rji, nj, dnj_drj, dnj_drnb1, dnj_drnb2, dnj_drnb3, drhosqji_drj,
      drhosqji_dri, drhosqji_drnb1, drhosqji_drnb2, drhosqji_drnb3);

  /* transverse decay function f(rho) and its derivative w.r.t. rhosq */
  tdij = td(C0, C2, C4, delta, rij, r, ni, &dtdij);
  tdji = td(C0, C2, C4, delta, rji, r, nj, &dtdji);

  V1 = exp(-lambda*(r-z0));
  V2 = C + tdij + tdji;

  phi = V1*V2;

  /* forces */
  if (buffer->comp_forces) {
    for (k=0; k<DIM; k++) {

      /* derivative of V1 */
      tmp = - HALF * lambda * phi * rij[k]/r;
      forces[i][k] += tmp;
      forces[j][k] -= tmp;

      /* derivative of V2 contribute to atoms i, j */
      forces[i][k] -= HALF * V1 * (dtdij*drhosqij_dri[k] + dtdji*drhosqji_dri[k]);
      forces[j][k] += HALF * V1 * (dtdij*drhosqij_drj[k] + dtdji*drhosqji_drj[k]);

      /* derivative of V2 contribute to neighs of atom i */
      forces[nbi1][k] -= HALF * V1 * dtdij*drhosqij_drnb1[k];
      forces[nbi2][k] -= HALF * V1 * dtdij*drhosqij_drnb2[k];
      forces[nbi3][k] -= HALF * V1 * dtdij*drhosqij_drnb3[k];

      /* derivative of V2 contribute to neighs of atom j */
      forces[nbj1][k] -= HALF * V1 * dtdji*drhosqji_drnb1[k];
      forces[nbj2][k] -= HALF * V1 * dtdji*drhosqji_drnb2[k];
      forces[nbj3][k] -= HALF * V1 * dtdji*drhosqji_drnb3[k];
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


/* compute normal and its derivatives  w.r.t atom ri, and its 3 nearest neighs k, l, m.*/
static void normal(model_buffer *const buffer, const int i, int *const nbi1,
    int* const nbi2, int* const nbi3, double normal[DIM], double dn_dri[DIM][DIM],
    double dn_drk[DIM][DIM], double dn_drl[DIM][DIM], double dn_drm[DIM][DIM])
{
  /* local variables */
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
  int (*nearest3neigh)[3];
  vectorDIM* coords;
  int j,k;

  /* unpack data from buffer */
  nearest3neigh = buffer->nearest3neigh;
  coords = buffer->coords;

  *nbi1 = nearest3neigh[i][0];
  *nbi2 = nearest3neigh[i][1];
  *nbi3 = nearest3neigh[i][2];

  /* get x = rkl and y = rkm */
  for (j=0; j<DIM; j++) {
    x[j] = coords[*nbi2][j] - coords[*nbi1][j];
    y[j] = coords[*nbi3][j] - coords[*nbi1][j];
    for (k=0; k<DIM; k++) {
      /*dn/dri transposed */
      dn_dri[j][k] = 0.0; /* normal does not depend on i*/
    }
  }

  /* the following is standard method to compute the derivatives of x cross y */

  /* cross product */
  cross(x, y, p);
  q = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);

  /* normal */
  normal[0] = p[0]/q;
  normal[1] = p[1]/q;
  normal[2] = p[2]/q;

  /* compute derivatives */
  /* derivative of inverse q (i.e. 1/q) w.r.t x and y*/
  q_cubic = q*q*q;
  d_invq_d_x0 = (           + p[1]*y[2] - p[2]*y[1])/q_cubic;
  d_invq_d_x1 = (-p[0]*y[2]             + p[2]*y[0])/q_cubic;
  d_invq_d_x2 = ( p[0]*y[1] - p[1]*y[0]            )/q_cubic;
  d_invq_d_y0 = (           - p[1]*x[2] + p[2]*x[1])/q_cubic;
  d_invq_d_y1 = ( p[0]*x[2]             - p[2]*x[0])/q_cubic;
  d_invq_d_y2 = (-p[0]*x[1] + p[1]*x[0]            )/q_cubic;

  /* dn/drl transposed */
  dn_drl[0][0] =           p[0]*d_invq_d_x0;
  dn_drl[0][1] = -y[2]/q + p[1]*d_invq_d_x0;
  dn_drl[0][2] =  y[1]/q + p[2]*d_invq_d_x0;

  dn_drl[1][0] =  y[2]/q + p[0]*d_invq_d_x1;
  dn_drl[1][1] =           p[1]*d_invq_d_x1;
  dn_drl[1][2] = -y[0]/q + p[2]*d_invq_d_x1;

  dn_drl[2][0] = -y[1]/q + p[0]*d_invq_d_x2;
  dn_drl[2][1] =  y[0]/q + p[1]*d_invq_d_x2;
  dn_drl[2][2] =           p[2]*d_invq_d_x2;

  /* dn/drm transposed */
  dn_drm[0][0] =           p[0]*d_invq_d_y0;
  dn_drm[0][1] =  x[2]/q + p[1]*d_invq_d_y0;
  dn_drm[0][2] = -x[1]/q + p[2]*d_invq_d_y0;

  dn_drm[1][0] = -x[2]/q + p[0]*d_invq_d_y1;
  dn_drm[1][1] =           p[1]*d_invq_d_y1;
  dn_drm[1][2] =  x[0]/q + p[2]*d_invq_d_y1;

  dn_drm[2][0] =  x[1]/q + p[0]*d_invq_d_y2;
  dn_drm[2][1] = -x[0]/q + p[1]*d_invq_d_y2;
  dn_drm[2][2] =           p[2]*d_invq_d_y2;

  /* dn/drk transposed */
  for (j=0; j<DIM; j++) {
    for (k=0; k<DIM; k++) {
      dn_drk[j][k] = - (dn_drl[j][k] + dn_drm[j][k]);
    }
  }

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
    const double* rvec, double r, const double* n, double *const dtd)
{
  double n_dot_r;
  double rho_sq;
  double del_sq;
  double rod_sq;
  double td;

  n_dot_r = dot(n, rvec);
  rho_sq = r*r - n_dot_r*n_dot_r;
  del_sq = delta*delta;
  rod_sq = rho_sq/del_sq;
  td = exp(-rod_sq) * (C0 + (rod_sq*(C2 + rod_sq*C4)));
  *dtd = td/del_sq + exp(-rod_sq) * (C2 + 2*C4*rod_sq)/del_sq;

  return td;
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

/* cross product, z = x cross y */
static void cross(const double x[DIM], const double y[DIM], double z[DIM])
{
  z[0] = x[1]*y[2] - x[2]*y[1];
  z[1] = x[2]*y[0] - x[0]*y[2];
  z[2] = x[0]*y[1] - x[1]*y[0];
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
  double* model_C0;
  double* model_C2;
  double* model_C4;
  double* model_C;
  double* model_delta;
  double* model_lambda;
  double* model_A;
  double* model_z0;

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
  model_cutoff = (double*) malloc(nInteract*sizeof(double));
  if (model_C0==NULL
      || model_C2==NULL
      || model_C4==NULL
      || model_C==NULL
      || model_delta==NULL
      || model_lambda==NULL
      || model_A==NULL
      || model_z0==NULL
      || model_cutoff==NULL) {
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

  /* destroy the buffer */
  free(buffer);

  ier = KIM_STATUS_OK;
  return ier;
}
