//
// CDDL HEADER START
//
// The contents of this file are subject to the terms of the Common Development
// and Distribution License Version 1.0 (the "License").
//
// You can obtain a copy of the license at
// http://www.opensource.org/licenses/CDDL-1.0.  See the License for the
// specific language governing permissions and limitations under the License.
//
// When distributing Covered Code, include this CDDL HEADER in each file and
// include the License file in a prominent location with the name LICENSE.CDDL.
// If applicable, add the following below this CDDL HEADER, with the fields
// enclosed by brackets "[]" replaced with your own identifying information:
//
// Portions Copyright (c) [yyyy] [name of copyright owner]. All rights reserved.
//
// CDDL HEADER END
//

//
// Copyright (c) 2018, Regents of the University of Minnesota.
// All rights reserved.
//
// Contributors:
//    Mingjian Wen


#ifndef RDP_IMPLEMENTATION_HPP_
#define RDP_IMPLEMENTATION_HPP_

#include <vector>
#include "KIM_LogVerbosity.hpp"
#include "RDP.hpp"
#include "helper.hpp"

#define DIM 3
#define ONE 1.0
#define HALF 0.5

#define MAX_PARAMETER_FILES 1


//==============================================================================
//
// Type definitions, enumerations, and helper function prototypes
//
//==============================================================================

// type declaration for get neighbor functions
typedef int (GetNeighborFunction)(void const * const, int const,
                                  int * const, int const ** const);

// type declaration for vector of constant dimension
typedef double VectorOfSizeDIM[DIM];
typedef int VectorOfSize3Int[DIM];


//==============================================================================
//
// Declaration of RDPImplementation class
//
//==============================================================================

//******************************************************************************
class RDPImplementation
{
 public:
  RDPImplementation(
      KIM::ModelDriverCreate * const modelDriverCreate,
      KIM::LengthUnit const requestedLengthUnit,
      KIM::EnergyUnit const requestedEnergyUnit,
      KIM::ChargeUnit const requestedChargeUnit,
      KIM::TemperatureUnit const requestedTemperatureUnit,
      KIM::TimeUnit const requestedTimeUnit,
      int * const ier);
  ~RDPImplementation();  // no explicit Destroy() needed here

  int Refresh(KIM::ModelRefresh * const modelRefresh);
  int Compute(KIM::ModelCompute const * const modelCompute,
      KIM::ModelComputeArguments const * const modelComputeArguments);
  int ComputeArgumentsCreate(
      KIM::ModelComputeArgumentsCreate * const modelComputeArgumentsCreate)
      const;
  int ComputeArgumentsDestroy(
      KIM::ModelComputeArgumentsDestroy * const modelComputeArgumentsDestroy)
      const;


 private:
  // Constant values that never change
  //   Set in constructor (via SetConstantValues)
  //
  //
  // RDPImplementation: constants
  int numberModelSpecies_;
  std::vector<int> modelSpeciesCodeList_;
  int numberUniqueSpeciesPairs_;


  // Constant values that are read from the input files and never change
  //   Set in constructor (via functions listed below)
  //
  //
  // KIM API: Model Fixed Parameters
  //   Memory allocated in   AllocateFixedParameterMemory()
  //   Memory deallocated in destructor
  //   Data set in ReadParameterFile routines
  // none
  //
  // KIM API: Model Free Parameters whose (pointer) values never change
  //   Memory allocated in   AllocateFreeParameterMemory() (from constructor)
  //   Memory deallocated in destructor
  //   Data set in ReadParameterFile routines OR by KIM Simulator
  double* cutoff_;
  double* rhocut_;
  double* C0_;
  double* C2_;
  double* C4_;
  double* C_;
  double* delta_;
  double* lambda_;
  double* A_;
  double* z0_;
  double* B_;
  double* eta_;


  // Mutable values that only change when Refresh() executes
  //   Set in Refresh (via SetRefreshMutableValues)
  //
  //
  // KIM API: Model Fixed Parameters
  // none
  //
  // KIM API: Model Free Parameters
  // none
  //
  // RDPImplementation: values
  double influenceDistance_;
  double** cutoffSq_2D_;
  double** rhocutSq_2D_;
  double** C0_2D_;
  double** C2_2D_;
  double** C4_2D_;
  double** C_2D_;
  double** delta_2D_;
  double** lambda_2D_;
  double** A_2D_;
  double** z0_2D_;
  double** B_2D_;
  double** eta_2D_;


  // Mutable values that can change with each call to Refresh() and Compute()
  //   Memory may be reallocated on each call
  //
  //
  // RDPImplementation: values that change
  int cachedNumberOfParticles_;
  int cachedIsComputeEnergy_;
  int cachedIsComputeParticleEnergy_;
  int cachedIsComputeForces_;
  int cachedIsComputeProcess_dEdr_;
  int cachedIsComputeProcess_d2Edr2_;


  // Helper methods
  //
  //
  // Related to constructor
  void AllocateFreeParameterMemory();
  static int OpenParameterFiles(
      KIM::ModelDriverCreate * const modelDriverCreate,
      int const numberParameterFiles,
      FILE* parameterFilePointers[MAX_PARAMETER_FILES]);
  static void CloseParameterFiles(
      int const numberParameterFiles,
      FILE* const parameterFilePointers[MAX_PARAMETER_FILES]);
  int ProcessParameterFiles(
      KIM::ModelDriverCreate * const modelDriverCreate,
      int const numberParameterFiles,
      FILE* const parameterFilePointers[MAX_PARAMETER_FILES]);
  void getNextDataLine(FILE* const filePtr, char* const nextLine,
                       int const maxSize, int* endOfFileFlag);
  int ConvertUnits(
      KIM::ModelDriverCreate * const modelDriverCreate,
      KIM::LengthUnit const requestedLengthUnit,
      KIM::EnergyUnit const requestedEnergyUnit,
      KIM::ChargeUnit const requestedChargeUnit,
      KIM::TemperatureUnit const requestedTemperatureUnit,
      KIM::TimeUnit const requestedTimeUnit);
  int RegisterKIMModelSettings(KIM::ModelDriverCreate * const modelDriverCreate) const;
  int RegisterKIMComputeArgumentsSettings(
      KIM::ModelComputeArgumentsCreate * const modelComputeArgumentsCreate) const;
  int RegisterKIMParameters(KIM::ModelDriverCreate * const modelDriverCreate);
  int RegisterKIMFunctions(KIM::ModelDriverCreate * const modelDriverCreate) const;
  //
  // Related to Refresh()
  template<class ModelObj>
  int SetRefreshMutableValues(ModelObj * const modelObj);

  //
  // Related to Compute()
  int SetComputeMutableValues(
      KIM::ModelComputeArguments const * const modelCompute,
      bool& isComputeProcess_dEdr,
      bool& isComputeProcess_d2Edr2,
      bool& isComputeEnergy,
      bool& isComputeForces,
      bool& isComputeParticleEnergy,
      int const*& particleSpecies,
      int const*& particleContributing,
      VectorOfSizeDIM const*& coordinates,
      double*& energy,
      double*& particleEnergy,
      VectorOfSizeDIM*& forces);
  int CheckParticleSpecies(KIM::ModelCompute const * const modelCompute,
                           int const* const particleSpecies) const;
  int GetComputeIndex(const bool& isComputeProcess_dEdr,
                      const bool& isComputeProcess_d2Edr2,
                      const bool& isComputeEnergy,
                      const bool& isComputeForces,
                      const bool& isComputeParticleEnergy) const;

  // compute functions
  template< bool isComputeProcess_dEdr, bool isComputeProcess_d2Edr2,
            bool isComputeEnergy, bool isComputeForces,
            bool isComputeParticleEnergy>
  int Compute(KIM::ModelCompute const * const modelCompute,
              KIM::ModelComputeArguments const * const modelComputeArguments,
              const int* const particleSpecies,
              const int* const particleContributing,
              const VectorOfSizeDIM* const coordinates,
              double* const energy,
              VectorOfSizeDIM* const forces,
              double* const particleEnergy);


  // RDP functions
 int create_layers(KIM::ModelCompute const * const modelCompute,
    KIM::ModelComputeArguments const * const modelComputeArguments,
    VectorOfSizeDIM const * const coordinates,
    int * const in_layer,
    VectorOfSize3Int * const nearest3neigh);

  double calc_attractive(int const i, int const j,
     int const iSpecies, int const jSpecies, double const * const rij,
     double const r, VectorOfSizeDIM * const forces);

  double calc_repulsive(int const i, int const j,
      int const * const particleSpecies,
      VectorOfSizeDIM const * const coordinates,
      VectorOfSize3Int const * const nearest3neigh,
      const double * const rij,
      double const r, const int nbj1, const int nbj2, const int nbj3,
      double const * const ni,  VectorOfSizeDIM const * const dni_dri,
      VectorOfSizeDIM const * const dni_drnb1, VectorOfSizeDIM const * const dni_drnb2,
      VectorOfSizeDIM const * const dni_drnb3, VectorOfSizeDIM * const forces);

  void normal(const int i, VectorOfSizeDIM const * const coordinates,
      VectorOfSize3Int const * const nearest3neigh,
      int & k1, int & k2, int & k3, double * const normal,
      VectorOfSizeDIM * const dn_dri, VectorOfSizeDIM * const dn_drk1,
      VectorOfSizeDIM * const dn_drk2, VectorOfSizeDIM * const dn_drk3);

  double td(double C0, double C2, double C4, double delta,
      double const * const rvec, double r, const double* const n, double & rho_sq,
      double & dtd);

  void get_drhosqij(double const * const rij,
      double const * const ni,
      VectorOfSizeDIM const * const dni_dri,
      VectorOfSizeDIM const * const dni_drn1,
      VectorOfSizeDIM const * const dni_drn2,
      VectorOfSizeDIM const * const dni_drn3,
      double * const drhosq_dri, double * const drhosq_drj,
      double * const drhosq_drn1, double * const drhosq_drn2,
      double * const drhosq_drn3);

  double dihedral(const int i, const int j,
      int const * const particleSpecies,
      VectorOfSizeDIM const * const coordinates,
      VectorOfSize3Int const * const nearest3neigh,
      double const rhosq,
      double & d_drhosq, double d_dri[DIM], double d_drj[DIM],
      double d_drk1[DIM], double d_drk2[DIM], double d_drk3[DIM],
      double d_drl1[DIM], double d_drl2[DIM], double d_drl3[DIM]);

  double deriv_cos_omega(double const * const rk,
      double const * const ri, double const * const rj, double const * const rl,
      double * const dcos_drk, double * const dcos_dri, double * const dcos_drj,
      double * const dcos_drl);

  double tap(double r, double cutoff, double & dtap);

  double tap_rho(double const rhosq, double cut_rhosq, double & drhosq);

  // helper
  double dot(double const * const x, double const * const y);

  void deriv_cross(double const * const rk, double const * const rl,
      double const * const rm, double * const cross, VectorOfSizeDIM * const dcross_drk,
      VectorOfSizeDIM * const dcross_drl, VectorOfSizeDIM * const dcross_drm);

  void mat_dot_vec(VectorOfSizeDIM const * const X, double const * const y,
        double * const z);

};

//==============================================================================
//
// Definition of RDPImplementation::Compute functions
//
// NOTE: Here we rely on the compiler optimizations to prune dead code
//       after the template expansions.  This provides high efficiency
//       and easy maintenance.
//
//==============================================================================

#include "KIM_ModelComputeLogMacros.hpp"
template< bool isComputeProcess_dEdr, bool isComputeProcess_d2Edr2,
          bool isComputeEnergy, bool isComputeForces,
          bool isComputeParticleEnergy>
int RDPImplementation::Compute(
    KIM::ModelCompute const * const modelCompute,
    KIM::ModelComputeArguments const * const modelComputeArguments,
    int const * const particleSpecies,
    int const * const particleContributing,
    VectorOfSizeDIM const * const coordinates,
    double * const energy,
    VectorOfSizeDIM * const forces,
    double * const particleEnergy)
{
  int ier = false;
  const int Natoms = cachedNumberOfParticles_;

  if ((isComputeEnergy == false) &&
      (isComputeParticleEnergy == false) &&
      (isComputeForces == false) &&
      (isComputeProcess_dEdr == false) &&
      (isComputeProcess_d2Edr2 == false))
    return ier;

  // put them global
  cachedIsComputeEnergy_ = isComputeEnergy;
  cachedIsComputeParticleEnergy_ = isComputeParticleEnergy;
  cachedIsComputeForces_ = isComputeForces;
  cachedIsComputeProcess_dEdr_ = isComputeProcess_dEdr;
  cachedIsComputeProcess_d2Edr2_ = isComputeProcess_d2Edr2;

  // initialize energy and forces
  if (isComputeEnergy == true) {
    *energy = 0.0;
  }

  if (isComputeParticleEnergy == true) {
    for (int i = 0; i < Natoms; ++i) {
      particleEnergy[i] = 0.0;
    }
  }

  if (isComputeForces == true) {
    for (int i = 0; i < Natoms; ++i) {
      for (int j = 0; j < DIM; ++j)
        forces[i][j] = 0.0;
    }
  }


  // allocate atoms into layers
  int * in_layer; // atoms in which layer , -1: not classified in any layer
  int ** n3n; // nearest 3 neighbors of atom
  AllocateAndInitialize1DArray<int>(in_layer, Natoms);
  AllocateAndInitialize2DArray<int>(n3n, Natoms, 3);
  VectorOfSize3Int * const nearest3neigh = (VectorOfSize3Int*) n3n[0];

  ier = create_layers(modelCompute, modelComputeArguments, coordinates,
      in_layer, nearest3neigh);
  if (ier) {
    LOG_ERROR("create_layers");
    return ier;
  }


  // calculate contribution from pair function
  //
  // Setup loop over contributing particles
  int i = 0;
  int numnei = 0;
  int const * n1atom = 0;

  for (i = 0; i < cachedNumberOfParticles_; ++i) {
    if (particleContributing[i]) {
      modelComputeArguments->GetNeighborList(0, i, &numnei, &n1atom);
      int const iSpecies = particleSpecies[i];
      int const ilayer = in_layer[i];

      // three nearest neighbors of atom i
      int nbi1;
      int nbi2;
      int nbi3;
      VectorOfSizeDIM ni; // normal at atom i
      VectorOfSizeDIM dni_dri[DIM];
      VectorOfSizeDIM dni_drnb1[DIM];
      VectorOfSizeDIM dni_drnb2[DIM];
      VectorOfSizeDIM dni_drnb3[DIM];

      // get 3 nearest neighs of atom i and compute the derivative of ni w.r.t. them
      normal(i, coordinates, nearest3neigh, nbi1, nbi2, nbi3, ni, dni_dri,
          dni_drnb1, dni_drnb2, dni_drnb3);

      // Setup loop over neighbors of current particle
      for (int jj = 0; jj < numnei; ++jj) {
        int const j = n1atom[jj];
        int const jSpecies = particleSpecies[j];
        int const jlayer = in_layer[j];
        if (ilayer == jlayer) continue;

        // Compute rij
        double rij[DIM];
        for (int dim = 0; dim < DIM; ++dim)
          rij[dim] = coordinates[j][dim] - coordinates[i][dim];

        // compute distance squared
        double const rij_sq = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
        double const rij_mag = sqrt(rij_sq);

        if (rij_sq <= cutoffSq_2D_[iSpecies][jSpecies]) {

          double phi_attr = calc_attractive(i, j, iSpecies, jSpecies, rij,
              rij_mag, forces);
          double phi_repul = calc_repulsive(i, j, particleSpecies, coordinates,
              nearest3neigh, rij, rij_mag, nbi1, nbi2, nbi3, ni,
              dni_dri, dni_drnb1, dni_drnb2, dni_drnb3, forces);

          if (isComputeEnergy) {
            // HALF because of full neighborlist
            *energy += HALF * (phi_repul + phi_attr);
          }
          if (isComputeParticleEnergy) {
            for (int k=0; k<Natoms; k++) {
              particleEnergy[k] +=  HALF * (phi_repul + phi_attr);
            }
          }


        }  // if particleContributing
      }  // if particles i and j interact
    }  // end of first neighbor loop
  }  // end of loop over contributing particles


  // deallocate memory
  Deallocate1DArray<int> (in_layer);
  Deallocate2DArray<int> (n3n);

  // everything is good
  ier = false;
  return ier;
}

#endif  // RDP_IMPLEMENTATION_HPP_
