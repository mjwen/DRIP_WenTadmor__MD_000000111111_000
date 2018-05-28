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

#define DIMENSION 3
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
typedef double VectorOfSizeDIM[DIMENSION];


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
  // TODO add functions declaration of your implementation

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
    const int* const particleSpecies,
    const int* const particleContributing,
    const VectorOfSizeDIM* const coordinates,
    double* const energy,
    VectorOfSizeDIM* const forces,
    double* const particleEnergy)
{
  int ier = false;

  if ((isComputeEnergy == false) &&
      (isComputeParticleEnergy == false) &&
      (isComputeForces == false) &&
      (isComputeProcess_dEdr == false) &&
      (isComputeProcess_d2Edr2 == false))
    return ier;

  // initialize energy and forces
  if (isComputeEnergy == true) {
    *energy = 0.0;
  }

  if (isComputeParticleEnergy == true) {
    for (int i = 0; i < cachedNumberOfParticles_; ++i) {
      particleEnergy[i] = 0.0;
    }
  }

  if (isComputeForces == true) {
    for (int i = 0; i < cachedNumberOfParticles_; ++i) {
      for (int j = 0; j < DIMENSION; ++j)
        forces[i][j] = 0.0;
    }
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

      // Setup loop over neighbors of current particle
      for (int jj = 0; jj < numnei; ++jj) {
        int const j = n1atom[jj];
        int const jSpecies = particleSpecies[j];

        // Compute rij
        double rij[DIMENSION];
        for (int dim = 0; dim < DIMENSION; ++dim)
          rij[dim] = coordinates[j][dim] - coordinates[i][dim];

        // compute distance squared
        double const rij_sq = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
        double const rij_mag = sqrt(rij_sq);

        if (rij_sq <= cutoffSq_2D_[iSpecies][jSpecies]) {

          // two-body contributions
          // TODO add your implementation




        }  // if particleContributing
      }  // if particles i and j interact
    }  // end of first neighbor loop
  }  // end of loop over contributing particles

  // everything is good
  ier = false;
  return ier;
}

#endif  // RDP_IMPLEMENTATION_HPP_
