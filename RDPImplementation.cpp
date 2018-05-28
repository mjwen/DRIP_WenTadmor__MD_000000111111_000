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
// Copyright (c) 2013--2015, Regents of the University of Minnesota.
// All rights reserved.
//
// Contributors:
//    Mingjian Wen
//    Ryan S. Elliott
//    Stephen M. Whalen
//    Andrew Akerson
//


#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>

#include "RDPImplementation.hpp"
#include "KIM_Numbering.hpp"
#include "KIM_LanguageName.hpp"
#include "KIM_SpeciesName.hpp"
#include "KIM_SupportStatus.hpp"
#include "KIM_ComputeArgumentName.hpp"
#include "KIM_ComputeCallbackName.hpp"

#define MAXLINE 1024


//==============================================================================
//
// Implementation of RDPImplementation public member functions
//
//==============================================================================

//******************************************************************************
#include "KIM_ModelDriverCreateLogMacros.hpp"
RDPImplementation::RDPImplementation(
    KIM::ModelDriverCreate * const modelDriverCreate,
    KIM::LengthUnit const requestedLengthUnit,
    KIM::EnergyUnit const requestedEnergyUnit,
    KIM::ChargeUnit const requestedChargeUnit,
    KIM::TemperatureUnit const requestedTemperatureUnit,
    KIM::TimeUnit const requestedTimeUnit,
    int * const ier)
: numberModelSpecies_(0),
  numberUniqueSpeciesPairs_(0),
  cutoff_(0),
  // TODO add potential parameters
  A_(0),
  influenceDistance_(0.0),
  cutoffSq_2D_(0),
  // TODO add potential parameters in 2D
  A_2D_(0),
  cachedNumberOfParticles_(0)
{
  FILE* parameterFilePointers[MAX_PARAMETER_FILES];
  int numberParameterFiles;
  modelDriverCreate->GetNumberOfParameterFiles(&numberParameterFiles);
  *ier = OpenParameterFiles(modelDriverCreate, numberParameterFiles,
                            parameterFilePointers);
  if (*ier) return;

  *ier = ProcessParameterFiles(modelDriverCreate, numberParameterFiles,
                               parameterFilePointers);
  CloseParameterFiles(numberParameterFiles, parameterFilePointers);
  if (*ier) return;

  *ier = ConvertUnits(modelDriverCreate,
                      requestedLengthUnit,
                      requestedEnergyUnit,
                      requestedChargeUnit,
                      requestedTemperatureUnit,
                      requestedTimeUnit);
  if (*ier) return;

  *ier = SetRefreshMutableValues(modelDriverCreate);
  if (*ier) return;

  *ier = RegisterKIMModelSettings(modelDriverCreate);
  if (*ier) return;

  *ier = RegisterKIMParameters(modelDriverCreate);
  if (*ier) return;

  *ier = RegisterKIMFunctions(modelDriverCreate);
  if (*ier) return;

  // everything is good
  *ier = false;
  return;
}

//******************************************************************************
RDPImplementation::~RDPImplementation()
{ // note: it is ok to delete a null pointer and we have ensured that
  // everything is initialized to null

  Deallocate1DArray(cutoff_);
  //TODO deallocate memory of your parameters
  Deallocate1DArray(A_);

  Deallocate2DArray(cutoffSq_2D_);
  //TODO deallocate memory of your parameters
  Deallocate2DArray(A_2D_);
}

//******************************************************************************
int RDPImplementation::Refresh(
    KIM::ModelRefresh * const modelRefresh)
{
  int ier;

  ier = SetRefreshMutableValues(modelRefresh);
  if (ier) return ier;

  //TODO you may want to recompute cutoffSq_2D_ and others related to parameters


  // everything is good
  ier = false;
  return ier;
}

//******************************************************************************
int RDPImplementation::Compute(
    KIM::ModelCompute const * const modelCompute,
    KIM::ModelComputeArguments const * const modelComputeArguments)
{
  int ier;

  // KIM API Model Input compute flags
  bool isComputeProcess_dEdr = false;
  bool isComputeProcess_d2Edr2 = false;
  //
  // KIM API Model Output compute flags
  bool isComputeEnergy = false;
  bool isComputeForces = false;
  bool isComputeParticleEnergy = false;
  //
  // KIM API Model Input
  int const* particleSpecies = 0;
  int const* particleContributing = 0;
  VectorOfSizeDIM const* coordinates = 0;
  //
  // KIM API Model Output
  double* energy = 0;
  double* particleEnergy = 0;
  VectorOfSizeDIM* forces = 0;
  ier = SetComputeMutableValues(modelComputeArguments,
                                isComputeProcess_dEdr,
                                isComputeProcess_d2Edr2, isComputeEnergy,
                                isComputeForces, isComputeParticleEnergy,
                                particleSpecies, particleContributing,
                                coordinates, energy, particleEnergy, forces);
  if (ier) return ier;

  // Skip this check for efficiency
  //
  // ier = CheckParticleSpecies(pkim, particleSpecies);
  // if (ier) return ier;


#include "RDPImplementationComputeDispatch.cpp"
  return ier;
}



//******************************************************************************
int RDPImplementation::ComputeArgumentsCreate(
    KIM::ModelComputeArgumentsCreate * const modelComputeArgumentsCreate) const
{
  int ier;

  ier = RegisterKIMComputeArgumentsSettings(modelComputeArgumentsCreate);
  if (ier) return ier;

  // nothing else to do for this case

  // everything is good
  ier = false;
  return ier;
}

//******************************************************************************
int RDPImplementation::ComputeArgumentsDestroy(
    KIM::ModelComputeArgumentsDestroy * const modelComputeArgumentsDestroy)
  const
{
  int ier;

  // nothing else to do for this case

  // everything is good
  ier = false;
  return ier;
}




//==============================================================================
//
// Implementation of RDPImplementation private member functions
//
//==============================================================================

//******************************************************************************
#include "KIM_ModelDriverCreateLogMacros.hpp"
int RDPImplementation::OpenParameterFiles(
    KIM::ModelDriverCreate * const modelDriverCreate,
    int const numberParameterFiles,
    FILE* parameterFilePointers[MAX_PARAMETER_FILES])
{
  int ier;

  if (numberParameterFiles > MAX_PARAMETER_FILES) {
    ier = true;
    LOG_ERROR("RDP given too many parameter files");
    return ier;
  }

  for (int i = 0; i < numberParameterFiles; ++i) {
    std::string const * paramFileName;
    ier = modelDriverCreate->GetParameterFileName(i, &paramFileName);
    if (ier) {
      LOG_ERROR("Unable to get parameter file name");
      return ier;
    }

    parameterFilePointers[i] = fopen(paramFileName->c_str(), "r");
    if (parameterFilePointers[i] == 0) {
      char message[MAXLINE];
      sprintf(message,
              "RDP parameter file number %d cannot be opened",
              i);
      ier = true;
      LOG_ERROR(message);
      for (int j = i - 1; i <= 0; --i) {
        fclose(parameterFilePointers[j]);
      }
      return ier;
    }
  }

  // everything is good
  ier = false;
  return ier;
}

//******************************************************************************
int RDPImplementation::ProcessParameterFiles(
    KIM::ModelDriverCreate * const modelDriverCreate,
    int const numberParameterFiles,
    FILE* const parameterFilePointers[MAX_PARAMETER_FILES])
{
  int N, ier;
  int endOfFileFlag = 0;
  char spec1[MAXLINE], spec2[MAXLINE], nextLine[MAXLINE];
  int iIndex, jIndex , indx;
  double next_A, next_B, next_p, next_q, next_sigma, next_lambda, next_gamma;
  double next_costheta0, next_cutoff;

  getNextDataLine(parameterFilePointers[0], nextLine, MAXLINE, &endOfFileFlag);
  ier = sscanf(nextLine, "%d", &N);
  if (ier != 1) {
    sprintf(nextLine, "unable to read first line of the parameter file");
    ier = true;
    LOG_ERROR(nextLine);
    fclose(parameterFilePointers[0]);
    return ier;
  }
  numberModelSpecies_ = N;
  numberUniqueSpeciesPairs_ = ((numberModelSpecies_+1)*numberModelSpecies_)/2;

  // allocate memory for all parameters
  AllocateFreeParameterMemory();

  // set all values of p_ to -1.1e10 for later check that we have read all params
  for (int i = 0; i < ((N+1)*N/2); i++) {
    p_[i]  = -1.1e10;
  }

  // keep track of known species
  std::map<KIM::SpeciesName const, int, KIM::SPECIES_NAME::Comparator>
      modelSpeciesMap;
  int index = 0;   // species code integer code starting from 0

  // Read and process data lines
  getNextDataLine(parameterFilePointers[0], nextLine, MAXLINE, &endOfFileFlag);
  while (endOfFileFlag == 0)
  {
    ier = sscanf(nextLine, "%s %s %lg %lg %lg %lg %lg %lg %lg %lg %lg",
                 spec1, spec2, &next_A, &next_B, &next_p, &next_q, &next_sigma,
                 &next_lambda, &next_gamma, &next_costheta0, &next_cutoff);
    if (ier != 11) {
      sprintf(nextLine, "error reading lines of the parameter file");
      LOG_ERROR(nextLine);
      return true;
    }

    // convert species strings to proper type instances
    KIM::SpeciesName const specName1(spec1);
    KIM::SpeciesName const specName2(spec2);

    // check for new species
    std::map<KIM::SpeciesName const, int, KIM::SPECIES_NAME::Comparator>::
        const_iterator iIter = modelSpeciesMap.find(specName1);
    if (iIter == modelSpeciesMap.end()) {
      modelSpeciesMap[specName1] = index;
      modelSpeciesCodeList_.push_back(index);

      ier = modelDriverCreate->SetSpeciesCode(specName1, index);
      if (ier) return ier;
      iIndex = index;
      index++;
    }
    else {
      iIndex = modelSpeciesMap[specName1];
    }

    std::map<KIM::SpeciesName const, int, KIM::SPECIES_NAME::Comparator>::
        const_iterator jIter = modelSpeciesMap.find(specName2);
    if (jIter == modelSpeciesMap.end()) {
      modelSpeciesMap[specName2] = index;
      modelSpeciesCodeList_.push_back(index);

      ier = modelDriverCreate->SetSpeciesCode(specName2, index);
      if (ier) return ier;
      jIndex = index;
      index++;
    }
    else {
      jIndex = modelSpeciesMap[specName2];
    }

    if (iIndex >= jIndex) {
      indx = jIndex*N + iIndex - (jIndex*jIndex + jIndex)/2;
    }
    else {
      indx = iIndex*N + jIndex - (iIndex*iIndex + iIndex)/2;
    }
    A_[indx] = next_A;
    B_[indx] = next_B;
    p_[indx] = next_p;
    q_[indx] = next_q;
    sigma_[indx] = next_sigma;
    lambda_[indx] = next_lambda;
    gamma_[indx] = next_gamma;
    costheta0_[indx] = next_costheta0;
    cutoff_[indx] = next_cutoff;

    getNextDataLine(parameterFilePointers[0], nextLine, MAXLINE, &endOfFileFlag);
  }

  // check we have read all parameters
  for (int i = 0; i < ((N+1)*N/2); i++) {
    if (p_[i] < -1e10) {
      sprintf(nextLine, "error: not enough parameter data.\n");
      sprintf(nextLine, "%d species requires %d data lines.", N, (N+1)*N/2);
      LOG_ERROR(nextLine);
      return true;
    }
  }

  // everything is good
  ier = false;
  return ier;
}

//******************************************************************************
void RDPImplementation::getNextDataLine(
    FILE* const filePtr, char* nextLinePtr, int const maxSize,
    int *endOfFileFlag)
{
  do
  {
    if(fgets(nextLinePtr, maxSize, filePtr) == NULL) {
      *endOfFileFlag = 1;
      break;
    }

    while ((nextLinePtr[0] == ' ' || nextLinePtr[0] == '\t') ||
           (nextLinePtr[0] == '\n' || nextLinePtr[0] == '\r' ))
    {
      nextLinePtr = (nextLinePtr + 1);
    }
  }

  while ((strncmp("#", nextLinePtr, 1) == 0) || (strlen(nextLinePtr) == 0));
}

//******************************************************************************
void RDPImplementation::CloseParameterFiles(
    int const numberParameterFiles,
    FILE* const parameterFilePointers[MAX_PARAMETER_FILES])
{
  for (int i = 0; i < numberParameterFiles; ++i)
    fclose(parameterFilePointers[i]);
}

//******************************************************************************
void RDPImplementation::AllocateFreeParameterMemory()
{ // allocate memory for data
  AllocateAndInitialize1DArray(cutoff_, numberUniqueSpeciesPairs_);
  // TODO add parameters
  AllocateAndInitialize1DArray(A_, numberUniqueSpeciesPairs_);

  AllocateAndInitialize2DArray(cutoffSq_2D_, numberModelSpecies_, numberModelSpecies_);
  // TODO add parameters
  AllocateAndInitialize2DArray(A_2D_, numberModelSpecies_, numberModelSpecies_);

}

//******************************************************************************
#include "KIM_ModelDriverCreateLogMacros.hpp"
int RDPImplementation::ConvertUnits(
    KIM::ModelDriverCreate * const modelDriverCreate,
    KIM::LengthUnit const requestedLengthUnit,
    KIM::EnergyUnit const requestedEnergyUnit,
    KIM::ChargeUnit const requestedChargeUnit,
    KIM::TemperatureUnit const requestedTemperatureUnit,
    KIM::TimeUnit const requestedTimeUnit)
{
  int ier;

  // define default base units
  KIM::LengthUnit fromLength = KIM::LENGTH_UNIT::A;
  KIM::EnergyUnit fromEnergy = KIM::ENERGY_UNIT::eV;
  KIM::ChargeUnit fromCharge = KIM::CHARGE_UNIT::e;
  KIM::TemperatureUnit fromTemperature = KIM::TEMPERATURE_UNIT::K;
  KIM::TimeUnit fromTime = KIM::TIME_UNIT::ps;

  // changing length units
  double convertLength = 1.0;
  ier = modelDriverCreate->ConvertUnit(
      fromLength, fromEnergy, fromCharge, fromTemperature, fromTime,
      requestedLengthUnit, requestedEnergyUnit, requestedChargeUnit,
      requestedTemperatureUnit, requestedTimeUnit,
      1.0, 0.0, 0.0, 0.0, 0.0,
      &convertLength);
  if (ier) {
    LOG_ERROR("Unable to convert length unit");
    return ier;
  }
  // convert to active units
  if (convertLength != ONE) {
    for (int i = 0; i < numberUniqueSpeciesPairs_; ++i) {
      cutoff_[i] *= convertLength;
      // TODO add more parameters
    }
  }

  // changing energy units
  double convertEnergy = 1.0;
  ier = modelDriverCreate->ConvertUnit(
      fromLength, fromEnergy, fromCharge, fromTemperature, fromTime,
      requestedLengthUnit, requestedEnergyUnit, requestedChargeUnit,
      requestedTemperatureUnit, requestedTimeUnit,
      0.0, 1.0, 0.0, 0.0, 0.0,
      &convertEnergy);
  if (ier) {
    LOG_ERROR("Unable to convert energy unit");
    return ier;
  }
  // convert to active units
  if (convertLength != ONE) {
    for (int i = 0; i < numberUniqueSpeciesPairs_; ++i) {
      A_[i] *= convertEnergy;
      // TODO add more parameters
    }
  }

  // register units
  ier = modelDriverCreate->SetUnits(
      requestedLengthUnit,
      requestedEnergyUnit,
      KIM::CHARGE_UNIT::unused,
      KIM::TEMPERATURE_UNIT::unused,
      KIM::TIME_UNIT::unused);
  if (ier) {
    LOG_ERROR("Unable to set units to requested values");
    return ier;
  }

  // everything is good
  ier = false;
  return ier;
}

//******************************************************************************
int RDPImplementation::RegisterKIMModelSettings(
    KIM::ModelDriverCreate * const modelDriverCreate) const
{
  // register numbering
  int error = modelDriverCreate->SetModelNumbering(KIM::NUMBERING::zeroBased);

  return error;
}

//******************************************************************************
#include "KIM_ModelComputeArgumentsCreateLogMacros.hpp"
int RDPImplementation::RegisterKIMComputeArgumentsSettings(
    KIM::ModelComputeArgumentsCreate * const modelComputeArgumentsCreate) const
{
  // TODO modifiy accordingly
  // register arguments
  LOG_INFORMATION("Register argument supportStatus");
  int error =
      modelComputeArgumentsCreate->SetArgumentSupportStatus(
          KIM::COMPUTE_ARGUMENT_NAME::partialEnergy,
          KIM::SUPPORT_STATUS::optional)
      || modelComputeArgumentsCreate->SetArgumentSupportStatus(
          KIM::COMPUTE_ARGUMENT_NAME::partialForces,
          KIM::SUPPORT_STATUS::optional)
      || modelComputeArgumentsCreate->SetArgumentSupportStatus(
          KIM::COMPUTE_ARGUMENT_NAME::partialParticleEnergy,
          KIM::SUPPORT_STATUS::optional);

  // TODO modifiy accordingly
  // register callbacks
  LOG_INFORMATION("Register callback supportStatus");
  error = error
      || modelComputeArgumentsCreate->SetCallbackSupportStatus(
          KIM::COMPUTE_CALLBACK_NAME::ProcessDEDrTerm,
          KIM::SUPPORT_STATUS::optional)
      || modelComputeArgumentsCreate->SetCallbackSupportStatus(
          KIM::COMPUTE_CALLBACK_NAME::ProcessD2EDr2Term,
          KIM::SUPPORT_STATUS::optional);

  return error;
}

//******************************************************************************
#include "KIM_ModelDriverCreateLogMacros.hpp"
int RDPImplementation::RegisterKIMParameters(
    KIM::ModelDriverCreate * const modelDriverCreate)
{
  int ier = false;

  // publish parameters (order is important)
  // TODO add parameters
  ier = modelDriverCreate->SetParameterPointer(numberUniqueSpeciesPairs_, A_, "A")
     || modelDriverCreate->SetParameterPointer(numberUniqueSpeciesPairs_, cutoff_, "cutoff");
  if (ier) {
    LOG_ERROR("set_parameters");
    return ier;
  }

  // everything is good
  ier = false;
  return ier;
}

//******************************************************************************
int RDPImplementation::RegisterKIMFunctions(
    KIM::ModelDriverCreate * const modelDriverCreate)
    const
{
  int error;

  // register the Destroy(), Refresh(), and Compute() functions
  error = modelDriverCreate->SetDestroyPointer(
      KIM::LANGUAGE_NAME::cpp, (KIM::func*) &(RDP::Destroy))
      || modelDriverCreate->SetRefreshPointer(
          KIM::LANGUAGE_NAME::cpp, (KIM::func*) &(RDP::Refresh))
      || modelDriverCreate->SetComputePointer(
          KIM::LANGUAGE_NAME::cpp, (KIM::func*) &(RDP::Compute))
      || modelDriverCreate->SetComputeArgumentsCreatePointer(
          KIM::LANGUAGE_NAME::cpp,
          (KIM::func*) &(RDP::ComputeArgumentsCreate))
      || modelDriverCreate->SetComputeArgumentsDestroyPointer(
          KIM::LANGUAGE_NAME::cpp,
          (KIM::func*) &(RDP::ComputeArgumentsDestroy));

  return error;
}

//******************************************************************************
template<class ModelObj>
int RDPImplementation::SetRefreshMutableValues(
    ModelObj * const modelObj)
{ // use (possibly) new values of free parameters to compute other quantities
  int ier;

  // update parameters
  for (int i = 0; i < numberModelSpecies_; ++i) {
    for (int j = 0; j <= i ; ++j) {
      int const index = j*numberModelSpecies_ + i - (j*j + j)/2;
      cutoffSq_2D_[i][j] = cutoffSq_2D_[j][i] = cutoff_[index] * cutoff_[index];
      A_2D_[i][j] = A_2D_[j][i] = A_[index];
      //TODO add parameters
    }
  }

  // update cutoff value in KIM API object
  influenceDistance_ = 0.0;

  for (int i = 0; i < numberModelSpecies_; i++) {
    int indexI = modelSpeciesCodeList_[i];

    for (int j = 0; j < numberModelSpecies_; j++) {
      int indexJ = modelSpeciesCodeList_[j];

      if (influenceDistance_ < cutoffSq_2D_[indexI][indexJ]) {
        influenceDistance_ = cutoffSq_2D_[indexI][indexJ];
      }
    }
  }

  influenceDistance_ = sqrt(influenceDistance_);
  modelObj->SetInfluenceDistancePointer(&influenceDistance_);
  modelObj->SetNeighborListCutoffsPointer(1, &influenceDistance_);


  // everything is good
  ier = false;
  return ier;
}

//******************************************************************************
#include "KIM_ModelComputeArgumentsLogMacros.hpp"
int RDPImplementation::SetComputeMutableValues(
    KIM::ModelComputeArguments const * const modelComputeArguments,
    bool& isComputeProcess_dEdr,
    bool& isComputeProcess_d2Edr2,
    bool& isComputeEnergy,
    bool& isComputeForces,
    bool& isComputeParticleEnergy,
    int const*& particleSpeciesCodes,
    int const*& particleContributing,
    VectorOfSizeDIM const*& coordinates,
    double*& energy,
    double*& particleEnergy,
    VectorOfSizeDIM*& forces)
{
  int ier = true;

  // get compute flags
  int compProcess_dEdr;
  int compProcess_d2Edr2;

  modelComputeArguments->IsCallbackPresent(
      KIM::COMPUTE_CALLBACK_NAME::ProcessDEDrTerm,
      &compProcess_dEdr);
  modelComputeArguments->IsCallbackPresent(
      KIM::COMPUTE_CALLBACK_NAME::ProcessD2EDr2Term,
      &compProcess_d2Edr2);

  isComputeProcess_dEdr = compProcess_dEdr;
  isComputeProcess_d2Edr2 = compProcess_d2Edr2;


  int const* numberOfParticles;
  ier =
      modelComputeArguments->GetArgumentPointer(
          KIM::COMPUTE_ARGUMENT_NAME::numberOfParticles,
          &numberOfParticles)
      || modelComputeArguments->GetArgumentPointer(
          KIM::COMPUTE_ARGUMENT_NAME::particleSpeciesCodes,
          &particleSpeciesCodes)
      || modelComputeArguments->GetArgumentPointer(
          KIM::COMPUTE_ARGUMENT_NAME::particleContributing,
          &particleContributing)
      || modelComputeArguments->GetArgumentPointer(
          KIM::COMPUTE_ARGUMENT_NAME::coordinates,
          (double const ** const) &coordinates)
      || modelComputeArguments->GetArgumentPointer(
          KIM::COMPUTE_ARGUMENT_NAME::partialEnergy,
          &energy)
      || modelComputeArguments->GetArgumentPointer(
          KIM::COMPUTE_ARGUMENT_NAME::partialParticleEnergy,
          &particleEnergy)
      || modelComputeArguments->GetArgumentPointer(
          KIM::COMPUTE_ARGUMENT_NAME::partialForces,
          (double const ** const) &forces);
  if (ier)
  {
    LOG_ERROR("GetArgumentPointer");
    return ier;
  }

  isComputeEnergy = (energy != 0);
  isComputeParticleEnergy = (particleEnergy != 0);
  isComputeForces = (forces != 0);

  // update values
  cachedNumberOfParticles_ = *numberOfParticles;

  // everything is good
  ier = false;
  return ier;
}

//******************************************************************************
// Assume that the particle species interge code starts from 0
#include "KIM_ModelComputeLogMacros.hpp"
int RDPImplementation::CheckParticleSpecies(
    KIM::ModelCompute const * const modelCompute,
    int const* const particleSpecies)
    const
{
  int ier;
  for (int i = 0; i < cachedNumberOfParticles_; ++i) {
    if ((particleSpecies[i] < 0) || (particleSpecies[i] >= numberModelSpecies_)) {
      ier = true;
      LOG_ERROR("unsupported particle species detected");
      return ier;
    }
  }

  // everything is good
  ier = false;
  return ier;
}

//******************************************************************************
int RDPImplementation::GetComputeIndex(
    const bool& isComputeProcess_dEdr,
    const bool& isComputeProcess_d2Edr2,
    const bool& isComputeEnergy,
    const bool& isComputeForces,
    const bool& isComputeParticleEnergy) const
{
  const int processd2E = 2;
  const int energy = 2;
  const int force = 2;
  const int particleEnergy = 2;


  int index = 0;

  // processdE
  index += (int(isComputeProcess_dEdr)) * processd2E * energy * force * particleEnergy;

  // processd2E
  index += (int(isComputeProcess_d2Edr2)) * energy * force * particleEnergy;

  // energy
  index += (int(isComputeEnergy)) * force * particleEnergy;

  // force
  index += (int(isComputeForces)) * particleEnergy;

  // particleEnergy
  index += (int(isComputeParticleEnergy));


  return index;
}


//==============================================================================
//
// RDP functions
//
//==============================================================================
//TODO add implementation of your functions


