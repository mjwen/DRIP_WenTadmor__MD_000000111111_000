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
// Copyright (c) 2013--2018, Regents of the University of Minnesota.
// All rights reserved.
//
// Contributors:
//    Mingjian Wen
//


#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>

#include "RDPImplementation.hpp"
#include "KIM_ModelDriverHeaders.hpp"

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
  cutoff_(NULL),
  rhocut_(NULL),
  C0_(NULL),
  C2_(NULL),
  C4_(NULL),
  C_(NULL),
  delta_(NULL),
  lambda_(NULL),
  A_(NULL),
  z0_(NULL),
  B_(NULL),
  eta_(NULL),
  influenceDistance_(0.0),
  paddingNeighborHints_(0),
  halfListHints_(0),
  cutoffSq_2D_(NULL),
  rhocutSq_2D_(NULL),
  C0_2D_(NULL),
  C2_2D_(NULL),
  C4_2D_(NULL),
  C_2D_(NULL),
  delta_2D_(NULL),
  lambda_2D_(NULL),
  A_2D_(NULL),
  z0_2D_(NULL),
  B_2D_(NULL),
  eta_2D_(NULL),
  cachedNumberOfParticles_(0),
  cachedIsComputeEnergy_(0),
  cachedIsComputeParticleEnergy_(0),
  cachedIsComputeForces_(0),
  cachedIsComputeProcess_dEdr_(0),
  cachedIsComputeProcess_d2Edr2_(0)
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

  Deallocate1DArray<double> (cutoff_);
  Deallocate1DArray<double> (rhocut_);
  Deallocate1DArray<double> (C0_);
  Deallocate1DArray<double> (C2_);
  Deallocate1DArray<double> (C4_);
  Deallocate1DArray<double> (C_);
  Deallocate1DArray<double> (delta_);
  Deallocate1DArray<double> (lambda_);
  Deallocate1DArray<double> (A_);
  Deallocate1DArray<double> (z0_);
  Deallocate1DArray<double> (B_);
  Deallocate1DArray<double> (eta_);

  Deallocate2DArray<double> (cutoffSq_2D_);
  Deallocate2DArray<double> (rhocutSq_2D_);
  Deallocate2DArray<double> (C0_2D_);
  Deallocate2DArray<double> (C2_2D_);
  Deallocate2DArray<double> (C4_2D_);
  Deallocate2DArray<double> (C_2D_);
  Deallocate2DArray<double> (delta_2D_);
  Deallocate2DArray<double> (lambda_2D_);
  Deallocate2DArray<double> (A_2D_);
  Deallocate2DArray<double> (z0_2D_);
  Deallocate2DArray<double> (B_2D_);
  Deallocate2DArray<double> (eta_2D_);
}

//******************************************************************************
int RDPImplementation::Refresh( KIM::ModelRefresh * const modelRefresh)
{
  int ier;

  ier = SetRefreshMutableValues(modelRefresh);
  if (ier) return ier;

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
  bool isComputeVirial = false;
  bool isComputeParticleVirial = false;
  //
  // KIM API Model Input
  int const* particleSpeciesCodes = NULL;
  int const* particleContributing = NULL;
  VectorOfSizeDIM const* coordinates = NULL;
  //
  // KIM API Model Output
  double* energy = NULL;
  double* particleEnergy = NULL;
  VectorOfSizeDIM* forces = NULL;
  VectorOfSizeSix* virial = NULL;
  VectorOfSizeSix* particleVirial = NULL;
  ier = SetComputeMutableValues(modelComputeArguments,
      isComputeProcess_dEdr, isComputeProcess_d2Edr2,
      isComputeEnergy, isComputeForces, isComputeParticleEnergy,
      isComputeVirial, isComputeParticleVirial,
      particleSpeciesCodes, particleContributing, coordinates,
      energy, forces, particleEnergy, virial, particleVirial);
  if (ier) return ier;

  // Skip this check for efficiency
  //
  //ier = CheckParticleSpecies(modelComputeArguments, particleSpeciesCodes);
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

///******************************************************************************
void RDPImplementation::AllocatePrivateParameterMemory()
{
  // nothing to do for this case
}


//******************************************************************************
void RDPImplementation::AllocateParameterMemory()
{ // allocate memory for data
  AllocateAndInitialize1DArray<double> (cutoff_, numberUniqueSpeciesPairs_);
  AllocateAndInitialize1DArray<double> (rhocut_, numberUniqueSpeciesPairs_);
  AllocateAndInitialize1DArray<double> (C0_, numberUniqueSpeciesPairs_);
  AllocateAndInitialize1DArray<double> (C2_, numberUniqueSpeciesPairs_);
  AllocateAndInitialize1DArray<double> (C4_, numberUniqueSpeciesPairs_);
  AllocateAndInitialize1DArray<double> (C_, numberUniqueSpeciesPairs_);
  AllocateAndInitialize1DArray<double> (delta_, numberUniqueSpeciesPairs_);
  AllocateAndInitialize1DArray<double> (lambda_, numberUniqueSpeciesPairs_);
  AllocateAndInitialize1DArray<double> (A_, numberUniqueSpeciesPairs_);
  AllocateAndInitialize1DArray<double> (z0_, numberUniqueSpeciesPairs_);
  AllocateAndInitialize1DArray<double> (B_, numberUniqueSpeciesPairs_);
  AllocateAndInitialize1DArray<double> (eta_, numberUniqueSpeciesPairs_);

  AllocateAndInitialize2DArray<double> (cutoffSq_2D_, numberModelSpecies_, numberModelSpecies_);
  AllocateAndInitialize2DArray<double> (rhocutSq_2D_, numberModelSpecies_, numberModelSpecies_);
  AllocateAndInitialize2DArray<double> (C0_2D_, numberModelSpecies_, numberModelSpecies_);
  AllocateAndInitialize2DArray<double> (C2_2D_, numberModelSpecies_, numberModelSpecies_);
  AllocateAndInitialize2DArray<double> (C4_2D_, numberModelSpecies_, numberModelSpecies_);
  AllocateAndInitialize2DArray<double> (C_2D_, numberModelSpecies_, numberModelSpecies_);
  AllocateAndInitialize2DArray<double> (delta_2D_, numberModelSpecies_, numberModelSpecies_);
  AllocateAndInitialize2DArray<double> (lambda_2D_, numberModelSpecies_, numberModelSpecies_);
  AllocateAndInitialize2DArray<double> (A_2D_, numberModelSpecies_, numberModelSpecies_);
  AllocateAndInitialize2DArray<double> (z0_2D_, numberModelSpecies_, numberModelSpecies_);
  AllocateAndInitialize2DArray<double> (B_2D_, numberModelSpecies_, numberModelSpecies_);
  AllocateAndInitialize2DArray<double> (eta_2D_, numberModelSpecies_, numberModelSpecies_);

}

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
#include "KIM_ModelDriverCreateLogMacros.hpp"
int RDPImplementation::ProcessParameterFiles(
    KIM::ModelDriverCreate * const modelDriverCreate,
    int const numberParameterFiles,
    FILE* const parameterFilePointers[MAX_PARAMETER_FILES])
{
  int N, ier;
  int numberOfLinesRead;
  int endOfFileFlag = 0;
  char spec1[MAXLINE], nextLine[MAXLINE];
  int iIndex, jIndex , indx;
  double next_C0, next_C2, next_C4, next_C, next_delta, next_lambda, next_A;
  double next_z0, next_B, next_eta, next_rhocut, next_cutoff;


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
  AllocateParameterMemory();


  // keep track of known species
  std::map<KIM::SpeciesName const, int, KIM::SPECIES_NAME::Comparator> modelSpeciesMap;
  int index = 0;   // species code integer code starting from 0

  // Read and process data lines
  numberOfLinesRead = 0;
  getNextDataLine(parameterFilePointers[0], nextLine, MAXLINE, &endOfFileFlag);
  while (endOfFileFlag == 0)
  {
    ier = sscanf(nextLine, "%s %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg",
                 spec1, &next_C0, &next_C2, &next_C4, &next_C, &next_delta,
                 &next_lambda, &next_A, &next_z0, &next_B, &next_eta,
                 &next_rhocut, &next_cutoff);
    if (ier != 13) {
      sprintf(nextLine, "error reading lines of the parameter file");
      LOG_ERROR(nextLine);
      return true;
    }

    // convert species strings to proper type instances
    KIM::SpeciesName const specName1(spec1);
    KIM::SpeciesName const specName2(spec1);

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
    C0_[indx] = next_C0;
    C2_[indx] = next_C2;
    C4_[indx] = next_C4;
    C_[indx] = next_C;
    delta_[indx] = next_delta;
    lambda_[indx] = next_lambda;
    A_[indx] = next_A;
    z0_[indx] = next_z0;
    B_[indx] = next_B;
    eta_[indx] = next_eta;
    rhocut_[indx] = next_rhocut;
    cutoff_[indx] = next_cutoff;

    numberOfLinesRead += 1;
    getNextDataLine(parameterFilePointers[0], nextLine, MAXLINE, &endOfFileFlag);
  }

  // check we have read all parameters
  if ( (N+1)*N/2 != numberOfLinesRead ) {
    sprintf(nextLine, "error in parameter file.\n");
    LOG_ERROR(nextLine);
    return true;

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
      rhocut_[i] *= convertLength;
      delta_[i] *= convertLength;
      z0_[i] *= convertLength;
      lambda_[i] /= convertLength;
      eta_[i] /= convertLength;
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
      C0_[i] *= convertEnergy;
      C2_[i] *= convertEnergy;
      C4_[i] *= convertEnergy;
      C_[i] *= convertEnergy;
      B_[i] *= convertEnergy;
      A_[i] *= convertEnergy;
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
  // register arguments
  LOG_INFORMATION("Register argument supportStatus");


  int error =
    modelComputeArgumentsCreate->SetArgumentSupportStatus(
        KIM::COMPUTE_ARGUMENT_NAME::partialEnergy,
        KIM::SUPPORT_STATUS::optional) ||
    modelComputeArgumentsCreate->SetArgumentSupportStatus(
        KIM::COMPUTE_ARGUMENT_NAME::partialForces,
        KIM::SUPPORT_STATUS::optional) ||
    modelComputeArgumentsCreate->SetArgumentSupportStatus(
        KIM::COMPUTE_ARGUMENT_NAME::partialParticleEnergy,
        KIM::SUPPORT_STATUS::optional) ||
    modelComputeArgumentsCreate->SetArgumentSupportStatus(
        KIM::COMPUTE_ARGUMENT_NAME::partialVirial,
        KIM::SUPPORT_STATUS::optional) ||
    modelComputeArgumentsCreate->SetArgumentSupportStatus(
        KIM::COMPUTE_ARGUMENT_NAME::partialParticleVirial,
        KIM::SUPPORT_STATUS::optional);

  // register callbacks
  LOG_INFORMATION("Register callback supportStatus");
  error = error
      || modelComputeArgumentsCreate->SetCallbackSupportStatus(
          KIM::COMPUTE_CALLBACK_NAME::ProcessDEDrTerm,
          KIM::SUPPORT_STATUS::notSupported)
      || modelComputeArgumentsCreate->SetCallbackSupportStatus(
          KIM::COMPUTE_CALLBACK_NAME::ProcessD2EDr2Term,
          KIM::SUPPORT_STATUS::notSupported);

  return error;
}

//******************************************************************************
#include "KIM_ModelDriverCreateLogMacros.hpp"
int RDPImplementation::RegisterKIMParameters(
    KIM::ModelDriverCreate * const modelDriverCreate)
{
  int ier = false;

  // publish parameters (order is important)
  ier = modelDriverCreate->SetParameterPointer(
  numberUniqueSpeciesPairs_, C0_, "C0")
     || modelDriverCreate->SetParameterPointer(
     numberUniqueSpeciesPairs_, C2_, "C2")
     || modelDriverCreate->SetParameterPointer(
     numberUniqueSpeciesPairs_, C4_, "C4")
     || modelDriverCreate->SetParameterPointer(
     numberUniqueSpeciesPairs_, C_, "C")
     || modelDriverCreate->SetParameterPointer(
     numberUniqueSpeciesPairs_, delta_, "delta")
     || modelDriverCreate->SetParameterPointer(
     numberUniqueSpeciesPairs_, lambda_, "lambda")
     || modelDriverCreate->SetParameterPointer(
     numberUniqueSpeciesPairs_, A_, "A")
     || modelDriverCreate->SetParameterPointer(
     numberUniqueSpeciesPairs_, z0_, "z0")
     || modelDriverCreate->SetParameterPointer(
     numberUniqueSpeciesPairs_, B_, "B")
     || modelDriverCreate->SetParameterPointer(
     numberUniqueSpeciesPairs_, eta_, "eta")
     || modelDriverCreate->SetParameterPointer(
     numberUniqueSpeciesPairs_, rhocut_, "rhocut")
     || modelDriverCreate->SetParameterPointer(
     numberUniqueSpeciesPairs_, cutoff_, "cutoff");
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
  error =
      modelDriverCreate->SetDestroyPointer(
      KIM::LANGUAGE_NAME::cpp,
      (KIM::func*) &(RDP::Destroy))
      || modelDriverCreate->SetRefreshPointer(
          KIM::LANGUAGE_NAME::cpp,
          (KIM::func*) &(RDP::Refresh))
      || modelDriverCreate->SetComputePointer(
          KIM::LANGUAGE_NAME::cpp,
          (KIM::func*) &(RDP::Compute))
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
{ // use (possibly) new values of parameters to compute other quantities
  // NOTE: This function is templated because it's called with both a
  //       modelDriverCreate object during initialization and with a
  //       modelRefresh object when the Model's parameters have been altered
  int ier;

  // update parameters
  for (int i = 0; i < numberModelSpecies_; ++i) {
    for (int j = 0; j <= i ; ++j) {
      int const index = j*numberModelSpecies_ + i - (j*j + j)/2;
      cutoffSq_2D_[i][j] = cutoffSq_2D_[j][i] = cutoff_[index] * cutoff_[index];
      rhocutSq_2D_[i][j] = rhocutSq_2D_[j][i] = rhocut_[index] * rhocut_[index];
      C0_2D_[i][j] = C0_2D_[j][i] = C0_[index];
      C2_2D_[i][j] = C2_2D_[j][i] = C2_[index];
      C4_2D_[i][j] = C4_2D_[j][i] = C4_[index];
      C_2D_[i][j] = C_2D_[j][i] = C_[index];
      delta_2D_[i][j] = delta_2D_[j][i] = delta_[index];
      lambda_2D_[i][j] = lambda_2D_[j][i] = lambda_[index];
      A_2D_[i][j] = A_2D_[j][i] = A_[index];
      z0_2D_[i][j] = z0_2D_[j][i] = z0_[index];
      B_2D_[i][j] = B_2D_[j][i] = B_[index];
      eta_2D_[i][j] = eta_2D_[j][i] = eta_[index];
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
  modelObj->SetNeighborListPointers(1,
      &influenceDistance_, &paddingNeighborHints_, &halfListHints_);


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
    bool& isComputeVirial,
    bool& isComputeParticleVirial,
    int const*& particleSpeciesCodes,
    int const*& particleContributing,
    VectorOfSizeDIM const*& coordinates,
    double*& energy,
    VectorOfSizeDIM*& forces,
    double*& particleEnergy,
    VectorOfSizeSix*& virial,
    VectorOfSizeSix*& particleVirial)
{
  int ier = true;


  // callback not supported
  isComputeProcess_dEdr = false;
  isComputeProcess_d2Edr2 = false;


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
      ||  modelComputeArguments->GetArgumentPointer(
        KIM::COMPUTE_ARGUMENT_NAME::partialForces,
        (double const** const)&forces) ||
    modelComputeArguments->GetArgumentPointer(
        KIM::COMPUTE_ARGUMENT_NAME::partialParticleEnergy,
        &particleEnergy) ||
    modelComputeArguments->GetArgumentPointer(
        KIM::COMPUTE_ARGUMENT_NAME::partialVirial,
        (double const** const)&virial) ||
    modelComputeArguments->GetArgumentPointer(
        KIM::COMPUTE_ARGUMENT_NAME::partialParticleVirial,
        (double const** const)&particleVirial);
  if (ier) {
    LOG_ERROR("GetArgumentPointer");
    return ier;
  }

  isComputeEnergy = (energy != NULL);
  isComputeForces = (forces != NULL);
  isComputeParticleEnergy = (particleEnergy != NULL);
  isComputeVirial = (virial != NULL);
  isComputeParticleVirial = (particleVirial != NULL);

  // update values
  cachedNumberOfParticles_ = *numberOfParticles;

  // everything is good
  ier = false;
  return ier;
}

//******************************************************************************
// Assume that the particle species interge code starts from 0
#include "KIM_ModelComputeLogMacros.hpp"
int RDPImplementation::CheckParticleSpeciesCodes(
    KIM::ModelCompute const * const modelCompute,
    int const* const particleSpeciesCodes) const
{
  int ier;
  for (int i = 0; i < cachedNumberOfParticles_; ++i) {
    if ((particleSpeciesCodes[i] < 0) || (particleSpeciesCodes[i] >= numberModelSpecies_)) {
      ier = true;
      LOG_ERROR("unsupported particle species codes detected");
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
    const bool& isComputeParticleEnergy,
    const bool& isComputeVirial,
    const bool& isComputeParticleVirial) const
{
  //const int processdE = 2;
  const int processd2E = 2;
  const int energy = 2;
  const int force = 2;
  const int particleEnergy = 2;
  const int virial = 2;
  const int particleVirial = 2;


  int index = 0;

  // processdE
  index += (int(isComputeProcess_dEdr))
           * processd2E * energy * force * particleEnergy * virial * particleVirial;

  // processd2E
  index += (int(isComputeProcess_d2Edr2))
           * energy * force * particleEnergy * virial * particleVirial;

  // energy
  index += (int(isComputeEnergy))
           * force * particleEnergy * virial * particleVirial;

  // force
  index += (int(isComputeForces))
           * particleEnergy * virial * particleVirial;

  // particleEnergy
  index += (int(isComputeParticleEnergy))
           * virial * particleVirial;

  // virial
  index += (int(isComputeVirial))
           * particleVirial;

  // particleVirial
  index += (int(isComputeParticleVirial));


  return index;
}


//==============================================================================
//
// RDP functions
//
//==============================================================================

// create layers and find the nearest 3 neighbors of each atom
int RDPImplementation::create_layers(KIM::ModelCompute const * const modelCompute,
    KIM::ModelComputeArguments const * const modelComputeArguments,
    VectorOfSizeDIM const * const coordinates,
    int * const in_layer,
    VectorOfSizeThreeInt * const nearest3neigh) const
{


  int const Natoms = cachedNumberOfParticles_;
  int ier = true;

  // no atoms, which can happen when doing domain decompisition
  if (Natoms <= 0) return ier;

  // cutoff to separate atoms into layers (may need to be adjusted)
  double cutsq_layer = (0.72*3.44)*(0.72*3.44);

  // layer contains which atoms
  int* layer;
  AllocateAndInitialize1DArray<int> (layer, Natoms);

  // init atoms layer to -1 (not in any layer)
  for (int k=0; k<Natoms; k++) {
    in_layer[k] = -1;
  }

  int nlayers = 1;
  int nremain = Natoms; // number of atoms not included in any layer

  while(true) {  // to create all layers

    // init all atoms in current layer have -1: not classified in any layer
    // We do not always have Natoms (except the very first loop), we allocate layer
    // to the size of nAtom just because we do not want to reallocate it each time.
    for (int k=0; k<Natoms; k++) {
      layer[k] = -1;
    }

    // find an atom not incldued in any layer and start with it
    int currentLayer = nlayers - 1;
    for (int k=0; k<Natoms; k++) {
      if (in_layer[k] == -1) {
        in_layer[k] = currentLayer;
        layer[0] = k; // first atom in current layer
        break;
      }
    }

    int nin = 1; // number of atoms in current layer
    int ii = 0;
    while(true) { // to find all atoms in currentLayer

      const int i = layer[ii];

      // get neighbors of atom i
      int num_neigh = 0;
      int const * ilist = 0;
      ier = modelComputeArguments->GetNeighborList(0, i, &num_neigh, &ilist);
      if (ier) {
        LOG_ERROR("modelComputeArguments->GetNeighborList");
        return ier;
      }

      // init nb1 to be the 1st nearest neigh, nb3 the 3rd nearest
      int nb1 = -1;
      int nb2 = -1;
      int nb3 = -1;
      double nb1_rsq = 1.1e10;
      double nb2_rsq = 2.0e10;
      double nb3_rsq = 3.0e10;

      // loop over the neighbors of atom i
      for (int jj=0; jj<num_neigh; ++jj) {

        const int j = ilist[jj];

        // compute relative position vector and squared distance
        double rij[DIM];
        for (int dim = 0; dim < DIM; ++dim) {
          rij[dim] = coordinates[j][dim] - coordinates[i][dim];
        }
        // compute distance squared
        double const rsq = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];


        if (rsq < cutsq_layer) { // belongs to current layer
          if (in_layer[j] == -1) { // has not been included in some layer
            nin += 1;
            layer[nin-1] = j;
            in_layer[j] = currentLayer;
          } else {
            if(in_layer[j] != currentLayer) { // If the choice of cutsq_layer is appropriate,
                                              // this should not happen.

              LOG_ERROR("Atom included in two layers.");
              return ier;
            }
          }
        }

        // find the 3 nearest neigh
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

      } // loop on jj


      // check whether we have enough neighs to compute normal
      if (nb1_rsq >= 1.0e10 || nb2_rsq >= 1.0e10 ||  nb3_rsq >= 1.0e10) {
        ier = true;
        LOG_ERROR("no enough neigh to construct normal");
        return ier;
      }

      // store nearest 3 neigh. This is for latter use to comptue the derivatives
      // of normal; Naively, we can compute the derivatives here, but each atom has
      // 3 nerest neighs, and the derivative with respect to each neigh has 9
      // component, resulting in storing an array of size N*(1+3)*9. Seems not good.
      nearest3neigh[i][0] = nb1;
      nearest3neigh[i][1] = nb2;
      nearest3neigh[i][2] = nb3;

      // get to the next atom in current layer
      ii++;
      if (ii == nin) break;

    } // end while finding one layer

      nremain -= nin;
      if (nremain == 0) break;
      nlayers += 1;
  } // end while finding all layers

  // whether all atoms in the same layer?
  int more_than_one_layer = false;
  int currentLayer = 0;
  for (int i=0; i<Natoms; i++) {
    if (in_layer[i] != currentLayer) { // find a atom not in current layer
      more_than_one_layer = true;
      break;
    }
  }
  if (! more_than_one_layer) {
    ier = true;
    LOG_ERROR("Only one layer detected. The layer separation seems too small. "
    "You can either use larger layer separation or tune `cutsq_layer` in the Model.");
    return ier;
  }

  // dealloate
  Deallocate1DArray<int> (layer);

  return false;
}


// the r^{-6} attractive part
double RDPImplementation::calc_attractive(int const i, int const j,
    int const iSpecies, int const jSpecies, double const * const rij,
    double const r, VectorOfSizeDIM * const forces) const
{
  // compute flags
  int const comp_forces = cachedIsComputeForces_;

  // compute params
  double const z0 = z0_2D_[iSpecies][jSpecies];
  double const A = A_2D_[iSpecies][jSpecies];
  double const cutoff = sqrt(cutoffSq_2D_[iSpecies][jSpecies]);

  double roz0_sq = r*r/(z0*z0);
  double dtp;
  double tp = tap(r, cutoff, dtp);
  double r6 = A/(roz0_sq*roz0_sq*roz0_sq);
  double dr6 = -6*r6/r;
  double phi = - r6*tp;

  /* forces */
  if (comp_forces) {
    double fpair = - HALF * (r6*dtp + dr6*tp);
    for (int k=0; k<DIM; k++) {
      forces[i][k] += rij[k]*fpair/r;
      forces[j][k] -= rij[k]*fpair/r;
    }
  }

  return phi;
}


// the repulsive part
double RDPImplementation::calc_repulsive(int const i, int const j,
    int const * const particleSpeciesCodes,
    VectorOfSizeDIM const * const coordinates,
    VectorOfSizeThreeInt const * const nearest3neigh,
    const double * const rij,
    double const r, const int nbi1, const int nbi2, const int nbi3,
    double const * const ni,  VectorOfSizeDIM const * const dni_dri,
    VectorOfSizeDIM const * const dni_drnb1, VectorOfSizeDIM const * const dni_drnb2,
    VectorOfSizeDIM const * const dni_drnb3, VectorOfSizeDIM * const forces) const
{
  // compute flags
  int const comp_forces = cachedIsComputeForces_;

  // parameters
  int iSpecies = particleSpeciesCodes[i];
  int jSpecies = particleSpeciesCodes[j];
  double C0 = C0_2D_[iSpecies][jSpecies];
  double C2 = C2_2D_[iSpecies][jSpecies];
  double C4 = C4_2D_[iSpecies][jSpecies];
  double C = C_2D_[iSpecies][jSpecies];
  double delta = delta_2D_[iSpecies][jSpecies];
  double lambda = lambda_2D_[iSpecies][jSpecies];
  double z0 = z0_2D_[iSpecies][jSpecies];
  double cutoff = sqrt(cutoffSq_2D_[iSpecies][jSpecies]);

  // nearest 3 neighbors of atom j
  int nbj1 = nearest3neigh[j][0];
  int nbj2 = nearest3neigh[j][1];
  int nbj3 = nearest3neigh[j][2];


  VectorOfSizeDIM dgij_dri;
  VectorOfSizeDIM dgij_drj;
  VectorOfSizeDIM dgij_drk1;
  VectorOfSizeDIM dgij_drk2;
  VectorOfSizeDIM dgij_drk3;
  VectorOfSizeDIM dgij_drl1;
  VectorOfSizeDIM dgij_drl2;
  VectorOfSizeDIM dgij_drl3;

  VectorOfSizeDIM drhosqij_dri;
  VectorOfSizeDIM drhosqij_drj;
  VectorOfSizeDIM drhosqij_drnb1;
  VectorOfSizeDIM drhosqij_drnb2;
  VectorOfSizeDIM drhosqij_drnb3;


  // derivative of rhosq w.r.t coordinates of atoms i, j, and the nearests 3 neighs of i
  get_drhosqij(rij, ni, dni_dri, dni_drnb1, dni_drnb2, dni_drnb3, drhosqij_dri,
      drhosqij_drj, drhosqij_drnb1, drhosqij_drnb2, drhosqij_drnb3);

  // transverse decay function f(rho) and its derivative w.r.t. rhosq
  double rhosqij;
  double dtdij;
  double tdij = td(C0, C2, C4, delta, rij, r, ni, rhosqij, dtdij);

  // dihedral angle function and its derivateives
  double dgij_drhosq;
  double gij = dihedral(i, j, particleSpeciesCodes, coordinates, nearest3neigh,
      rhosqij, dgij_drhosq, dgij_dri, dgij_drj,
      dgij_drk1, dgij_drk2, dgij_drk3, dgij_drl1, dgij_drl2, dgij_drl3);

  double V2 = C + tdij + gij;

  // tap part
  double dtp;
  double tp = tap(r, cutoff, dtp);

  /* exponential part */
  double V1 = exp(-lambda*(r-z0));
  double dV1 = -V1*lambda;

  double phi = tp*V1*V2;

  // forces
  if (comp_forces) {
    for (int k=0; k<DIM; k++) {

      // forces due to derivatives of tap and V1
      double tmp =  HALF* (dtp*V1 + tp*dV1) *V2 *rij[k]/r;
      forces[i][k] += tmp;
      forces[j][k] -= tmp;

      // the following incldue the transverse decay part tdij and the dihedral part gij
      // derivative of V2 contribute to atoms i, j
      forces[i][k] -= HALF*tp*V1* ( (dtdij+dgij_drhosq)*drhosqij_dri[k] + dgij_dri[k]);
      forces[j][k] -= HALF*tp*V1* ( (dtdij+dgij_drhosq)*drhosqij_drj[k] + dgij_drj[k]);


      // derivative of V2 contribute to neighs of atom i
      forces[nbi1][k] -= HALF*tp*V1* ( (dtdij+dgij_drhosq)*drhosqij_drnb1[k]
          + dgij_drk1[k]);
      forces[nbi2][k] -= HALF*tp*V1* ( (dtdij+dgij_drhosq)*drhosqij_drnb2[k]
          + dgij_drk2[k]);
      forces[nbi3][k] -= HALF*tp*V1* ( (dtdij+dgij_drhosq)*drhosqij_drnb3[k]
          + dgij_drk3[k]);

      // derivative of V2 contribute to neighs of atom j
      forces[nbj1][k] -= HALF*tp*V1* dgij_drl1[k];
      forces[nbj2][k] -= HALF*tp*V1* dgij_drl2[k];
      forces[nbj3][k] -= HALF*tp*V1* dgij_drl3[k];

    }
  }

  return phi;
}


// compute normal and its derivatives w.r.t atom ri, and its 3 nearest neighs k1, k2, k3
void RDPImplementation::normal(const int i, VectorOfSizeDIM const * const coordinates,
    VectorOfSizeThreeInt const * const nearest3neigh,
    int & k1, int & k2, int & k3, double * const normal,
    VectorOfSizeDIM * const dn_dri, VectorOfSizeDIM * const dn_drk1,
    VectorOfSizeDIM * const dn_drk2, VectorOfSizeDIM * const dn_drk3) const
{

  k1 = nearest3neigh[i][0];
  k2 = nearest3neigh[i][1];
  k3 = nearest3neigh[i][2];

  // normal does not depend on i, setting to zero
  for (int j=0; j<DIM; j++) {
    for (int k=0; k<DIM; k++) {
      dn_dri[j][k] = 0.0;
    }
  }

  // get normal and derives of normal w.r.t to neighbors of
  deriv_cross(coordinates[k1], coordinates[k2], coordinates[k3],
      normal, dn_drk1, dn_drk2, dn_drk3);
}


void RDPImplementation::get_drhosqij(double const * const rij,
    double const * const ni,
    VectorOfSizeDIM const * const dni_dri,
    VectorOfSizeDIM const * const dni_drn1,
    VectorOfSizeDIM const * const dni_drn2,
    VectorOfSizeDIM const * const dni_drn3,
    double * const drhosq_dri, double * const drhosq_drj,
    double * const drhosq_drn1, double * const drhosq_drn2,
    double * const drhosq_drn3) const
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


// derivartive of transverse decay function f(rho) w.r.t rho
double RDPImplementation::td(double C0, double C2, double C4, double delta,
    double const * const rvec, double r, const double* const n, double & rho_sq,
    double & dtd) const
{

  double n_dot_r = dot(n, rvec);
  rho_sq = r*r - n_dot_r*n_dot_r;

  if(rho_sq<0) {   // in case n is [0, 0, 1] and rho_sq is negative due to numerical error
    rho_sq = 0;
  }

  double del_sq = delta*delta;
  double rod_sq = rho_sq/del_sq;
  double td = exp(-rod_sq) * (C0 + rod_sq*(C2 + rod_sq*C4));
  dtd = -td/del_sq + exp(-rod_sq) * (C2 + 2*C4*rod_sq)/del_sq;

  return td;
}


// derivartive of dihedral angle func gij w.r.t rho, and atom positions
double RDPImplementation::dihedral(const int i, const int j,
    int const * const particleSpeciesCodes,
    VectorOfSizeDIM const * const coordinates,
    VectorOfSizeThreeInt const * const nearest3neigh,
    double const rhosq,
    double & d_drhosq, double d_dri[DIM], double d_drj[DIM],
    double d_drk1[DIM], double d_drk2[DIM], double d_drk3[DIM],
    double d_drl1[DIM], double d_drl2[DIM], double d_drl3[DIM]) const
{
  // get parameter
  int ispec = particleSpeciesCodes[i];
  int jspec = particleSpeciesCodes[j];
  double B = B_2D_[ispec][jspec];
  double eta = eta_2D_[ispec][jspec];
  double cut_rhosq = rhocutSq_2D_[ispec][jspec];

  // local vars
  double cos_kl[3][3];          // cos_omega_k1ijl1, cos_omega_k1ijl2 ...
  double d_dcos_kl[3][3];       // deriv of dihedral w.r.t to cos_omega_kijl
  double dcos_kl[3][3][4][DIM]; // 4 indicates k, i, j, l, e.g. dcoskl[0][1][0] means
                                // dcos_omega_k1ijl2 / drk


  // if larger than cutoff of rho, return 0
  if (rhosq >= cut_rhosq) {
    d_drhosq = 0;
    for (int dim=0; dim<DIM; dim++) {
      d_dri[dim] = 0;
      d_drj[dim] = 0;
      d_drk1[dim] = 0;
      d_drk2[dim] = 0;
      d_drk3[dim] = 0;
      d_drl1[dim] = 0;
      d_drl2[dim] = 0;
      d_drl3[dim] = 0;
    }
    double dihe = 0.0;
    return dihe;
  }

  // 3 neighs of atoms i and j
  int k[3];
  int l[3];
  for (int m=0; m<3; m++) {
    k[m] = nearest3neigh[i][m];
    l[m] = nearest3neigh[j][m];
  }

  // cos_omega_kijl and the derivatives w.r.t coordinates
  for (int m=0; m<3; m++) {
    for (int n=0; n<3; n++) {
      cos_kl[m][n] = deriv_cos_omega(
          coordinates[k[m]], coordinates[i], coordinates[j], coordinates[l[n]],
          dcos_kl[m][n][0], dcos_kl[m][n][1], dcos_kl[m][n][2], dcos_kl[m][n][3]);
    }
  }

  double epart1 = exp(-eta * cos_kl[0][0]*cos_kl[0][1]*cos_kl[0][2]);
  double epart2 = exp(-eta * cos_kl[1][0]*cos_kl[1][1]*cos_kl[1][2]);
  double epart3 = exp(-eta * cos_kl[2][0]*cos_kl[2][1]*cos_kl[2][2]);
  double  D2 = epart1 + epart2 + epart3;

  // cutoff function
  double d_drhosq_tap;
  double D0 = B * tap_rho(rhosq, cut_rhosq, d_drhosq_tap);

  // dihedral energy
  double dihe = D0*D2;

  // deriv of dihedral w.r.t rhosq
  d_drhosq = B*d_drhosq_tap*D2;

  // deriv of dihedral w.r.t cos_omega_kijl
  d_dcos_kl[0][0] = -D0* epart1 *eta *cos_kl[0][1]*cos_kl[0][2];
  d_dcos_kl[0][1] = -D0* epart1 *eta *cos_kl[0][0]*cos_kl[0][2];
  d_dcos_kl[0][2] = -D0* epart1 *eta *cos_kl[0][0]*cos_kl[0][1];
  d_dcos_kl[1][0] = -D0* epart2 *eta *cos_kl[1][1]*cos_kl[1][2];
  d_dcos_kl[1][1] = -D0* epart2 *eta *cos_kl[1][0]*cos_kl[1][2];
  d_dcos_kl[1][2] = -D0* epart2 *eta *cos_kl[1][0]*cos_kl[1][1];
  d_dcos_kl[2][0] = -D0* epart3 *eta *cos_kl[2][1]*cos_kl[2][2];
  d_dcos_kl[2][1] = -D0* epart3 *eta *cos_kl[2][0]*cos_kl[2][2];
  d_dcos_kl[2][2] = -D0* epart3 *eta *cos_kl[2][0]*cos_kl[2][1];

  // initialization to be zero and later add values
  for (int dim=0; dim<DIM; dim++) {
    d_drk1[dim] = 0.;
    d_drk2[dim] = 0.;
    d_drk3[dim] = 0.;
    d_dri[dim] = 0.;
    d_drj[dim] = 0.;
    d_drl1[dim] = 0.;
    d_drl2[dim] = 0.;
    d_drl3[dim] = 0.;
  }

  for (int m=0; m<3; m++) {
    for (int dim=0; dim<3; dim++) {
      d_drk1[dim] += d_dcos_kl[0][m]*dcos_kl[0][m][0][dim];
      d_drk2[dim] += d_dcos_kl[1][m]*dcos_kl[1][m][0][dim];
      d_drk3[dim] += d_dcos_kl[2][m]*dcos_kl[2][m][0][dim];
      d_drl1[dim] += d_dcos_kl[m][0]*dcos_kl[m][0][3][dim];
      d_drl2[dim] += d_dcos_kl[m][1]*dcos_kl[m][1][3][dim];
      d_drl3[dim] += d_dcos_kl[m][2]*dcos_kl[m][2][3][dim];
    }
    for (int n=0; n<3; n++) {
      for (int dim=0; dim<3; dim++) {
        d_dri[dim] += d_dcos_kl[m][n]*dcos_kl[m][n][1][dim];
        d_drj[dim] += d_dcos_kl[m][n]*dcos_kl[m][n][2][dim];
      }
    }
  }

  return dihe;
}



// compute cos(omega_kijl) and the derivateives
double RDPImplementation::deriv_cos_omega(double const * const rk,
    double const * const ri, double const * const rj, double const * const rl,
    double * const dcos_drk, double * const dcos_dri, double * const dcos_drj,
    double * const dcos_drl) const
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


  // ejik and derivatives (Note the dejik_dri ... returned are actually the transpose)
  deriv_cross(ri, rj, rk, ejik, dejik_dri, dejik_drj, dejik_drk);

  // flip sign because deriv_cross computes rij cross rik, here we need rji cross rik
  for (int m=0; m<DIM; m++) {
    ejik[m] = - ejik[m];
    for (int n=0; n<DIM; n++) {
      dejik_dri[m][n] = - dejik_dri[m][n];
      dejik_drj[m][n] = - dejik_drj[m][n];
      dejik_drk[m][n] = - dejik_drk[m][n];
    }
  }

  // eijl and derivatives
  deriv_cross(rj, ri, rl, eijl, deijl_drj, deijl_dri, deijl_drl);
  // flip sign
  for (int m=0; m<DIM; m++) {
    eijl[m] = - eijl[m];
    for (int n=0; n<DIM; n++) {
      deijl_drj[m][n] = - deijl_drj[m][n];
      deijl_dri[m][n] = - deijl_dri[m][n];
      deijl_drl[m][n] = - deijl_drl[m][n];
    }
  }

  // dcos_drk
  mat_dot_vec(dejik_drk, eijl, dcos_drk);
  // dcos_dri
  mat_dot_vec(dejik_dri, eijl, tmp1);
  mat_dot_vec(deijl_dri, ejik, tmp2);
  for (int m=0; m<DIM; m++) {
    dcos_dri[m] = tmp1[m] + tmp2[m];
  }
  // dcos_drj
  mat_dot_vec(dejik_drj, eijl, tmp1);
  mat_dot_vec(deijl_drj, ejik, tmp2);
  for (int m=0; m<DIM; m++) {
    dcos_drj[m] = tmp1[m] + tmp2[m];
  }
  // dcos drl
  mat_dot_vec(deijl_drl, ejik, dcos_drl);

  // cos_oemga_kijl
  double cos_omega = dot(ejik, eijl);

  return cos_omega;
}


// tap cutoff function
double RDPImplementation::tap(double r, double cutoff, double & dtap) const
{
  double t;
  double r_min = 0;

  if (r <= r_min) {
    t = 1;
    dtap = 0;
  } else {
    double roc = (r-r_min)/(cutoff-r_min);
    double roc_sq = roc*roc;
    t = roc_sq*roc_sq* (-35.0 + 84.0*roc + roc_sq* (-70.0 + 20.0*roc)) + 1;
    dtap = roc_sq*roc/(cutoff-r_min)* (-140.0 + 420.0*roc + roc_sq* (-420.0 + 140.0*roc));
  }

  return t;
}

// tap rho
double RDPImplementation::tap_rho(double rhosq, double cut_rhosq, double & drhosq) const
{
  double roc_sq;
  double roc;
  double t;

  roc_sq = rhosq/cut_rhosq;
  roc = sqrt(roc_sq);
  t = roc_sq*roc_sq* (-35.0 + 84.0*roc + roc_sq* (-70.0 + 20.0*roc)) + 1;
  // Note this dtap/drho_sq not dtap/drho
  drhosq = roc_sq/cut_rhosq* (-70.0 + 210.0*roc + roc_sq* (-210.0 + 70.0*roc));

  return t;
}


// helper functions

// dot product of two vector
double RDPImplementation::dot(double const * const x, double const * const y) const
{
  return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}

// matrix  product a vector, return a vector
void RDPImplementation::mat_dot_vec(VectorOfSizeDIM const * const X,
    double const * const y, double * const z) const
{
  int k;
  for (k=0; k<DIM; k++) {
    z[k] = X[k][0]*y[0] + X[k][1]*y[1] + X[k][2]*y[2];
  }
}

// Compute the normalized cross product of two vector rkl, rkm, and the derivates
//   w.r.t rk, rl, rm.
//  NOTE, the dcross_drk, dcross_drl, and dcross_drm is actually the transpose of the
//  actual one.
void RDPImplementation::deriv_cross(double const * const rk, double const * const rl,
    double const * const rm, double * const cross, VectorOfSizeDIM * const dcross_drk,
    VectorOfSizeDIM * const dcross_drl, VectorOfSizeDIM * const dcross_drm) const
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


  // get x = rkl and y = rkm
  for (i=0; i<DIM; i++) {
    x[i] = rl[i] - rk[i];
    y[i] = rm[i] - rk[i];
  }

  // cross product
  p[0] = x[1]*y[2] - x[2]*y[1];
  p[1] = x[2]*y[0] - x[0]*y[2];
  p[2] = x[0]*y[1] - x[1]*y[0];

  q = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);

  // normalized cross
  cross[0] = p[0]/q;
  cross[1] = p[1]/q;
  cross[2] = p[2]/q;

  // compute derivatives
  // derivative of inverse q (i.e. 1/q) w.r.t x and y
  q_cubic = q*q*q;
  d_invq_d_x0 = (           + p[1]*y[2] - p[2]*y[1])/q_cubic;
  d_invq_d_x1 = (-p[0]*y[2]             + p[2]*y[0])/q_cubic;
  d_invq_d_x2 = ( p[0]*y[1] - p[1]*y[0]            )/q_cubic;
  d_invq_d_y0 = (           - p[1]*x[2] + p[2]*x[1])/q_cubic;
  d_invq_d_y1 = ( p[0]*x[2]             - p[2]*x[0])/q_cubic;
  d_invq_d_y2 = (-p[0]*x[1] + p[1]*x[0]            )/q_cubic;

  // dcross/drl transposed
  dcross_drl[0][0] =           p[0]*d_invq_d_x0;
  dcross_drl[0][1] = -y[2]/q + p[1]*d_invq_d_x0;
  dcross_drl[0][2] =  y[1]/q + p[2]*d_invq_d_x0;

  dcross_drl[1][0] =  y[2]/q + p[0]*d_invq_d_x1;
  dcross_drl[1][1] =           p[1]*d_invq_d_x1;
  dcross_drl[1][2] = -y[0]/q + p[2]*d_invq_d_x1;

  dcross_drl[2][0] = -y[1]/q + p[0]*d_invq_d_x2;
  dcross_drl[2][1] =  y[0]/q + p[1]*d_invq_d_x2;
  dcross_drl[2][2] =           p[2]*d_invq_d_x2;

  // dcross/drm transposed
  dcross_drm[0][0] =           p[0]*d_invq_d_y0;
  dcross_drm[0][1] =  x[2]/q + p[1]*d_invq_d_y0;
  dcross_drm[0][2] = -x[1]/q + p[2]*d_invq_d_y0;

  dcross_drm[1][0] = -x[2]/q + p[0]*d_invq_d_y1;
  dcross_drm[1][1] =           p[1]*d_invq_d_y1;
  dcross_drm[1][2] =  x[0]/q + p[2]*d_invq_d_y1;

  dcross_drm[2][0] =  x[1]/q + p[0]*d_invq_d_y2;
  dcross_drm[2][1] = -x[0]/q + p[1]*d_invq_d_y2;
  dcross_drm[2][2] =           p[2]*d_invq_d_y2;

  // dcross/drk transposed
  for (i=0; i<DIM; i++) {
    for (j=0; j<DIM; j++) {
      dcross_drk[i][j] = - (dcross_drl[i][j] + dcross_drm[i][j]);
    }
  }

}



