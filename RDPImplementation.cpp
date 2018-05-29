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


#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>

#include "KIM_Numbering.hpp"
#include "KIM_LanguageName.hpp"
#include "KIM_SpeciesName.hpp"
#include "KIM_SupportStatus.hpp"
#include "KIM_ComputeArgumentName.hpp"
#include "KIM_ComputeCallbackName.hpp"
#include "RDPImplementation.hpp"
#include "helper.hpp"

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
  C0_(0),
  C2_(0),
  C4_(0),
  C_(0),
  delta_(0),
  lambda_(0),
  A_(0),
  z0_(0),
  B_(0),
  eta_(0),
  influenceDistance_(0.0),
  cutoffSq_2D_(0),
  C0_2D_(0),
  C2_2D_(0),
  C4_2D_(0),
  C_2D_(0),
  delta_2D_(0),
  lambda_2D_(0),
  A_2D_(0),
  z0_2D_(0),
  B_2D_(0),
  eta_2D_(0),
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
int RDPImplementation::Refresh(
    KIM::ModelRefresh * const modelRefresh)
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
  int numberOfLinesRead;
  int endOfFileFlag = 0;
  char spec1[MAXLINE], nextLine[MAXLINE];
  int iIndex, jIndex , indx;
  double next_C0, next_C2, next_C4, next_C, next_delta, next_lambda, next_A;
  double next_z0, next_B, next_eta, next_cutoff;


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


  // keep track of known species
  std::map<KIM::SpeciesName const, int, KIM::SPECIES_NAME::Comparator>
      modelSpeciesMap;
  int index = 0;   // species code integer code starting from 0

  // Read and process data lines
  numberOfLinesRead = 0;
  getNextDataLine(parameterFilePointers[0], nextLine, MAXLINE, &endOfFileFlag);
  while (endOfFileFlag == 0)
  {
    ier = sscanf(nextLine, "%s %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg",
                 &spec1, &next_C0, &next_C2, &next_C4, &next_C, &next_delta,
                 &next_lambda, &next_A, &next_z0, &next_B, &next_eta, &next_cutoff);
    if (ier != 12) {
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
void RDPImplementation::AllocateFreeParameterMemory()
{ // allocate memory for data
  AllocateAndInitialize1DArray<double> (cutoff_, numberUniqueSpeciesPairs_);
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
          KIM::SUPPORT_STATUS::optional)
      || modelComputeArgumentsCreate->SetArgumentSupportStatus(
          KIM::COMPUTE_ARGUMENT_NAME::partialForces,
          KIM::SUPPORT_STATUS::optional)
      || modelComputeArgumentsCreate->SetArgumentSupportStatus(
          KIM::COMPUTE_ARGUMENT_NAME::partialParticleEnergy,
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
  ier = modelDriverCreate->SetParameterPointer(numberUniqueSpeciesPairs_, C0_, "C0")
     || modelDriverCreate->SetParameterPointer(numberUniqueSpeciesPairs_, C2_, "C2")
     || modelDriverCreate->SetParameterPointer(numberUniqueSpeciesPairs_, C4_, "C4")
     || modelDriverCreate->SetParameterPointer(numberUniqueSpeciesPairs_, C_, "C")
     || modelDriverCreate->SetParameterPointer(numberUniqueSpeciesPairs_, delta_, "delta")
     || modelDriverCreate->SetParameterPointer(numberUniqueSpeciesPairs_, lambda_, "lambda")
     || modelDriverCreate->SetParameterPointer(numberUniqueSpeciesPairs_, A_, "A")
     || modelDriverCreate->SetParameterPointer(numberUniqueSpeciesPairs_, z0_, "z0")
     || modelDriverCreate->SetParameterPointer(numberUniqueSpeciesPairs_, B_, "B")
     || modelDriverCreate->SetParameterPointer(numberUniqueSpeciesPairs_, eta_, "eta")
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

// create layers and find the nearest 3 neighbors of each atom
int RDPImplementation::create_layers( KIM::ModelCompute const * const modelCompute,
    KIM::ModelComputeArguments const * const modelComputeArguments,
    VectorOfSizeDIM const * const coordinates,
    int * const in_layer, // atoms in which layer , -1: not classified in any layer
    VectorOfSize3Int * const nearest3neigh)
{


  int const nAtoms = cachedNumberOfParticles_;
  int ier = true;

  // no atoms, which can happen when doing domain decompisition
  if (nAtoms <= 0) return ier;

  // cutoff to separate atoms into layers (may need to be adjusted)
  double cutsq_layer = (0.72*3.44)*(0.72*3.44);

  // layer contains which atoms
  int* layer;
  AllocateAndInitialize1DArray<int> (layer, nAtoms);

  // init atoms layer to -1 (not in any layer)
  for (k=0; k<nAtoms; k++) {
    in_layer[k] = -1;
  }

  int nlayers = 1;
  int nremain = nAtoms; // number of atoms not included in any layer

  while(true) {  // to create all layers

    // init all atoms in current layer have -1: not classified in any layer
    // We do not always have nAtoms (except the very first loop), we allocate layer
    // to the size of nAtom just because we do not want to reallocate it each time.
    for (k=0; k<nAtoms; k++) {
      layer[k] = -1;
    }

    // find an atom not incldued in any layer and start with it
    int currentLayer = nlayers - 1;
    for (int k=0; k<nAtoms; k++) {
      if (in_layer[k] == -1) {
        in_layer[k] = currentLayer;
        layer[0] = k; // first atom in current layer
        break;
      }
    }

    int nin = 1; // number of atoms in current layer
    int ii = 0;
    while(true) { // to find all atoms in currentLayer

      int i = layer[ii];

      // get neighbors of atom i
      int num_nei = 0;
      int const * ilist = 0;
      ier = modelComputeArguments->GetNeighborList(0, i, &num_nei, &ilist);
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

        j = ilist[jj];

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
  for (int i=0; i<nAtoms; i++) {
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
double RDPImplementation::calc_attractive(const int i, const int j,
    const int iSpecies, const int jSpecies, double const * const rij,
    const double r, VectorOfSizeDIM* const forces)
{
  // compute flags
  int const comp_forces = cachedIsComputeForces_;
  // compute params
  double const z0 = z0_2D_[iSpecies][jSpecies];

  double roz0_sq = r*r/(z0*z0);
  double dtp;
  double tp = tap(buffer, i, j, r, dtp);
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


// compute normal and its derivatives w.r.t atom ri, and its 3 nearest neighs k1, k2, k3
void RDPImplementation::normal(const int i, VectorOfSizeDIM const * const coordinates,
    VectorOfSize3Int const * const nearest3neigh,
    int & k1, int & k2, int & k3, double * const normal,
    VectorOfSizeDIM * const dn_dri, VectorOfSizeDIM * const dn_drk1,
    VectorOfSizeDIM * const dn_drk2, VectorOfSizeDIM * const dn_drk3)
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


// tap cutoff function
double RDPImplementation::tap(int i, int j, double r, double & dtap)
{
  double t;
  double r_min = 0;

  double cutoff = sqrt(cutoffSq_2D_[iSpecies][jSpecies]);

  if (r <= r_min) {
    t = 1;
    *dtap = 0;
  } else {
    double roc = (r-r_min)/(cutoff-r_min);
    double roc_sq = roc*roc;
    t = roc_sq*roc_sq* (-35.0 + 84.0*roc + roc_sq* (-70.0 + 20.0*roc)) + 1;
    dtap = roc_sq*roc/(cutoff-r_min)* (-140.0 + 420.0*roc + roc_sq* (-420.0 + 140.0*roc));
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


// Compute the normalized cross product of two vector rkl, rkm, and the derivates
//   w.r.t rk, rl, rm.
//  NOTE, the dcross_drk, dcross_drl, and dcross_drm is actually the transpose of the
//  actual one.
void RDPImplementation::deriv_cross(double const * const rk, double const * const rl,
    double const * const rm, double * const cross, VectorOfSizeDIM * const dcross_drk,
    VectorOfSizeDIM * const dcross_drl, VectorOfSizeDIM * const dcross_drm)
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


