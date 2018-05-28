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


#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>

#include "RDP.hpp"
#include "RDPImplementation.hpp"


//==============================================================================
//
// This is the standard interface to KIM Model Drivers
//
//==============================================================================

//******************************************************************************
extern "C"
{
int model_driver_create(
    KIM::ModelDriverCreate * const modelDriverCreate,
    KIM::LengthUnit const requestedLengthUnit,
    KIM::EnergyUnit const requestedEnergyUnit,
    KIM::ChargeUnit const requestedChargeUnit,
    KIM::TemperatureUnit const requestedTemperatureUnit,
    KIM::TimeUnit const requestedTimeUnit)
{
  int ier;

  // read input files, convert units if needed, compute
  // interpolation coefficients, set cutoff, and publish parameters
  RDP* const modelObject
      = new RDP(modelDriverCreate,
                            requestedLengthUnit,
                            requestedEnergyUnit,
                            requestedChargeUnit,
                            requestedTemperatureUnit,
                            requestedTimeUnit,
                            &ier);
  if (ier)
  {
    // constructor already reported the error
    delete modelObject;
    return ier;
  }

  // register pointer to RDP object in KIM object
  modelDriverCreate->SetModelBufferPointer(
      static_cast<void*>(modelObject));

  // everything is good
  ier = false;
  return ier;
}
}

//==============================================================================
//
// Implementation of RDP public wrapper functions
//
//==============================================================================

//******************************************************************************
RDP::RDP(
    KIM::ModelDriverCreate * const modelDriverCreate,
    KIM::LengthUnit const requestedLengthUnit,
    KIM::EnergyUnit const requestedEnergyUnit,
    KIM::ChargeUnit const requestedChargeUnit,
    KIM::TemperatureUnit const requestedTemperatureUnit,
    KIM::TimeUnit const requestedTimeUnit,
    int* const ier)
{
  implementation_ = new RDPImplementation(
      modelDriverCreate,
      requestedLengthUnit,
      requestedEnergyUnit,
      requestedChargeUnit,
      requestedTemperatureUnit,
      requestedTimeUnit,
      ier);
}

//******************************************************************************
RDP::~RDP()
{
  delete implementation_;
}

//******************************************************************************
// static member function
int RDP::Destroy(KIM::ModelDestroy * const modelDestroy)
{
  RDP * modelObject;
  modelDestroy->GetModelBufferPointer(reinterpret_cast<void**>(&modelObject));

  if (modelObject != NULL)
  {
    // delete object itself
    delete modelObject;
  }

  // everything is good
  return false;
}

//******************************************************************************
// static member function
int RDP::Refresh(
    KIM::ModelRefresh * const modelRefresh)
{
  RDP * modelObject;
  modelRefresh->GetModelBufferPointer(
      reinterpret_cast<void**>(&modelObject));

  return modelObject->implementation_->Refresh(modelRefresh);
}

//******************************************************************************
// static member function
int RDP::Compute(
    KIM::ModelCompute const * const modelCompute,
    KIM::ModelComputeArguments const * const modelComputeArguments)
{
  RDP * modelObject;
  modelCompute->GetModelBufferPointer(reinterpret_cast<void**>(&modelObject));

  return modelObject->implementation_->Compute(modelCompute, modelComputeArguments);

}

//******************************************************************************
// static member function
int RDP::ComputeArgumentsCreate(
    KIM::ModelCompute const * const modelCompute,
    KIM::ModelComputeArgumentsCreate * const modelComputeArgumentsCreate)
{
  RDP * modelObject;
  modelCompute->GetModelBufferPointer(reinterpret_cast<void**>(&modelObject));

  return modelObject->implementation_
      ->ComputeArgumentsCreate(modelComputeArgumentsCreate);
}

//******************************************************************************
// static member function
int RDP::ComputeArgumentsDestroy(
    KIM::ModelCompute const * const modelCompute,
    KIM::ModelComputeArgumentsDestroy * const modelComputeArgumentsDestroy)
{
  RDP * modelObject;
  modelCompute->GetModelBufferPointer(reinterpret_cast<void**>(&modelObject));

  return modelObject->implementation_
      ->ComputeArgumentsDestroy(modelComputeArgumentsDestroy);
}

