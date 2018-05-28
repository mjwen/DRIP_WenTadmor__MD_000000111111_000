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
#include "helper.h"


// allocate memory and set pointers
void AllocateAndInitialize1DArray(double*& arrayPtr, int const extent)
{
  arrayPtr = new double[extent];
  for (int i = 0; i < extent; ++i) {
    arrayPtr[i] = 0.0;
  }
}

// deallocate memory
void Deallocate1DArray(double*& arrayPtr) {
	delete [] arrayPtr;
	// nullify pointer
	arrayPtr = 0;
}

// allocate memory and set pointers
void AllocateAndInitialize2DArray(double**& arrayPtr, int const extentZero,
		int const extentOne)
{
	arrayPtr = new double*[extentZero];
	arrayPtr[0] = new double[extentZero * extentOne];
	for (int i = 1; i < extentZero; ++i) {
		arrayPtr[i] = arrayPtr[i-1] + extentOne;
	}

	// initialize
	for (int i = 0; i < extentZero; ++i) {
		for (int j = 0; j < extentOne; ++j) {
			arrayPtr[i][j] = 0.0;
		}
	}
}

// deallocate memory
void Deallocate2DArray(double**& arrayPtr) {
	if (arrayPtr != 0) delete [] arrayPtr[0];
	delete [] arrayPtr;

	// nullify pointer
	arrayPtr = 0;
}

// allocate memory and set pointers
void AllocateAndInitialize3DArray(double***& arrayPtr, int const extentZero,
		int const extentOne, int const extentTwo)
{
  arrayPtr = new double**[extentZero];
  arrayPtr[0] = new double*[extentZero * extentOne];
  arrayPtr[0][0] = new double[extentZero * extentOne * extentTwo];

  for (int i = 1; i < extentZero; ++i) {
    arrayPtr[i] = arrayPtr[i-1] + extentOne;
    arrayPtr[i][0] = arrayPtr[i-1][0] + extentOne * extentTwo;
  }

  for (int i = 0; i < extentZero; ++i) {
    for (int j = 1; j < extentOne; ++j) {
      arrayPtr[i][j] = arrayPtr[i][j-1] + extentTwo;
    }
  }

  // initialize
  for (int i = 0; i < extentZero; ++i) {
    for (int j = 0; j < extentOne; ++j) {
      for (int k = 0; k < extentTwo; ++k) {
        arrayPtr[i][j][k] = 0.0;
      }
    }
  }
}

// deallocate memory
void Deallocate3DArray(double***& arrayPtr) {
  if (arrayPtr != 0) {
    if (arrayPtr[0] != 0) {
      delete [] arrayPtr[0][0];
    }
    delete [] arrayPtr[0];
  }
  delete [] arrayPtr;

	// nullify pointer
	arrayPtr = 0;
}


