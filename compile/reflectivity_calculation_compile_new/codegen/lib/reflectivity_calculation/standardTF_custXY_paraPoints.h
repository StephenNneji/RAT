//
// Non-Degree Granting Education License -- for use at non-degree
// granting, nonprofit, educational organizations only. Not for
// government, commercial, or other organizational use.
//
// standardTF_custXY_paraPoints.h
//
// Code generation for function 'standardTF_custXY_paraPoints'
//

#ifndef STANDARDTF_CUSTXY_PARAPOINTS_H
#define STANDARDTF_CUSTXY_PARAPOINTS_H

// Include files
#include "rtwtypes.h"
#include "coder_array.h"
#include "omp.h"
#include <cstddef>
#include <cstdlib>

// Type Declarations
struct struct0_T;

struct cell_14;

struct cell_wrap_10;

struct cell_wrap_1;

struct cell_wrap_8;

struct cell_wrap_23;

// Function Declarations
void standardTF_custXY_paraPoints(
    const struct0_T *problemDef, const cell_14 *problemDef_cells,
    coder::array<double, 1U> &outSsubs, coder::array<double, 1U> &backgs,
    coder::array<double, 1U> &qshifts, coder::array<double, 1U> &sfs,
    coder::array<double, 1U> &nbas, coder::array<double, 1U> &nbss,
    coder::array<double, 1U> &resols, coder::array<double, 1U> &chis,
    coder::array<cell_wrap_10, 1U> &reflectivity,
    coder::array<cell_wrap_10, 1U> &Simulation,
    coder::array<cell_wrap_1, 1U> &shifted_data,
    coder::array<cell_wrap_8, 1U> &layerSlds,
    coder::array<cell_wrap_1, 1U> &sldProfiles,
    coder::array<cell_wrap_23, 1U> &allLayers,
    coder::array<double, 1U> &allRoughs);

#endif
// End of code generation (standardTF_custXY_paraPoints.h)