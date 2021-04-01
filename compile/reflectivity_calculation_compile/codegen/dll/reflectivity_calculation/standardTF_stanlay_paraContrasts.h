//
// Non-Degree Granting Education License -- for use at non-degree
// granting, nonprofit, educational organizations only. Not for
// government, commercial, or other organizational use.
// File: standardTF_stanlay_paraContrasts.h
//
// MATLAB Coder version            : 5.0
// C/C++ source code generated on  : 24-Feb-2021 09:21:20
//
#ifndef STANDARDTF_STANLAY_PARACONTRASTS_H
#define STANDARDTF_STANLAY_PARACONTRASTS_H

// Include Files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "omp.h"
#include "reflectivity_calculation_types.h"
#define MAX_THREADS                    omp_get_max_threads()

// Function Declarations
extern void c_standardTF_stanlay_paraContra(const coder::array<double, 2U>
  &resample, double numberOfContrasts, const coder::array<char, 2U> &geometry,
  const coder::array<cell_wrap_0, 2U> &repeatLayers, const coder::array<double,
  2U> &cBacks, const coder::array<double, 2U> &cShifts, const coder::array<
  double, 2U> &cScales, const coder::array<double, 2U> &cNbas, const coder::
  array<double, 2U> &cNbss, const coder::array<double, 2U> &cRes, const coder::
  array<double, 2U> &backs, const coder::array<double, 2U> &shifts, const coder::
  array<double, 2U> &sf, const coder::array<double, 2U> &nba, const coder::array<
  double, 2U> &nbs, const coder::array<double, 2U> &res, const coder::array<
  double, 2U> &dataPresent, const coder::array<cell_wrap_1, 2U> &allData, const
  coder::array<cell_wrap_0, 2U> &dataLimits, const coder::array<cell_wrap_0, 2U>
  &simLimits, double nParams, const coder::array<double, 2U> &params, const
  coder::array<cell_wrap_11, 2U> &contrastLayers, const coder::array<cell_wrap_3,
  1U> &layersDetails, const coder::array<double, 2U> &backsType, double calcSld,
  coder::array<double, 1U> &outSsubs, coder::array<double, 1U> &backgs, coder::
  array<double, 1U> &qshifts, coder::array<double, 1U> &sfs, coder::array<double,
  1U> &nbas, coder::array<double, 1U> &nbss, coder::array<double, 1U> &resols,
  coder::array<double, 1U> &chis, coder::array<cell_wrap_6, 1U> &reflectivity,
  coder::array<cell_wrap_6, 1U> &Simulation, coder::array<cell_wrap_1, 1U>
  &shifted_data, coder::array<cell_wrap_1, 1U> &layerSlds, coder::array<
  cell_wrap_11, 1U> &sldProfiles, coder::array<cell_wrap_10, 1U> &allLayers,
  coder::array<double, 1U> &allRoughs);

#endif

//
// File trailer for standardTF_stanlay_paraContrasts.h
//
// [EOF]
//