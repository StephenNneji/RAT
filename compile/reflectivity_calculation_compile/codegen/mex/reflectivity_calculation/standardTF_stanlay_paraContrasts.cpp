//
//  Non-Degree Granting Education License -- for use at non-degree
//  granting, nonprofit, educational organizations only. Not for
//  government, commercial, or other organizational use.
//
//  standardTF_stanlay_paraContrasts.cpp
//
//  Code generation for function 'standardTF_stanlay_paraContrasts'
//


// Include files
#include "standardTF_stanlay_paraContrasts.h"
#include "allocateLayersForContrast.h"
#include "allocateParamsToLayers.h"
#include "backSort.h"
#include "groupLayers_Mod.h"
#include "mwmathutil.h"
#include "reflectivity_calculation.h"
#include "reflectivity_calculation_data.h"
#include "rt_nonfinite.h"
#include "standardTF_layers_core.h"
#include "standardTF_stanlay_paraAll.h"
#include "sum.h"

// Variable Definitions
static emlrtRSInfo mg_emlrtRSI = { 58, // lineNo
  "standardTF_stanlay_paraContrasts",  // fcnName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m"// pathName 
};

static emlrtRSInfo ng_emlrtRSI = { 67, // lineNo
  "standardTF_stanlay_paraContrasts",  // fcnName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m"// pathName 
};

static emlrtRSInfo og_emlrtRSI = { 72, // lineNo
  "standardTF_stanlay_paraContrasts",  // fcnName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m"// pathName 
};

static emlrtRSInfo pg_emlrtRSI = { 91, // lineNo
  "standardTF_stanlay_paraContrasts",  // fcnName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m"// pathName 
};

static emlrtRTEInfo ac_emlrtRTEI = { 39,// lineNo
  9,                                   // colNo
  "standardTF_stanlay_paraContrasts",  // fName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m"// pName 
};

static emlrtRTEInfo bc_emlrtRTEI = { 44,// lineNo
  9,                                   // colNo
  "standardTF_stanlay_paraContrasts",  // fName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m"// pName 
};

static emlrtRTEInfo cc_emlrtRTEI = { 49,// lineNo
  9,                                   // colNo
  "standardTF_stanlay_paraContrasts",  // fName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m"// pName 
};

static emlrtDCInfo pb_emlrtDCI = { 61, // lineNo
  14,                                  // colNo
  "standardTF_stanlay_paraContrasts",  // fName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m",// pName 
  1                                    // checkKind
};

static emlrtRTEInfo dc_emlrtRTEI = { 61,// lineNo
  12,                                  // colNo
  "standardTF_stanlay_paraContrasts",  // fName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m"// pName 
};

static emlrtBCInfo nh_emlrtBCI = { -1, // iFirst
  -1,                                  // iLast
  76,                                  // lineNo
  24,                                  // colNo
  "params",                            // aName
  "standardTF_stanlay_paraContrasts",  // fName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m",// pName 
  0                                    // checkKind
};

static emlrtBCInfo oh_emlrtBCI = { -1, // iFirst
  -1,                                  // iLast
  78,                                  // lineNo
  29,                                  // colNo
  "resample",                          // aName
  "standardTF_stanlay_paraContrasts",  // fName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m",// pName 
  0                                    // checkKind
};

static emlrtBCInfo ph_emlrtBCI = { -1, // iFirst
  -1,                                  // iLast
  81,                                  // lineNo
  35,                                  // colNo
  "dataPresent",                       // aName
  "standardTF_stanlay_paraContrasts",  // fName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m",// pName 
  0                                    // checkKind
};

static emlrtBCInfo qh_emlrtBCI = { -1, // iFirst
  -1,                                  // iLast
  84,                                  // lineNo
  31,                                  // colNo
  "backsType",                         // aName
  "standardTF_stanlay_paraContrasts",  // fName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m",// pName 
  0                                    // checkKind
};

static emlrtBCInfo rh_emlrtBCI = { -1, // iFirst
  -1,                                  // iLast
  77,                                  // lineNo
  37,                                  // colNo
  "repeatLayers",                      // aName
  "standardTF_stanlay_paraContrasts",  // fName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m",// pName 
  0                                    // checkKind
};

static emlrtBCInfo sh_emlrtBCI = { -1, // iFirst
  -1,                                  // iLast
  80,                                  // lineNo
  24,                                  // colNo
  "allData",                           // aName
  "standardTF_stanlay_paraContrasts",  // fName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m",// pName 
  0                                    // checkKind
};

static emlrtBCInfo th_emlrtBCI = { -1, // iFirst
  -1,                                  // iLast
  82,                                  // lineNo
  33,                                  // colNo
  "dataLimits",                        // aName
  "standardTF_stanlay_paraContrasts",  // fName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m",// pName 
  0                                    // checkKind
};

static emlrtBCInfo uh_emlrtBCI = { -1, // iFirst
  -1,                                  // iLast
  83,                                  // lineNo
  31,                                  // colNo
  "simLimits",                         // aName
  "standardTF_stanlay_paraContrasts",  // fName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m",// pName 
  0                                    // checkKind
};

static emlrtDCInfo qb_emlrtDCI = { 26, // lineNo
  17,                                  // colNo
  "standardTF_stanlay_paraContrasts",  // fName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m",// pName 
  1                                    // checkKind
};

static emlrtDCInfo rb_emlrtDCI = { 25, // lineNo
  16,                                  // colNo
  "standardTF_stanlay_paraContrasts",  // fName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m",// pName 
  1                                    // checkKind
};

static emlrtDCInfo sb_emlrtDCI = { 25, // lineNo
  16,                                  // colNo
  "standardTF_stanlay_paraContrasts",  // fName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m",// pName 
  4                                    // checkKind
};

static emlrtBCInfo vh_emlrtBCI = { -1, // iFirst
  -1,                                  // iLast
  40,                                  // lineNo
  5,                                   // colNo
  "reflectivity",                      // aName
  "standardTF_stanlay_paraContrasts",  // fName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m",// pName 
  0                                    // checkKind
};

static emlrtBCInfo wh_emlrtBCI = { -1, // iFirst
  -1,                                  // iLast
  45,                                  // lineNo
  5,                                   // colNo
  "Simulation",                        // aName
  "standardTF_stanlay_paraContrasts",  // fName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m",// pName 
  0                                    // checkKind
};

static emlrtBCInfo xh_emlrtBCI = { -1, // iFirst
  -1,                                  // iLast
  50,                                  // lineNo
  5,                                   // colNo
  "allLayers",                         // aName
  "standardTF_stanlay_paraContrasts",  // fName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m",// pName 
  0                                    // checkKind
};

static emlrtBCInfo yh_emlrtBCI = { -1, // iFirst
  -1,                                  // iLast
  67,                                  // lineNo
  77,                                  // colNo
  "cBacks",                            // aName
  "standardTF_stanlay_paraContrasts",  // fName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m",// pName 
  0                                    // checkKind
};

static emlrtBCInfo ai_emlrtBCI = { -1, // iFirst
  -1,                                  // iLast
  67,                                  // lineNo
  87,                                  // colNo
  "cShifts",                           // aName
  "standardTF_stanlay_paraContrasts",  // fName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m",// pName 
  0                                    // checkKind
};

static emlrtBCInfo bi_emlrtBCI = { -1, // iFirst
  -1,                                  // iLast
  67,                                  // lineNo
  98,                                  // colNo
  "cScales",                           // aName
  "standardTF_stanlay_paraContrasts",  // fName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m",// pName 
  0                                    // checkKind
};

static emlrtBCInfo ci_emlrtBCI = { -1, // iFirst
  -1,                                  // iLast
  67,                                  // lineNo
  109,                                 // colNo
  "cNbas",                             // aName
  "standardTF_stanlay_paraContrasts",  // fName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m",// pName 
  0                                    // checkKind
};

static emlrtBCInfo di_emlrtBCI = { -1, // iFirst
  -1,                                  // iLast
  67,                                  // lineNo
  118,                                 // colNo
  "cNbss",                             // aName
  "standardTF_stanlay_paraContrasts",  // fName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m",// pName 
  0                                    // checkKind
};

static emlrtBCInfo ei_emlrtBCI = { -1, // iFirst
  -1,                                  // iLast
  67,                                  // lineNo
  127,                                 // colNo
  "cRes",                              // aName
  "standardTF_stanlay_paraContrasts",  // fName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m",// pName 
  0                                    // checkKind
};

static emlrtBCInfo fi_emlrtBCI = { -1, // iFirst
  -1,                                  // iLast
  72,                                  // lineNo
  52,                                  // colNo
  "contrastLayers",                    // aName
  "standardTF_stanlay_paraContrasts",  // fName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m",// pName 
  0                                    // checkKind
};

static emlrtBCInfo gi_emlrtBCI = { -1, // iFirst
  -1,                                  // iLast
  101,                                 // lineNo
  5,                                   // colNo
  "outSsubs",                          // aName
  "standardTF_stanlay_paraContrasts",  // fName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m",// pName 
  0                                    // checkKind
};

static emlrtBCInfo hi_emlrtBCI = { -1, // iFirst
  -1,                                  // iLast
  102,                                 // lineNo
  5,                                   // colNo
  "sldProfiles",                       // aName
  "standardTF_stanlay_paraContrasts",  // fName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m",// pName 
  0                                    // checkKind
};

static emlrtBCInfo ii_emlrtBCI = { -1, // iFirst
  -1,                                  // iLast
  103,                                 // lineNo
  5,                                   // colNo
  "reflectivity",                      // aName
  "standardTF_stanlay_paraContrasts",  // fName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m",// pName 
  0                                    // checkKind
};

static emlrtBCInfo ji_emlrtBCI = { -1, // iFirst
  -1,                                  // iLast
  104,                                 // lineNo
  5,                                   // colNo
  "Simulation",                        // aName
  "standardTF_stanlay_paraContrasts",  // fName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m",// pName 
  0                                    // checkKind
};

static emlrtBCInfo ki_emlrtBCI = { -1, // iFirst
  -1,                                  // iLast
  105,                                 // lineNo
  5,                                   // colNo
  "shifted_data",                      // aName
  "standardTF_stanlay_paraContrasts",  // fName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m",// pName 
  0                                    // checkKind
};

static emlrtBCInfo li_emlrtBCI = { -1, // iFirst
  -1,                                  // iLast
  106,                                 // lineNo
  5,                                   // colNo
  "layerSlds",                         // aName
  "standardTF_stanlay_paraContrasts",  // fName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m",// pName 
  0                                    // checkKind
};

static emlrtBCInfo mi_emlrtBCI = { -1, // iFirst
  -1,                                  // iLast
  107,                                 // lineNo
  5,                                   // colNo
  "chis",                              // aName
  "standardTF_stanlay_paraContrasts",  // fName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m",// pName 
  0                                    // checkKind
};

static emlrtBCInfo ni_emlrtBCI = { -1, // iFirst
  -1,                                  // iLast
  108,                                 // lineNo
  5,                                   // colNo
  "backgs",                            // aName
  "standardTF_stanlay_paraContrasts",  // fName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m",// pName 
  0                                    // checkKind
};

static emlrtBCInfo oi_emlrtBCI = { -1, // iFirst
  -1,                                  // iLast
  109,                                 // lineNo
  5,                                   // colNo
  "qshifts",                           // aName
  "standardTF_stanlay_paraContrasts",  // fName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m",// pName 
  0                                    // checkKind
};

static emlrtBCInfo pi_emlrtBCI = { -1, // iFirst
  -1,                                  // iLast
  110,                                 // lineNo
  5,                                   // colNo
  "sfs",                               // aName
  "standardTF_stanlay_paraContrasts",  // fName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m",// pName 
  0                                    // checkKind
};

static emlrtBCInfo qi_emlrtBCI = { -1, // iFirst
  -1,                                  // iLast
  111,                                 // lineNo
  5,                                   // colNo
  "nbas",                              // aName
  "standardTF_stanlay_paraContrasts",  // fName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m",// pName 
  0                                    // checkKind
};

static emlrtBCInfo ri_emlrtBCI = { -1, // iFirst
  -1,                                  // iLast
  112,                                 // lineNo
  5,                                   // colNo
  "nbss",                              // aName
  "standardTF_stanlay_paraContrasts",  // fName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m",// pName 
  0                                    // checkKind
};

static emlrtBCInfo si_emlrtBCI = { -1, // iFirst
  -1,                                  // iLast
  113,                                 // lineNo
  5,                                   // colNo
  "resols",                            // aName
  "standardTF_stanlay_paraContrasts",  // fName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m",// pName 
  0                                    // checkKind
};

static emlrtBCInfo ti_emlrtBCI = { -1, // iFirst
  -1,                                  // iLast
  114,                                 // lineNo
  5,                                   // colNo
  "allRoughs",                         // aName
  "standardTF_stanlay_paraContrasts",  // fName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m",// pName 
  0                                    // checkKind
};

static emlrtRTEInfo yi_emlrtRTEI = { 25,// lineNo
  10,                                  // colNo
  "standardTF_stanlay_paraContrasts",  // fName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m"// pName 
};

static emlrtRTEInfo aj_emlrtRTEI = { 4,// lineNo
  18,                                  // colNo
  "standardTF_stanlay_paraContrasts",  // fName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m"// pName 
};

static emlrtRTEInfo bj_emlrtRTEI = { 40,// lineNo
  5,                                   // colNo
  "standardTF_stanlay_paraContrasts",  // fName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m"// pName 
};

static emlrtRTEInfo cj_emlrtRTEI = { 45,// lineNo
  5,                                   // colNo
  "standardTF_stanlay_paraContrasts",  // fName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m"// pName 
};

static emlrtRTEInfo dj_emlrtRTEI = { 50,// lineNo
  5,                                   // colNo
  "standardTF_stanlay_paraContrasts",  // fName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m"// pName 
};

static emlrtRTEInfo ej_emlrtRTEI = { 32,// lineNo
  12,                                  // colNo
  "standardTF_stanlay_paraContrasts",  // fName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m"// pName 
};

static emlrtRTEInfo fj_emlrtRTEI = { 102,// lineNo
  5,                                   // colNo
  "standardTF_stanlay_paraContrasts",  // fName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m"// pName 
};

static emlrtRTEInfo gj_emlrtRTEI = { 103,// lineNo
  5,                                   // colNo
  "standardTF_stanlay_paraContrasts",  // fName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m"// pName 
};

static emlrtRTEInfo hj_emlrtRTEI = { 104,// lineNo
  5,                                   // colNo
  "standardTF_stanlay_paraContrasts",  // fName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m"// pName 
};

static emlrtRTEInfo ij_emlrtRTEI = { 105,// lineNo
  5,                                   // colNo
  "standardTF_stanlay_paraContrasts",  // fName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m"// pName 
};

static emlrtRTEInfo jj_emlrtRTEI = { 106,// lineNo
  5,                                   // colNo
  "standardTF_stanlay_paraContrasts",  // fName
  "/home/arwel/Documents/RascalDev/RAT/targetFunctions/standard_TF/standardTF_stanLay/standardTF_stanlay_paraContrasts.m"// pName 
};

// Function Definitions
void c_standardTF_stanlay_paraContra(const emlrtStack *sp, const struct0_T
  *problemDef, const cell_15 *problemDef_cells, const struct2_T *controls, coder::
  array<real_T, 1U> &outSsubs, coder::array<real_T, 1U> &backgs, coder::array<
  real_T, 1U> &qshifts, coder::array<real_T, 1U> &sfs, coder::array<real_T, 1U>
  &nbas, coder::array<real_T, 1U> &nbss, coder::array<real_T, 1U> &resols, coder::
  array<real_T, 1U> &chis, coder::array<cell_wrap_7, 1U> &reflectivity, coder::
  array<cell_wrap_7, 1U> &Simulation, coder::array<cell_wrap_1, 1U>
  &shifted_data, coder::array<cell_wrap_1, 1U> &layerSlds, coder::array<
  cell_wrap_12, 1U> &sldProfiles, coder::array<cell_wrap_11, 1U> &allLayers,
  coder::array<real_T, 1U> &allRoughs)
{
  coder::array<cell_wrap_0, 2U> repeatLayers;
  int32_T loop_ub;
  int32_T i;
  coder::array<cell_wrap_1, 2U> allData;
  coder::array<cell_wrap_0, 2U> dataLimits;
  coder::array<cell_wrap_0, 2U> simLimits;
  coder::array<cell_wrap_12, 2U> contrastLayers;
  coder::array<char_T, 2U> geometry;
  coder::array<real_T, 2U> cBacks;
  coder::array<real_T, 2U> cShifts;
  coder::array<real_T, 2U> cScales;
  coder::array<real_T, 2U> cNbas;
  coder::array<real_T, 2U> cNbss;
  coder::array<real_T, 2U> cRes;
  coder::array<real_T, 2U> backs;
  coder::array<real_T, 2U> shifts;
  coder::array<real_T, 2U> sf;
  coder::array<real_T, 2U> nba;
  coder::array<real_T, 2U> nbs;
  coder::array<real_T, 2U> res;
  coder::array<real_T, 2U> dataPresent;
  int32_T nParams;
  coder::array<real_T, 2U> params;
  coder::array<real_T, 2U> resample;
  coder::array<real_T, 2U> backsType;
  real_T calcSld;
  int32_T b_i;
  int32_T i1;
  int32_T i2;
  coder::array<cell_wrap_19, 2U> outParameterisedLayers;
  int32_T c_i;
  coder::array<real_T, 2U> thisContrastLayers;
  coder::array<real_T, 2U> sldProfile;
  coder::array<real_T, 2U> reflect;
  coder::array<real_T, 2U> Simul;
  coder::array<real_T, 2U> shifted_dat;
  coder::array<real_T, 2U> thisData;
  real_T thisSsubs;
  real_T thisChiSquared;
  real_T thisResol;
  real_T thisNbs;
  real_T thisNba;
  real_T thisSf;
  real_T thisQshift;
  real_T thisBackground;
  int32_T i3;
  int32_T b_loop_ub;
  int32_T i4;
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  boolean_T emlrtHadParallelError = false;
  st.prev = sp;
  st.tls = sp->tls;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);

  //  Extract individual cell arrays
  repeatLayers.set_size((&ci_emlrtRTEI), sp, 1, problemDef_cells->f1.size(1));
  loop_ub = problemDef_cells->f1.size(0) * problemDef_cells->f1.size(1);
  for (i = 0; i < loop_ub; i++) {
    repeatLayers[i] = problemDef_cells->f1[i];
  }

  allData.set_size((&di_emlrtRTEI), sp, 1, problemDef_cells->f2.size(1));
  loop_ub = problemDef_cells->f2.size(0) * problemDef_cells->f2.size(1);
  for (i = 0; i < loop_ub; i++) {
    allData[i] = problemDef_cells->f2[i];
  }

  dataLimits.set_size((&ei_emlrtRTEI), sp, 1, problemDef_cells->f3.size(1));
  loop_ub = problemDef_cells->f3.size(0) * problemDef_cells->f3.size(1);
  for (i = 0; i < loop_ub; i++) {
    dataLimits[i] = problemDef_cells->f3[i];
  }

  simLimits.set_size((&fi_emlrtRTEI), sp, 1, problemDef_cells->f4.size(1));
  loop_ub = problemDef_cells->f4.size(0) * problemDef_cells->f4.size(1);
  for (i = 0; i < loop_ub; i++) {
    simLimits[i] = problemDef_cells->f4[i];
  }

  contrastLayers.set_size((&gi_emlrtRTEI), sp, 1, problemDef_cells->f5.size(1));
  loop_ub = problemDef_cells->f5.size(0) * problemDef_cells->f5.size(1);
  for (i = 0; i < loop_ub; i++) {
    contrastLayers[i] = problemDef_cells->f5[i];
  }

  //  Extract individual parameters from problemDef
  // Extract individual parameters from problemDef
  geometry.set_size((&hi_emlrtRTEI), sp, 1, problemDef->geometry.size(1));
  loop_ub = problemDef->geometry.size(0) * problemDef->geometry.size(1);
  for (i = 0; i < loop_ub; i++) {
    geometry[i] = problemDef->geometry[i];
  }

  cBacks.set_size((&ii_emlrtRTEI), sp, 1, problemDef->contrastBacks.size(1));
  loop_ub = problemDef->contrastBacks.size(0) * problemDef->contrastBacks.size(1);
  for (i = 0; i < loop_ub; i++) {
    cBacks[i] = problemDef->contrastBacks[i];
  }

  cShifts.set_size((&ji_emlrtRTEI), sp, 1, problemDef->contrastShifts.size(1));
  loop_ub = problemDef->contrastShifts.size(0) * problemDef->contrastShifts.size
    (1);
  for (i = 0; i < loop_ub; i++) {
    cShifts[i] = problemDef->contrastShifts[i];
  }

  cScales.set_size((&ki_emlrtRTEI), sp, 1, problemDef->contrastScales.size(1));
  loop_ub = problemDef->contrastScales.size(0) * problemDef->contrastScales.size
    (1);
  for (i = 0; i < loop_ub; i++) {
    cScales[i] = problemDef->contrastScales[i];
  }

  cNbas.set_size((&li_emlrtRTEI), sp, 1, problemDef->contrastNbas.size(1));
  loop_ub = problemDef->contrastNbas.size(0) * problemDef->contrastNbas.size(1);
  for (i = 0; i < loop_ub; i++) {
    cNbas[i] = problemDef->contrastNbas[i];
  }

  cNbss.set_size((&mi_emlrtRTEI), sp, 1, problemDef->contrastNbss.size(1));
  loop_ub = problemDef->contrastNbss.size(0) * problemDef->contrastNbss.size(1);
  for (i = 0; i < loop_ub; i++) {
    cNbss[i] = problemDef->contrastNbss[i];
  }

  cRes.set_size((&ni_emlrtRTEI), sp, 1, problemDef->contrastRes.size(1));
  loop_ub = problemDef->contrastRes.size(0) * problemDef->contrastRes.size(1);
  for (i = 0; i < loop_ub; i++) {
    cRes[i] = problemDef->contrastRes[i];
  }

  backs.set_size((&oi_emlrtRTEI), sp, 1, problemDef->backs.size(1));
  loop_ub = problemDef->backs.size(0) * problemDef->backs.size(1);
  for (i = 0; i < loop_ub; i++) {
    backs[i] = problemDef->backs[i];
  }

  shifts.set_size((&pi_emlrtRTEI), sp, 1, problemDef->shifts.size(1));
  loop_ub = problemDef->shifts.size(0) * problemDef->shifts.size(1);
  for (i = 0; i < loop_ub; i++) {
    shifts[i] = problemDef->shifts[i];
  }

  sf.set_size((&qi_emlrtRTEI), sp, 1, problemDef->sf.size(1));
  loop_ub = problemDef->sf.size(0) * problemDef->sf.size(1);
  for (i = 0; i < loop_ub; i++) {
    sf[i] = problemDef->sf[i];
  }

  nba.set_size((&ri_emlrtRTEI), sp, 1, problemDef->nba.size(1));
  loop_ub = problemDef->nba.size(0) * problemDef->nba.size(1);
  for (i = 0; i < loop_ub; i++) {
    nba[i] = problemDef->nba[i];
  }

  nbs.set_size((&si_emlrtRTEI), sp, 1, problemDef->nbs.size(1));
  loop_ub = problemDef->nbs.size(0) * problemDef->nbs.size(1);
  for (i = 0; i < loop_ub; i++) {
    nbs[i] = problemDef->nbs[i];
  }

  res.set_size((&ti_emlrtRTEI), sp, 1, problemDef->res.size(1));
  loop_ub = problemDef->res.size(0) * problemDef->res.size(1);
  for (i = 0; i < loop_ub; i++) {
    res[i] = problemDef->res[i];
  }

  dataPresent.set_size((&ui_emlrtRTEI), sp, 1, problemDef->dataPresent.size(1));
  loop_ub = problemDef->dataPresent.size(0) * problemDef->dataPresent.size(1);
  for (i = 0; i < loop_ub; i++) {
    dataPresent[i] = problemDef->dataPresent[i];
  }

  nParams = problemDef->params.size(1);
  params.set_size((&vi_emlrtRTEI), sp, 1, problemDef->params.size(1));
  loop_ub = problemDef->params.size(0) * problemDef->params.size(1);
  for (i = 0; i < loop_ub; i++) {
    params[i] = problemDef->params[i];
  }

  resample.set_size((&wi_emlrtRTEI), sp, 1, problemDef->resample.size(1));
  loop_ub = problemDef->resample.size(0) * problemDef->resample.size(1);
  for (i = 0; i < loop_ub; i++) {
    resample[i] = problemDef->resample[i];
  }

  backsType.set_size((&xi_emlrtRTEI), sp, 1, problemDef->contrastBacksType.size
                     (1));
  loop_ub = problemDef->contrastBacksType.size(0) *
    problemDef->contrastBacksType.size(1);
  for (i = 0; i < loop_ub; i++) {
    backsType[i] = problemDef->contrastBacksType[i];
  }

  calcSld = controls->calcSld;

  //  Allocate the memory for the output arrays before the main loop
  //    --- Begin memory allocation ---
  if (!(problemDef->numberOfContrasts >= 0.0)) {
    emlrtNonNegativeCheckR2012b(problemDef->numberOfContrasts, &sb_emlrtDCI, sp);
  }

  i = static_cast<int32_T>(muDoubleScalarFloor(problemDef->numberOfContrasts));
  if (problemDef->numberOfContrasts != i) {
    emlrtIntegerCheckR2012b(problemDef->numberOfContrasts, &rb_emlrtDCI, sp);
  }

  loop_ub = static_cast<int32_T>(problemDef->numberOfContrasts);
  backgs.set_size((&yi_emlrtRTEI), sp, loop_ub);
  if (problemDef->numberOfContrasts != i) {
    emlrtIntegerCheckR2012b(problemDef->numberOfContrasts, &qb_emlrtDCI, sp);
  }

  emlrtForLoopVectorCheckR2012b(1.0, 1.0, problemDef->numberOfContrasts,
    mxDOUBLE_CLASS, static_cast<int32_T>(problemDef->numberOfContrasts),
    &ac_emlrtRTEI, sp);
  reflectivity.set_size((&aj_emlrtRTEI), sp, loop_ub);
  for (b_i = 0; b_i < loop_ub; b_i++) {
    i1 = reflectivity.size(0) - 1;
    if (b_i > i1) {
      emlrtDynamicBoundsCheckR2012b(b_i, 0, i1, &vh_emlrtBCI, sp);
    }

    reflectivity[b_i].f1.set_size((&bj_emlrtRTEI), sp, 2, reflectivity[b_i].
      f1.size(1));
    i1 = reflectivity.size(0) - 1;
    if (b_i > i1) {
      emlrtDynamicBoundsCheckR2012b(b_i, 0, i1, &vh_emlrtBCI, sp);
    }

    reflectivity[b_i].f1.set_size((&bj_emlrtRTEI), sp, reflectivity[b_i].f1.size
      (0), 2);
    i1 = reflectivity.size(0) - 1;
    i2 = reflectivity.size(0) - 1;
    if (b_i > i2) {
      emlrtDynamicBoundsCheckR2012b(b_i, 0, i2, &vh_emlrtBCI, sp);
    }

    i2 = reflectivity.size(0) - 1;
    if (b_i > i2) {
      emlrtDynamicBoundsCheckR2012b(b_i, 0, i2, &vh_emlrtBCI, sp);
    }

    if (b_i > i1) {
      emlrtDynamicBoundsCheckR2012b(b_i, 0, i1, &vh_emlrtBCI, sp);
    }

    reflectivity[b_i].f1[0] = 1.0;
    if (b_i > i1) {
      emlrtDynamicBoundsCheckR2012b(b_i, 0, i1, &vh_emlrtBCI, sp);
    }

    reflectivity[b_i].f1[1] = 1.0;
    if (b_i > i1) {
      emlrtDynamicBoundsCheckR2012b(b_i, 0, i1, &vh_emlrtBCI, sp);
    }

    reflectivity[b_i].f1[2] = 1.0;
    if (b_i > i1) {
      emlrtDynamicBoundsCheckR2012b(b_i, 0, i1, &vh_emlrtBCI, sp);
    }

    reflectivity[b_i].f1[3] = 1.0;
  }

  emlrtForLoopVectorCheckR2012b(1.0, 1.0, problemDef->numberOfContrasts,
    mxDOUBLE_CLASS, static_cast<int32_T>(problemDef->numberOfContrasts),
    &bc_emlrtRTEI, sp);
  Simulation.set_size((&aj_emlrtRTEI), sp, loop_ub);
  for (b_i = 0; b_i < loop_ub; b_i++) {
    i1 = Simulation.size(0) - 1;
    if (b_i > i1) {
      emlrtDynamicBoundsCheckR2012b(b_i, 0, i1, &wh_emlrtBCI, sp);
    }

    Simulation[b_i].f1.set_size((&cj_emlrtRTEI), sp, 2, Simulation[b_i].f1.size
      (1));
    i1 = Simulation.size(0) - 1;
    if (b_i > i1) {
      emlrtDynamicBoundsCheckR2012b(b_i, 0, i1, &wh_emlrtBCI, sp);
    }

    Simulation[b_i].f1.set_size((&cj_emlrtRTEI), sp, Simulation[b_i].f1.size(0),
      2);
    i1 = Simulation.size(0) - 1;
    i2 = Simulation.size(0) - 1;
    if (b_i > i2) {
      emlrtDynamicBoundsCheckR2012b(b_i, 0, i2, &wh_emlrtBCI, sp);
    }

    i2 = Simulation.size(0) - 1;
    if (b_i > i2) {
      emlrtDynamicBoundsCheckR2012b(b_i, 0, i2, &wh_emlrtBCI, sp);
    }

    if (b_i > i1) {
      emlrtDynamicBoundsCheckR2012b(b_i, 0, i1, &wh_emlrtBCI, sp);
    }

    Simulation[b_i].f1[0] = 1.0;
    if (b_i > i1) {
      emlrtDynamicBoundsCheckR2012b(b_i, 0, i1, &wh_emlrtBCI, sp);
    }

    Simulation[b_i].f1[1] = 1.0;
    if (b_i > i1) {
      emlrtDynamicBoundsCheckR2012b(b_i, 0, i1, &wh_emlrtBCI, sp);
    }

    Simulation[b_i].f1[2] = 1.0;
    if (b_i > i1) {
      emlrtDynamicBoundsCheckR2012b(b_i, 0, i1, &wh_emlrtBCI, sp);
    }

    Simulation[b_i].f1[3] = 1.0;
  }

  emlrtForLoopVectorCheckR2012b(1.0, 1.0, problemDef->numberOfContrasts,
    mxDOUBLE_CLASS, static_cast<int32_T>(problemDef->numberOfContrasts),
    &cc_emlrtRTEI, sp);
  allLayers.set_size((&aj_emlrtRTEI), sp, loop_ub);
  for (b_i = 0; b_i < loop_ub; b_i++) {
    i1 = allLayers.size(0) - 1;
    i2 = allLayers.size(0) - 1;
    if (b_i > i2) {
      emlrtDynamicBoundsCheckR2012b(b_i, 0, i2, &xh_emlrtBCI, sp);
    }

    allLayers[b_i].f1.set_size((&dj_emlrtRTEI), sp, 2);
    if (b_i > i1) {
      emlrtDynamicBoundsCheckR2012b(b_i, 0, i1, &xh_emlrtBCI, sp);
    }

    allLayers[b_i].f1[0] = 1.0;
    if (b_i > i1) {
      emlrtDynamicBoundsCheckR2012b(b_i, 0, i1, &xh_emlrtBCI, sp);
    }

    allLayers[b_i].f1[1] = 1.0;
  }

  //    --- End memory allocation ---
  //  First we need to allocate the absolute values of the input
  //  parameters to all the layers in the layers list. This only needs
  //  to be done once, and so is done outside the contrasts loop
  st.site = &mg_emlrtRSI;
  allocateParamsToLayers(&st, problemDef->params, problemDef_cells->f6,
    outParameterisedLayers);
  if (problemDef->numberOfContrasts != i) {
    emlrtIntegerCheckR2012b(problemDef->numberOfContrasts, &pb_emlrtDCI, sp);
  }

  emlrtForLoopVectorCheckR2012b(1.0, 1.0, problemDef->numberOfContrasts,
    mxDOUBLE_CLASS, static_cast<int32_T>(problemDef->numberOfContrasts),
    &dc_emlrtRTEI, sp);
  outSsubs.set_size((&ej_emlrtRTEI), sp, loop_ub);
  sldProfiles.set_size((&aj_emlrtRTEI), sp, loop_ub);
  reflectivity.set_size((&aj_emlrtRTEI), sp, loop_ub);
  Simulation.set_size((&aj_emlrtRTEI), sp, loop_ub);
  shifted_data.set_size((&aj_emlrtRTEI), sp, loop_ub);
  layerSlds.set_size((&aj_emlrtRTEI), sp, loop_ub);
  chis.set_size((&ej_emlrtRTEI), sp, loop_ub);
  qshifts.set_size((&ej_emlrtRTEI), sp, loop_ub);
  sfs.set_size((&ej_emlrtRTEI), sp, loop_ub);
  nbas.set_size((&ej_emlrtRTEI), sp, loop_ub);
  nbss.set_size((&ej_emlrtRTEI), sp, loop_ub);
  resols.set_size((&ej_emlrtRTEI), sp, loop_ub);
  allRoughs.set_size((&ej_emlrtRTEI), sp, loop_ub);
  loop_ub--;
  emlrtEnterParallelRegion(sp, omp_in_parallel());

#pragma omp parallel \
 num_threads(emlrtAllocRegionTLSs(sp->tls, omp_in_parallel(), omp_get_max_threads(), omp_get_num_procs())) \
 private(thisContrastLayers,sldProfile,reflect,Simul,shifted_dat,thisData,thisSsubs,thisChiSquared,thisResol,thisNbs,thisNba,thisSf,thisQshift,thisBackground,c_st,i3,c_i,b_loop_ub,i4) \
 firstprivate(b_st,emlrtHadParallelError)

  {
    try {
      b_st.prev = sp;
      b_st.tls = emlrtAllocTLS((emlrtStack *)sp, omp_get_thread_num());
      b_st.site = NULL;
      c_st.prev = &b_st;
      c_st.tls = b_st.tls;
    } catch (...) {
      emlrtHadParallelError = true;
    }

#pragma omp for nowait

    for (c_i = 0; c_i <= loop_ub; c_i++) {
      if (emlrtHadParallelError)
        continue;
      try {
        //  Parallel Loop over all the contrasts
        //  Extract the relevant parameter values for this contrast
        //  from the input arrays.
        //  First need to decide which values of the backrounds, scalefactors
        //  data shifts and bulk contrasts are associated with this contrast
        i3 = static_cast<int32_T>(c_i + 1U);
        if ((i3 < 1) || (i3 > cBacks.size(1))) {
          emlrtDynamicBoundsCheckR2012b(i3, 1, cBacks.size(1), &yh_emlrtBCI,
            &b_st);
        }

        if (i3 > cShifts.size(1)) {
          emlrtDynamicBoundsCheckR2012b(i3, 1, cShifts.size(1), &ai_emlrtBCI,
            &b_st);
        }

        if (i3 > cScales.size(1)) {
          emlrtDynamicBoundsCheckR2012b(i3, 1, cScales.size(1), &bi_emlrtBCI,
            &b_st);
        }

        if (i3 > cNbas.size(1)) {
          emlrtDynamicBoundsCheckR2012b(i3, 1, cNbas.size(1), &ci_emlrtBCI,
            &b_st);
        }

        if (i3 > cNbss.size(1)) {
          emlrtDynamicBoundsCheckR2012b(i3, 1, cNbss.size(1), &di_emlrtBCI,
            &b_st);
        }

        if (i3 > cRes.size(1)) {
          emlrtDynamicBoundsCheckR2012b(i3, 1, cRes.size(1), &ei_emlrtBCI, &b_st);
        }

        c_st.site = &ng_emlrtRSI;
        backSort(&c_st, cBacks[i3 - 1], cShifts[i3 - 1], cScales[i3 - 1],
                 cNbas[i3 - 1], cNbss[i3 - 1], cRes[i3 - 1], backs, shifts, sf,
                 nba, nbs, res, &thisBackground, &thisQshift, &thisSf, &thisNba,
                 &thisNbs, &thisResol);

        //  Also need to determine which layers from the overall layers list
        //  are required for this contrast, and put them in the correct order
        //  according to geometry
        i3 = contrastLayers.size(1) - 1;
        if (c_i > i3) {
          emlrtDynamicBoundsCheckR2012b(c_i, 0, i3, &fi_emlrtBCI, &b_st);
        }

        c_st.site = &og_emlrtRSI;
        allocateLayersForContrast(&c_st, contrastLayers[c_i].f1,
          outParameterisedLayers, thisContrastLayers);

        //  For the other parameters, we extract the correct ones from the input 
        //  arrays
        if (1 > params.size(1)) {
          emlrtDynamicBoundsCheckR2012b(1, 1, params.size(1), &nh_emlrtBCI,
            &b_st);
        }

        //  Substrate roughness is always first parameter for standard layers
        i3 = repeatLayers.size(1) - 1;
        if (c_i > i3) {
          emlrtDynamicBoundsCheckR2012b(c_i, 0, i3, &rh_emlrtBCI, &b_st);
        }

        i3 = c_i + 1;
        if ((i3 < 1) || (i3 > resample.size(1))) {
          emlrtDynamicBoundsCheckR2012b(i3, 1, resample.size(1), &oh_emlrtBCI,
            &b_st);
        }

        i3 = allData.size(1) - 1;
        if (c_i > i3) {
          emlrtDynamicBoundsCheckR2012b(c_i, 0, i3, &sh_emlrtBCI, &b_st);
        }

        i3 = c_i + 1;
        if ((i3 < 1) || (i3 > dataPresent.size(1))) {
          emlrtDynamicBoundsCheckR2012b(i3, 1, dataPresent.size(1), &ph_emlrtBCI,
            &b_st);
        }

        i3 = dataLimits.size(1) - 1;
        if (c_i > i3) {
          emlrtDynamicBoundsCheckR2012b(c_i, 0, i3, &th_emlrtBCI, &b_st);
        }

        i3 = simLimits.size(1) - 1;
        if (c_i > i3) {
          emlrtDynamicBoundsCheckR2012b(c_i, 0, i3, &uh_emlrtBCI, &b_st);
        }

        i3 = c_i + 1;
        if ((i3 < 1) || (i3 > backsType.size(1))) {
          emlrtDynamicBoundsCheckR2012b(i3, 1, backsType.size(1), &qh_emlrtBCI,
            &b_st);
        }

        //  Now call the core standardTF_stanlay reflectivity calculation
        //  In this case we are single cored  over points.
        //  Call the calculation
        c_st.site = &pg_emlrtRSI;
        standardTF_layers_core(&c_st, thisContrastLayers, params[0], geometry,
          thisNba, thisNbs, resample[c_i], calcSld, thisSf, thisQshift,
          dataPresent[c_i], allData[c_i].f1, dataLimits[c_i].f1, simLimits[c_i].
          f1, repeatLayers[c_i].f1, thisBackground, thisResol, backsType[c_i],
          static_cast<real_T>(nParams), sldProfile, reflect, Simul, shifted_dat,
          thisData, &thisChiSquared, &thisSsubs);

        //  Store returned values for this contrast in the output arrays.
        //  As well as the calculated profiles, we also store a record of
        //  the other values (background, scalefactors etc) for each contrast
        //  for future use.
        i3 = c_i + 1;
        if ((i3 < 1) || (i3 > outSsubs.size(0))) {
          emlrtDynamicBoundsCheckR2012b(i3, 1, outSsubs.size(0), &gi_emlrtBCI,
            &b_st);
        }

        outSsubs[i3 - 1] = thisSsubs;
        b_loop_ub = sldProfile.size(0) * sldProfile.size(1);
        i3 = sldProfiles.size(0) - 1;
        if (c_i > i3) {
          emlrtDynamicBoundsCheckR2012b(c_i, 0, i3, &hi_emlrtBCI, &b_st);
        }

        sldProfiles[c_i].f1.set_size((&fj_emlrtRTEI), (&b_st), sldProfile.size(0),
          sldProfiles[c_i].f1.size(1));
        i3 = sldProfiles.size(0) - 1;
        if (c_i > i3) {
          emlrtDynamicBoundsCheckR2012b(c_i, 0, i3, &hi_emlrtBCI, &b_st);
        }

        sldProfiles[c_i].f1.set_size((&fj_emlrtRTEI), (&b_st), sldProfiles[c_i].
          f1.size(0), sldProfile.size(1));
        i3 = sldProfiles.size(0) - 1;
        i4 = sldProfiles.size(0) - 1;
        if (c_i > i4) {
          emlrtDynamicBoundsCheckR2012b(c_i, 0, i4, &hi_emlrtBCI, &b_st);
        }

        i4 = sldProfiles.size(0) - 1;
        if (c_i > i4) {
          emlrtDynamicBoundsCheckR2012b(c_i, 0, i4, &hi_emlrtBCI, &b_st);
        }

        for (i4 = 0; i4 < b_loop_ub; i4++) {
          if (c_i > i3) {
            emlrtDynamicBoundsCheckR2012b(c_i, 0, i3, &hi_emlrtBCI, &b_st);
          }

          sldProfiles[c_i].f1[i4] = sldProfile[i4];
        }

        b_loop_ub = reflect.size(0) * reflect.size(1);
        i3 = reflectivity.size(0) - 1;
        if (c_i > i3) {
          emlrtDynamicBoundsCheckR2012b(c_i, 0, i3, &ii_emlrtBCI, &b_st);
        }

        reflectivity[c_i].f1.set_size((&gj_emlrtRTEI), (&b_st), reflect.size(0),
          reflectivity[c_i].f1.size(1));
        i3 = reflectivity.size(0) - 1;
        if (c_i > i3) {
          emlrtDynamicBoundsCheckR2012b(c_i, 0, i3, &ii_emlrtBCI, &b_st);
        }

        reflectivity[c_i].f1.set_size((&gj_emlrtRTEI), (&b_st), reflectivity[c_i]
          .f1.size(0), 2);
        i3 = reflectivity.size(0) - 1;
        i4 = reflectivity.size(0) - 1;
        if (c_i > i4) {
          emlrtDynamicBoundsCheckR2012b(c_i, 0, i4, &ii_emlrtBCI, &b_st);
        }

        i4 = reflectivity.size(0) - 1;
        if (c_i > i4) {
          emlrtDynamicBoundsCheckR2012b(c_i, 0, i4, &ii_emlrtBCI, &b_st);
        }

        for (i4 = 0; i4 < b_loop_ub; i4++) {
          if (c_i > i3) {
            emlrtDynamicBoundsCheckR2012b(c_i, 0, i3, &ii_emlrtBCI, &b_st);
          }

          reflectivity[c_i].f1[i4] = reflect[i4];
        }

        b_loop_ub = Simul.size(0) * Simul.size(1);
        i3 = Simulation.size(0) - 1;
        if (c_i > i3) {
          emlrtDynamicBoundsCheckR2012b(c_i, 0, i3, &ji_emlrtBCI, &b_st);
        }

        Simulation[c_i].f1.set_size((&hj_emlrtRTEI), (&b_st), Simul.size(0),
          Simulation[c_i].f1.size(1));
        i3 = Simulation.size(0) - 1;
        if (c_i > i3) {
          emlrtDynamicBoundsCheckR2012b(c_i, 0, i3, &ji_emlrtBCI, &b_st);
        }

        Simulation[c_i].f1.set_size((&hj_emlrtRTEI), (&b_st), Simulation[c_i].
          f1.size(0), 2);
        i3 = Simulation.size(0) - 1;
        i4 = Simulation.size(0) - 1;
        if (c_i > i4) {
          emlrtDynamicBoundsCheckR2012b(c_i, 0, i4, &ji_emlrtBCI, &b_st);
        }

        i4 = Simulation.size(0) - 1;
        if (c_i > i4) {
          emlrtDynamicBoundsCheckR2012b(c_i, 0, i4, &ji_emlrtBCI, &b_st);
        }

        for (i4 = 0; i4 < b_loop_ub; i4++) {
          if (c_i > i3) {
            emlrtDynamicBoundsCheckR2012b(c_i, 0, i3, &ji_emlrtBCI, &b_st);
          }

          Simulation[c_i].f1[i4] = Simul[i4];
        }

        b_loop_ub = shifted_dat.size(0) * shifted_dat.size(1);
        i3 = shifted_data.size(0) - 1;
        if (c_i > i3) {
          emlrtDynamicBoundsCheckR2012b(c_i, 0, i3, &ki_emlrtBCI, &b_st);
        }

        shifted_data[c_i].f1.set_size((&ij_emlrtRTEI), (&b_st), shifted_dat.size
          (0), shifted_data[c_i].f1.size(1));
        i3 = shifted_data.size(0) - 1;
        if (c_i > i3) {
          emlrtDynamicBoundsCheckR2012b(c_i, 0, i3, &ki_emlrtBCI, &b_st);
        }

        shifted_data[c_i].f1.set_size((&ij_emlrtRTEI), (&b_st), shifted_data[c_i]
          .f1.size(0), 3);
        i3 = shifted_data.size(0) - 1;
        i4 = shifted_data.size(0) - 1;
        if (c_i > i4) {
          emlrtDynamicBoundsCheckR2012b(c_i, 0, i4, &ki_emlrtBCI, &b_st);
        }

        i4 = shifted_data.size(0) - 1;
        if (c_i > i4) {
          emlrtDynamicBoundsCheckR2012b(c_i, 0, i4, &ki_emlrtBCI, &b_st);
        }

        for (i4 = 0; i4 < b_loop_ub; i4++) {
          if (c_i > i3) {
            emlrtDynamicBoundsCheckR2012b(c_i, 0, i3, &ki_emlrtBCI, &b_st);
          }

          shifted_data[c_i].f1[i4] = shifted_dat[i4];
        }

        b_loop_ub = thisData.size(0) * thisData.size(1);
        i3 = layerSlds.size(0) - 1;
        if (c_i > i3) {
          emlrtDynamicBoundsCheckR2012b(c_i, 0, i3, &li_emlrtBCI, &b_st);
        }

        layerSlds[c_i].f1.set_size((&jj_emlrtRTEI), (&b_st), thisData.size(0),
          layerSlds[c_i].f1.size(1));
        i3 = layerSlds.size(0) - 1;
        if (c_i > i3) {
          emlrtDynamicBoundsCheckR2012b(c_i, 0, i3, &li_emlrtBCI, &b_st);
        }

        layerSlds[c_i].f1.set_size((&jj_emlrtRTEI), (&b_st), layerSlds[c_i].
          f1.size(0), 3);
        i3 = layerSlds.size(0) - 1;
        i4 = layerSlds.size(0) - 1;
        if (c_i > i4) {
          emlrtDynamicBoundsCheckR2012b(c_i, 0, i4, &li_emlrtBCI, &b_st);
        }

        i4 = layerSlds.size(0) - 1;
        if (c_i > i4) {
          emlrtDynamicBoundsCheckR2012b(c_i, 0, i4, &li_emlrtBCI, &b_st);
        }

        for (i4 = 0; i4 < b_loop_ub; i4++) {
          if (c_i > i3) {
            emlrtDynamicBoundsCheckR2012b(c_i, 0, i3, &li_emlrtBCI, &b_st);
          }

          layerSlds[c_i].f1[i4] = thisData[i4];
        }

        i3 = c_i + 1;
        if ((i3 < 1) || (i3 > chis.size(0))) {
          emlrtDynamicBoundsCheckR2012b(i3, 1, chis.size(0), &mi_emlrtBCI, &b_st);
        }

        chis[i3 - 1] = thisChiSquared;
        i3 = c_i + 1;
        if ((i3 < 1) || (i3 > backgs.size(0))) {
          emlrtDynamicBoundsCheckR2012b(i3, 1, backgs.size(0), &ni_emlrtBCI,
            &b_st);
        }

        backgs[i3 - 1] = thisBackground;
        i3 = c_i + 1;
        if ((i3 < 1) || (i3 > qshifts.size(0))) {
          emlrtDynamicBoundsCheckR2012b(i3, 1, qshifts.size(0), &oi_emlrtBCI,
            &b_st);
        }

        qshifts[i3 - 1] = thisQshift;
        i3 = c_i + 1;
        if ((i3 < 1) || (i3 > sfs.size(0))) {
          emlrtDynamicBoundsCheckR2012b(i3, 1, sfs.size(0), &pi_emlrtBCI, &b_st);
        }

        sfs[i3 - 1] = thisSf;
        i3 = c_i + 1;
        if ((i3 < 1) || (i3 > nbas.size(0))) {
          emlrtDynamicBoundsCheckR2012b(i3, 1, nbas.size(0), &qi_emlrtBCI, &b_st);
        }

        nbas[i3 - 1] = thisNba;
        i3 = c_i + 1;
        if ((i3 < 1) || (i3 > nbss.size(0))) {
          emlrtDynamicBoundsCheckR2012b(i3, 1, nbss.size(0), &ri_emlrtBCI, &b_st);
        }

        nbss[i3 - 1] = thisNbs;
        i3 = c_i + 1;
        if ((i3 < 1) || (i3 > resols.size(0))) {
          emlrtDynamicBoundsCheckR2012b(i3, 1, resols.size(0), &si_emlrtBCI,
            &b_st);
        }

        resols[i3 - 1] = thisResol;
        i3 = c_i + 1;
        if ((i3 < 1) || (i3 > allRoughs.size(0))) {
          emlrtDynamicBoundsCheckR2012b(i3, 1, allRoughs.size(0), &ti_emlrtBCI,
            &b_st);
        }

        allRoughs[i3 - 1] = params[0];
      } catch (...) {
        emlrtHadParallelError = true;
      }
    }

    if (!emlrtHadParallelError) {
      emlrtHeapReferenceStackLeaveScope(&b_st, 6);
    }
  }

  emlrtExitParallelRegion(sp, omp_in_parallel());
  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

// End of code generation (standardTF_stanlay_paraContrasts.cpp)