//
//  Academic License - for use in teaching, academic research, and meeting
//  course requirements at degree granting institutions only.  Not for
//  government, commercial, or other organizational use.
//
//  c_butterworth.h
//
//  Code generation for function 'c_butterworth'
//


#ifndef C_BUTTERWORTH_H
#define C_BUTTERWORTH_H

// Include files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
#include "c_butterworth_types.h"

// Function Declarations
extern void c_butterworth(const coder::array<double, 1U> &x, double fc, double
  fs, coder::array<double, 1U> &y);
extern void c_butterworth_initialize();
extern void c_butterworth_terminate();

#endif

// End of code generation (c_butterworth.h)
