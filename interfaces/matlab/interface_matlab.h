/** Funtions to interface the C code to MATLAB.
Uses Matlab's "C Matrix API" to run on Matlab R2017b and
earlier.
*/
#ifndef INTERFACE_MATLAB_H_
#define INTERFACE_MATLAB_H_

#include "mex.h"
#include "electrode.h"
//#include <complex.h>

/** cast_electrode
Cast the matlab Electrode structure to the C struct.
@param matlab_elect pointer to Matlab structure that holds the electrode structure
@param electrode pointer that will be used as the C struct
*/
int
cast_electrode (const mxArray *matlab_elect, mwIndex index, Electrode *electrode);

/** get_complex
Check if argument is complex or real, do the appropriate call and return it
as a complex.
*/
_Complex double
get_complex(const mxArray *array_ptr);

#endif /* INTERFACE_MATLAB_H_ */
