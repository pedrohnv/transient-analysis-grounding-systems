#include "interface_matlab.h"
#include "mex.h"
#include "electrode.h"
#include <complex.h>
#include <string.h>

int
cast_electrode (const mxArray *matlab_elect, mwIndex index, Electrode* electrode)
{
    int number_of_fields, field_index;
    const char *field_name;
    const mxArray *field_array_ptr;
    mxDouble *p;
    #if MX_HAS_INTERLEAVED_COMPLEX
        mxComplexDouble *pc;
    #else
        double *pr, *pi;
    #endif
    number_of_fields = mxGetNumberOfFields(matlab_elect);
    if (number_of_fields != 6) {
        mexErrMsgTxt("Electrode structure should have exactly 6 fields\n");
    }
    for (field_index = 0; field_index < number_of_fields; field_index++) {
        field_name = mxGetFieldNameByNumber(matlab_elect, field_index);
        field_array_ptr = mxGetFieldByNumber(matlab_elect, index, field_index);
        if (field_array_ptr == NULL) {
            mexPrintf("\tEmpty Field\n");
            return 1;
        }
        if (strcmp(field_name, "start_point") == 0) {
            p = (double *) mxGetPr(field_array_ptr);
            electrode->start_point[0] = *p++;
            electrode->start_point[1] = *p++;
            electrode->start_point[2] = *p++;
        } else if (strcmp(field_name, "end_point") == 0) {
            p = (double *) mxGetPr(field_array_ptr);
            electrode->end_point[0] = *p++;
            electrode->end_point[1] = *p++;
            electrode->end_point[2] = *p++;
        } else if (strcmp(field_name, "middle_point") == 0) {
            p = (double *) mxGetPr(field_array_ptr);
            electrode->middle_point[0] = *p++;
            electrode->middle_point[1] = *p++;
            electrode->middle_point[2] = *p++;
        } else if (strcmp(field_name, "length") == 0) {
            electrode->length = (double) mxGetScalar(field_array_ptr);
        } else if (strcmp(field_name, "radius") == 0) {
            electrode->radius = (double) mxGetScalar(field_array_ptr);
        } else if (strcmp(field_name, "zi") == 0) {
            #if MX_HAS_INTERLEAVED_COMPLEX
                mxComplexDouble *pc;
                mxDouble *p;
                if (mxIsComplex(field_array_ptr)) {
                    pc = mxGetComplexDoubles(field_array_ptr);
                    electrode->zi = (*pc).real + (*pc).imag*I;
                } else {
                    p = mxGetDoubles(field_array_ptr);
                    electrode->zi = (*p) + I*0.0;
                }
            #else
                pr = (double *) mxGetPr(field_array_ptr);
                electrode->zi = (*pr);
                if (mxIsComplex(field_array_ptr)) {
                    pi = (double *) mxGetPi(field_array_ptr);
                    electrode->zi += (*pi)*I;
                }
            #endif
        } else {
            mexPrintf("Unrecognized Electrode structure field: %s\n", field_name);
            mexErrMsgTxt("");
        }
    }
    return 0;
}

_Complex double
get_complex(const mxArray *array_ptr)
{
    _Complex double value;
    #if MX_HAS_INTERLEAVED_COMPLEX
        mxComplexDouble *pc;
        mxDouble *p;
        if (mxIsComplex(array_ptr)) {
            pc = mxGetComplexDoubles(array_ptr);
            value = (*pc).real + I*(*pc).imag;
        }
        else {
            p = mxGetDoubles(array_ptr);
            value = (*p);
        }
    #else
        double *pr, *pi;
        pr = (double *) mxGetPr(array_ptr);
        if (mxIsComplex(array_ptr)) {
            pi = (double *) mxGetPi(array_ptr);
            value = (*pr) + (*pi)*I;
        } else {
            value = (*pr) + 0.0*I;
        }
    #endif
    return value;
}
