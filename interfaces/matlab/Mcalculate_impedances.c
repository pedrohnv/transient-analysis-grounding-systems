#include "mex.h"
#include "interface_matlab.h"
#include "electrode.h"
#include <complex.h>
#include <string.h>

/**
Calculate the impedance matrices of the electrode system.
@see calculate_impedances.m
*/
void
mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 10) {
        mexErrMsgTxt("Exactly 10 arguments are expected.\n");
    }
    if (mxGetClassID(prhs[0]) != mxSTRUCT_CLASS) {
        mexErrMsgTxt("1st argument must be an Electrode struct.\n");
    }
    if (mxGetClassID(prhs[1]) != mxDOUBLE_CLASS) {
        mexErrMsgTxt("2nd argument must be a complex number.\n");
    }
    if (mxGetClassID(prhs[2]) != mxDOUBLE_CLASS) {
        mexErrMsgTxt("3rd argument must be a complex number.\n");
    }
    if (mxGetClassID(prhs[3]) != mxDOUBLE_CLASS
        || mxIsComplex(prhs[3])) {
        mexErrMsgTxt("4th argument must be a real number.\n");
    }
    if (mxGetClassID(prhs[4]) != mxDOUBLE_CLASS) {
        mexErrMsgTxt("5th argument must be a complex number.\n");
    }
    if (mxGetClassID(prhs[5]) != mxDOUBLE_CLASS
        || mxIsComplex(prhs[5])) {
        mexErrMsgTxt("6th argument must be a real integer.\n");
    }
    if (mxGetClassID(prhs[6]) != mxDOUBLE_CLASS
        || mxIsComplex(prhs[6])) {
        mexErrMsgTxt("7th argument must be a real number.\n");
    }
    if (mxGetClassID(prhs[7]) != mxDOUBLE_CLASS
        || mxIsComplex(prhs[7])) {
        mexErrMsgTxt("8th argument must be a real number.\n");
    }
    // FIXME get Enumeration value
    if (mxGetClassID(prhs[8]) != mxINT16_CLASS
        || mxIsComplex(prhs[8])) {
        mexErrMsgTxt("9th argument must be an INT16.\n");
    }
    if (mxGetClassID(prhs[9]) != mxINT16_CLASS
        || mxIsComplex(prhs[9])) {
        mexErrMsgTxt("10th argument must be an INT16.\n");
    }
    // ================
    size_t num_electrodes = mxGetNumberOfElements(prhs[0]);
    Electrode *electrodes = malloc(sizeof(Electrode)*num_electrodes);
    for (size_t i = 0; i < num_electrodes; i++) {
        cast_electrode(prhs[0], i, &(electrodes[i]));
    }
    _Complex double gamma = get_complex(prhs[1]);
    _Complex double s = get_complex(prhs[2]);
    double mur = (double) mxGetScalar(prhs[3]);
    _Complex double kappa = get_complex(prhs[4]);
    size_t max_eval = (size_t) mxGetScalar(prhs[5]);
    double req_abs_error = (double) mxGetScalar(prhs[6]);
    double req_rel_error = (double) mxGetScalar(prhs[7]);
    // FIXME get Enumeration value
    int error_norm = (int) mxGetScalar(prhs[8]);
    int integration_type = (int) mxGetScalar(prhs[9]);
    size_t ne2 = num_electrodes*num_electrodes;
    _Complex double *zl = malloc(ne2*sizeof(_Complex double));
    _Complex double *zt = malloc(ne2*sizeof(_Complex double));
    calculate_impedances(electrodes, num_electrodes, zl, zt, gamma, s, mur,
                         kappa, max_eval, req_abs_error, req_rel_error,
                         error_norm, integration_type);
    plhs[0] = mxCreateDoubleMatrix((mwSize)num_electrodes, (mwSize)num_electrodes, mxCOMPLEX);
    plhs[1] = mxCreateDoubleMatrix((mwSize)num_electrodes, (mwSize)num_electrodes, mxCOMPLEX);
    #if MX_HAS_INTERLEAVED_COMPLEX
        mxComplexDouble *zldata, *ztdata;// FIXME unknown type
        zldata = mxGetComplexDoubles(plhs[0]);
        ztdata = mxGetComplexDoubles(plhs[1]);
        for (size_t k = 0; k < num_electrodes; k++) {
            for (size_t i = 0; i < num_electrodes; i++) {
                zldata[k*num_electrodes + i].real = creal(zl[k*num_electrodes + i]);
                zldata[k*num_electrodes + i].imag = cimag(zl[k*num_electrodes + i]);
                ztdata[k*num_electrodes + i].real = creal(zt[k*num_electrodes + i]);
                ztdata[k*num_electrodes + i].imag = cimag(zt[k*num_electrodes + i]);
            }
        }
    #else
        double *zlre, *zlim, *ztre, *ztim;
        zlre = (double *) mxGetPr(plhs[0]);
        zlim = (double *) mxGetPi(plhs[0]);
        ztre = (double *) mxGetPr(plhs[1]);
        ztim = (double *) mxGetPi(plhs[1]);
        for (size_t k = 0; k < num_electrodes; k++) {
            for (size_t i = 0; i < num_electrodes; i++) {
                zlre[k*num_electrodes + i] = creal(zl[k*num_electrodes + i]);
                zlim[k*num_electrodes + i] = cimag(zl[k*num_electrodes + i]);
                ztre[k*num_electrodes + i] = creal(zt[k*num_electrodes + i]);
                ztim[k*num_electrodes + i] = cimag(zt[k*num_electrodes + i]);
            }
        }
    #endif
}
