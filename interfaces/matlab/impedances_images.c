#include "mex.h"
#include "interface_matlab.h"
#include "electrode.h"
#include <complex.h>
#include <string.h>

/**
Add the image effects to the impedance matrices of the electrode system.
TODO describe args
*/
void
mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 15) {
        mexErrMsgTxt("Exactly 15 arguments are expected.\n");
    }
    if (mxGetClassID(prhs[0]) != mxSTRUCT_CLASS) {
        mexErrMsgTxt("1st argument must be an Electrode struct.\n");
    }
    if (mxGetClassID(prhs[1]) != mxSTRUCT_CLASS) {
        mexErrMsgTxt("2nd argument must be an Electrode struct.\n");
    }
    if (mxGetClassID(prhs[2]) != mxDOUBLE_CLASS) {
        mexErrMsgTxt("3rd argument must be a complex matrix.\n");
        // TODO check size
    }
    if (mxGetClassID(prhs[3]) != mxDOUBLE_CLASS) {
        mexErrMsgTxt("4th argument must be a complex matrix.\n");
        // TODO check size
    }
    if (mxGetClassID(prhs[4]) != mxDOUBLE_CLASS) {
        mexErrMsgTxt("5th argument must be a complex number.\n");
    }
    if (mxGetClassID(prhs[5]) != mxDOUBLE_CLASS) {
        mexErrMsgTxt("6th argument must be a complex number.\n");
    }
    if (mxGetClassID(prhs[6]) != mxDOUBLE_CLASS
        || mxIsComplex(prhs[6])) {
        mexErrMsgTxt("7th argument must be a real number.\n");
    }
    if (mxGetClassID(prhs[7]) != mxDOUBLE_CLASS) {
        mexErrMsgTxt("8th argument must be a complex number.\n");
    }
    if (mxGetClassID(prhs[8]) != mxDOUBLE_CLASS) {
        mexErrMsgTxt("9th argument must be a complex number.\n");
    }
    if (mxGetClassID(prhs[9]) != mxDOUBLE_CLASS) {
        mexErrMsgTxt("10th argument must be a complex number.\n");
    }
    if (mxGetClassID(prhs[10]) != mxDOUBLE_CLASS
        || mxIsComplex(prhs[10])) {
        mexErrMsgTxt("11th argument must be a real integer.\n");
    }
    if (mxGetClassID(prhs[11]) != mxDOUBLE_CLASS
        || mxIsComplex(prhs[11])) {
        mexErrMsgTxt("12th argument must be a real number.\n");
    }
    if (mxGetClassID(prhs[12]) != mxDOUBLE_CLASS
        || mxIsComplex(prhs[12])) {
        mexErrMsgTxt("13th argument must be a real number.\n");
    }
    // FIXME get Enumeration value
    if (mxGetClassID(prhs[13]) != mxINT16_CLASS
        || mxIsComplex(prhs[13])) {
        mexErrMsgTxt("14th argument must be an INT16.\n");
    }
    if (mxGetClassID(prhs[14]) != mxINT16_CLASS
        || mxIsComplex(prhs[14])) {
        mexErrMsgTxt("15th argument must be an INT16.\n");
    }
    // ================
    size_t num_electrodes = mxGetNumberOfElements(prhs[0]);
    Electrode *electrodes = malloc(sizeof(Electrode)*num_electrodes);
    Electrode *images = malloc(sizeof(Electrode)*num_electrodes);
    for (size_t i = 0; i < num_electrodes; i++) {
        cast_electrode(prhs[0], i, &(electrodes[i]));
        cast_electrode(prhs[1], i, &(images[i]));
    }
    /* Too bad MATLAB does not accept in-place editing of variables.
    Doing so may set the computer ablaze... Though it seems it can be fooled
    @see https://undocumentedmatlab.com/blog/matlab-mex-in-place-editing */
    size_t ne2 = num_electrodes*num_electrodes;
    _Complex double *zl = malloc(ne2*sizeof(_Complex double));
    _Complex double *zt = malloc(ne2*sizeof(_Complex double));
    #if MX_HAS_INTERLEAVED_COMPLEX
        mxComplexDouble *zldata = mxGetComplexDoubles(prhs[2]);
        mxComplexDouble *ztdata = mxGetComplexDoubles(prhs[3]);
        for (size_t k = 0; k < num_electrodes; k++) {
            for (size_t i = 0; i < num_electrodes; i++) {
                zl[k*num_electrodes + i] = zldata[k*num_electrodes + i].real
                                         + I*zldata[k*num_electrodes + i].imag;
                zt[k*num_electrodes + i] = ztdata[k*num_electrodes + i].real
                                         + I*ztdata[k*num_electrodes + i].imag;
            }
        }
    #else
        double *zlre, *zlim, *ztre, *ztim;
        zlre = (double *) mxGetPr(prhs[2]);
        zlim = (double *) mxGetPi(prhs[2]);
        ztre = (double *) mxGetPr(prhs[3]);
        ztim = (double *) mxGetPi(prhs[3]);
        for (size_t k = 0; k < num_electrodes; k++) {
            for (size_t i = 0; i < num_electrodes; i++) {
                zl[k*num_electrodes + i] = zlre[k*num_electrodes + i]
                                         + I*zlim[k*num_electrodes + i];
                zt[k*num_electrodes + i] = ztre[k*num_electrodes + i]
                                         + I*ztim[k*num_electrodes + i];
            }
        }
    #endif
    _Complex double gamma = get_complex(prhs[4]);
    _Complex double s = get_complex(prhs[5]);
    double mur = (double) mxGetScalar(prhs[6]);
    _Complex double kappa = get_complex(prhs[7]);
    _Complex double ref_l = get_complex(prhs[8]);
    _Complex double ref_t = get_complex(prhs[9]);
    size_t max_eval = (size_t) mxGetScalar(prhs[10]);
    double req_abs_error = (double) mxGetScalar(prhs[11]);
    double req_rel_error = (double) mxGetScalar(prhs[12]);
    int error_norm = (int) mxGetScalar(prhs[13]);
    int integration_type = (int) mxGetScalar(prhs[14]);
    impedances_images(electrodes, images, num_electrodes, zl, zt, gamma, s, mur,
                      kappa, ref_l, ref_t, max_eval, req_abs_error, req_rel_error,
                      error_norm, integration_type);
    plhs[0] = mxCreateDoubleMatrix((mwSize)num_electrodes, (mwSize)num_electrodes, mxCOMPLEX);
    plhs[1] = mxCreateDoubleMatrix((mwSize)num_electrodes, (mwSize)num_electrodes, mxCOMPLEX);
    #if MX_HAS_INTERLEAVED_COMPLEX
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
