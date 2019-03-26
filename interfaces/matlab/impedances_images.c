#include "interface_matlab.h"
#include "electrode.h"
#include "mex.h"
//mex impedances_images.c InterfaceMatlab.c ..\\..\\cubature\\hcubature.c -I. -I..\\..\\cubature

/**
Add the image effects to the impedance matrices of the electrode system.
TODO describe args
*/
void
mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //TODO check inputs and print useful error messages
    //mxGetNumberOfElements = N_array_elements + N_struct_fields
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
    _Complex double *zl = malloc(sizeof(_Complex double)*ne2);
    _Complex double *zt = malloc(sizeof(_Complex double)*ne2);
    #if MX_HAS_INTERLEAVED_COMPLEX
        mxComplexDoubles *zldata, *ztdata;
        zldata = mxGetComplexDoubles(prhs[0]);
        ztdata = mxGetComplexDoubles(prhs[1]);
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
    _Complex double gamma; //prhs[4]
    _Complex double s; //prhs[5]
    double mur = (double) mxGetScalar(prhs[6]);
    _Complex double kappa; //prhs[7]
    _Complex double ref_l; //prhs[8]
    _Complex double ref_t; //prhs[9]
    size_t max_eval = (size_t) mxGetScalar(prhs[10]);
    double req_abs_error = (double) mxGetScalar(prhs[11]);
    double req_rel_error = (double) mxGetScalar(prhs[12]);
    int error_norm = (int) mxGetScalar(prhs[13]);
    int integration_type = (int) mxGetScalar(prhs[14]);
    #if MX_HAS_INTERLEAVED_COMPLEX
        mxComplexDoubles *pc;
        pc = mxGetComplexDoubles(prhs[4]);
        gamma = (*pc).real + I*(*pc).imag);
        pc = mxGetComplexDoubles(prhs[5]);
        s = (*pc).real + I*(*pc).imag);
        pc = mxGetComplexDoubles(prhs[7]);
        kappa = (*pc).real + I*(*pc).imag);
    #else
        double *pr, *pi;
        //FIXME will crash if any of the arguments is REAL in Matlab side...
        pr = (double *) mxGetPr(prhs[4]);
        pi = (double *) mxGetPi(prhs[4]);
        gamma = (*pr) + I*(*pi);
        pr = (double *) mxGetPr(prhs[5]);
        pi = (double *) mxGetPi(prhs[5]);
        s = (*pr) + I*(*pi);
        pr = (double *) mxGetPr(prhs[7]);
        pi = (double *) mxGetPi(prhs[7]);
        kappa = (*pr) + I*(*pi);
        pr = (double *) mxGetPr(prhs[8]);
        pi = (double *) mxGetPi(prhs[8]);
        ref_l = (*pr) + I*(*pi);
        pr = (double *) mxGetPr(prhs[9]);
        pi = (double *) mxGetPi(prhs[9]);
        ref_t = (*pr) + I*(*pi);
    #endif
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
