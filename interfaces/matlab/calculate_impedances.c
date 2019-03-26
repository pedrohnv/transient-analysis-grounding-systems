#include "interface_matlab.h"
#include "electrode.h"
#include "mex.h"
//mex calculate_impedances.c InterfaceMatlab.c ..\\..\\cubature\\hcubature.c -I. -I..\\..\\cubature

/**
Calculate the impedance matrices of the electrode system.
TODO describe args
*/
void
mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //TODO check inputs and print useful error messages
    //mxGetNumberOfElements = N_array_elements + N_struct_fields
    size_t num_electrodes = mxGetNumberOfElements(prhs[0])/12;
    Electrode *electrodes = malloc(sizeof(Electrode)*num_electrodes);
    for (size_t i = 0; i < num_electrodes; i++) {
        cast_electrode(prhs[0], i, &(electrodes[i]));
    }
    _Complex double gamma; //prhs[1]
    _Complex double s; //prhs[2]
    double mur = (double) mxGetScalar(prhs[3]);
    _Complex double kappa; //prhs[4]
    size_t max_eval = (size_t) mxGetScalar(prhs[5]);
    double req_abs_error = (double) mxGetScalar(prhs[6]);
    double req_rel_error = (double) mxGetScalar(prhs[7]);
    int error_norm = (int) mxGetScalar(prhs[8]);
    int integration_type = (int) mxGetScalar(prhs[9]);
    #if MX_HAS_INTERLEAVED_COMPLEX
        mxComplexDoubles *pc;
        pc = mxGetComplexDoubles(prhs[1]);
        gamma = (*pc).real + I*(*pc).imag);
        pc = mxGetComplexDoubles(prhs[2]);
        s = (*pc).real + I*(*pc).imag);
        pc = mxGetComplexDoubles(prhs[4]);
        kappa = (*pc).real + I*(*pc).imag);
    #else
        double *pr, *pi;
        //FIXME will crash if any of the arguments is REAL in Matlab side...
        pr = (double *) mxGetPr(prhs[1]);
        pi = (double *) mxGetPi(prhs[1]);
        gamma = (*pr) + I*(*pi);
        pr = (double *) mxGetPr(prhs[2]);
        pi = (double *) mxGetPi(prhs[2]);
        s = (*pr) + I*(*pi);
        pr = (double *) mxGetPr(prhs[4]);
        pi = (double *) mxGetPi(prhs[4]);
        kappa = (*pr) + I*(*pi);
    #endif
    size_t ne2 = num_electrodes*num_electrodes;
    _Complex double *zl = malloc(sizeof(_Complex double)*ne2);
    _Complex double *zt = malloc(sizeof(_Complex double)*ne2);
    calculate_impedances(electrodes, num_electrodes, zl, zt, gamma, s, mur,
                         kappa, max_eval, req_abs_error, req_rel_error,
                         error_norm, integration_type);
    plhs[0] = mxCreateDoubleMatrix((mwSize)num_electrodes, (mwSize)num_electrodes, mxCOMPLEX);
    plhs[1] = mxCreateDoubleMatrix((mwSize)num_electrodes, (mwSize)num_electrodes, mxCOMPLEX);
    #if MX_HAS_INTERLEAVED_COMPLEX
        mxComplexDoubles *zldata, *ztdata;
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
