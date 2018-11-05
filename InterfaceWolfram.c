/**
Funtions to interface the C code to Wolfram Mathematica.
*/
#include <WolframLibrary.h>
#include <Electrode.h>
#include <complex.h>

DLLEXPORT mint WolframLibrary_getVersion(){return WolframLibraryVersion;}
DLLEXPORT int WolframLibrary_initialize( WolframLibraryData libData) {return 0;}
DLLEXPORT void WolframLibrary_uninitialize( WolframLibraryData libData) {}

/**
Calculates de impedances of the electrode systems.
@param Args[0] Rank 2 tensor `9*num_electrodes`, each element is a flat list
describing an electrode as "{x0, y0, z0, x1, y1, z1, radius, Re(zi), Im(zi)}".
@param Args[1] complex wave propagation constant
@param Args[2] angular frequency in rad/s
@param Args[3] magnetic permeability of the medium
@param Args[4] medium complex conductivity `(sigma + j*w*eps)` in S/m
@param Args[5] specifies a maximum number of function evaluations (0 for no
limit)
@param Args[6] the absolute error requested (0 to ignore)
@param Args[7] the relative error requested (0 to ignore)
@param Args[8] (enumeration defined in cubature.h) error checking scheme
@param Args[9] integration_type type of integration to be done.
*/
DLLEXPORT int Mcalculate_impedances(
    WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res)
{
    int err; // error code
    MTensor electensor = MArgument_getMTensor(Args[0]);
    mcomplex gamma = MArgument_getComplex(Args[1]);
    mreal w = MArgument_getReal(Args[2]);
    mreal mu = MArgument_getReal(Args[3]);
    mcomplex kappa = MArgument_getComplex(Args[4]);
    mint max_eval = MArgument_getInteger(Args[5]);
    mreal req_abs_error = MArgument_getReal(Args[6]);
    mreal req_rel_error = MArgument_getReal(Args[7]);
    mint error_norm = MArgument_getInteger(Args[8]);
    mint integration_type = MArgument_getInteger(Args[9]);
    mint const* dims;
    dims = libData->MTensor_getDimensions(electensor);
    mreal* elecdata;
    elecdata = libData->MTensor_getRealData(electensor);
    int num_electrodes = dims[0];
    Electrode* electrodes = (Electrode*) malloc(sizeof(Electrode)*num_electrodes);
    double start_point[3], end_point[3];
    double radius;
    _Complex double zi;
    for (int i = 0; i < num_electrodes; i++)
    {
        start_point[0] = elecdata[9*i + 0];
        start_point[1] = elecdata[9*i + 1];
        start_point[2] = elecdata[9*i + 2];
        end_point[0] = elecdata[9*i + 3];
        end_point[1] = elecdata[9*i + 4];
        end_point[2] = elecdata[9*i + 5];
        radius = elecdata[9*i + 6];
        zi = elecdata[9*i + 7] + I*elecdata[9*i + 8];
        populate_electrode(&(electrodes[i]), start_point, end_point, radius, zi);
    }
    int ne2 = num_electrodes*num_electrodes;
    _Complex double* zl = (_Complex double*) malloc(sizeof(_Complex double)*ne2);
    _Complex double* zt = (_Complex double*) malloc(sizeof(_Complex double)*ne2);
    _Complex double gamma1 = gamma.ri[0] + I*gamma.ri[1];
    _Complex double kappa1 = kappa.ri[0] + I*kappa.ri[1];
    double w1 = w;
    double mu1 = mu;
    int max_eval1 = max_eval;
    int error_norm1 = error_norm;
    int integration_type1 = integration_type;
    double req_abs_error1 = req_abs_error;
    double req_rel_error1 = req_rel_error;
    int ne = num_electrodes;
    calculate_impedances(electrodes, ne, zl, zt, gamma1, w1, mu1,
            kappa1, max_eval1, req_abs_error1, req_rel_error1, error_norm1,
            integration_type1);
    free(electrodes);
    free(zl);
    free(zt);

    MTensor zlzt;
    mreal *data;
    mint out_dim[] = {4*ne2};
    err = libData->MTensor_new(MType_Real, 1, out_dim, &zlzt);
    data = libData->MTensor_getRealData(zlzt);
    for (int i = 0; i < ne2; i++) {
        data[i] = creal(zl[i]);
    }
    for (int i = 0; i < ne2; i++) {
        data[i + ne2] = cimag(zl[i]);
    }
    for (int i = 0; i < ne2; i++) {
        data[i + 2*ne2] = creal(zt[i]);
    }
    for (int i = 0; i < ne2; i++) {
        data[i + 3*ne2] = cimag(zt[i]);
    }
    MArgument_setMTensor(Res, zlzt);
    return LIBRARY_NO_ERROR;
}
