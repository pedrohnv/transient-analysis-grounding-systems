/**
Funtions to interface the C code to Wolfram Mathematica.

FIXME in notebook files:
make sure Mathematica can see all dependencies (MKL and SLATEC)
*/
#include <WolframLibrary.h>
#include <Electrode.h>
#include <auxiliary.h>
#include <complex.h>

DLLEXPORT mint WolframLibrary_getVersion(){return WolframLibraryVersion;}
DLLEXPORT int WolframLibrary_initialize( WolframLibraryData libData) {return 0;}
DLLEXPORT void WolframLibrary_uninitialize( WolframLibraryData libData) {}

/** Mcalculate_impedances
Calculates the impedances of the electrode systems.
@param Args[0] Rank 2 tensor `9*num_electrodes`, each element is a flat list
describing an electrode as "{x0, y0, z0, x1, y1, z1, radius, Re(zi), Im(zi)}".
@param Args[1] complex wave propagation constant
@param Args[2] complex angular frequency `c + jw` (rad/s)
@param Args[3] relative magnetic permeability of the medium
@param Args[4] medium complex conductivity `(sigma + j*w*eps)` in S/m
@param Args[5] specifies a maximum number of function evaluations (0 for no
limit)
@param Args[6] the absolute error requested (0 to ignore)
@param Args[7] the relative error requested (0 to ignore)
@param Args[8] (enumeration defined in cubature.h) error checking scheme
@param Args[9] integration_type type of integration to be done.
*/
DLLEXPORT int Mcalculate_impedances(
    WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res) {
    int err; // error code
    MTensor elec_tensor = MArgument_getMTensor(Args[0]);
    mcomplex gamma = MArgument_getComplex(Args[1]);
    mcomplex s = MArgument_getComplex(Args[2]);
    mreal mur = MArgument_getReal(Args[3]);
    mcomplex kappa = MArgument_getComplex(Args[4]);
    mint max_eval = MArgument_getInteger(Args[5]);
    mreal req_abs_error = MArgument_getReal(Args[6]);
    mreal req_rel_error = MArgument_getReal(Args[7]);
    mint error_norm = MArgument_getInteger(Args[8]);
    mint integration_type = MArgument_getInteger(Args[9]);
    mint const* dims;
    dims = libData->MTensor_getDimensions(elec_tensor);
    mreal* elec_data;
    elec_data = libData->MTensor_getRealData(elec_tensor);
    size_t num_electrodes = dims[0];
    Electrode* electrodes = (Electrode*) malloc(num_electrodes*sizeof(Electrode));
    double start_point[3], end_point[3];
    double radius;
    _Complex double zi;
    for (size_t i = 0; i < num_electrodes; i++) {
        start_point[0] = elec_data[9*i + 0];
        start_point[1] = elec_data[9*i + 1];
        start_point[2] = elec_data[9*i + 2];
        end_point[0] = elec_data[9*i + 3];
        end_point[1] = elec_data[9*i + 4];
        end_point[2] = elec_data[9*i + 5];
        radius = elec_data[9*i + 6];
        zi = elec_data[9*i + 7] + I*elec_data[9*i + 8];
        populate_electrode(&(electrodes[i]), start_point, end_point, radius, zi);
    }
    size_t ne2 = num_electrodes*num_electrodes;
    _Complex double* zl = malloc(ne2*sizeof(_Complex double));
    _Complex double* zt = malloc(ne2*sizeof(_Complex double));
    _Complex double gamma1 = gamma.ri[0] + I*gamma.ri[1];
    _Complex double kappa1 = kappa.ri[0] + I*kappa.ri[1];
    _Complex double s1 = s.ri[0] + I*s.ri[1];
    calculate_impedances(electrodes, num_electrodes, zl, zt, gamma1, s1, mur,
            kappa1, max_eval, req_abs_error, req_rel_error, error_norm,
            integration_type);
    MTensor zlzt;
    mcomplex *data;
    mint out_dim[] = {2*ne2};
    err = libData->MTensor_new(MType_Complex, 1, out_dim, &zlzt);
    data = libData->MTensor_getComplexData(zlzt);
    for (size_t i = 0; i < ne2; i++) {
        data[i].ri[0] = creal(zl[i]);
        data[i].ri[1] = cimag(zl[i]);
    }
    for (size_t i = 0; i < ne2; i++) {
        data[i + ne2].ri[0] = creal(zt[i]);
        data[i + ne2].ri[1] = cimag(zt[i]);
    }
    MArgument_setMTensor(Res, zlzt);
    free(electrodes);
    free(zl);
    free(zt);
    return LIBRARY_NO_ERROR;
}

/** Mimpedances_images
Calculates the mutual impedances of the electrode systems and its image.
@param Args[0] Rank 2 tensor `9*num_electrodes`, each element is a flat list
describing an electrode as "{x0, y0, z0, x1, y1, z1, radius, Re(zi), Im(zi)}".
@param Args[1] complex wave propagation constant
@param Args[2] complex angular frequency `c + jw` (rad/s)
@param Args[3] relative magnetic permeability of the medium
@param Args[4] medium complex conductivity `(sigma + j*w*eps)` in S/m
@param Args[5] specifies a maximum number of function evaluations (0 for no
limit)
@param Args[6] the absolute error requested (0 to ignore)
@param Args[7] the relative error requested (0 to ignore)
@param Args[8] (enumeration defined in cubature.h) error checking scheme
@param Args[9] integration_type type of integration to be done.
@param Args[10] longitudinal current reflection coefficient
@param Args[11] transversal current reflection coefficient
@param Args[12] Rank 2 tensor `9*num_electrodes`, each element is a flat list
describing the image of an electrode as
"{x0, y0, z0, x1, y1, z1, radius, Re(zi), Im(zi)}".
*/
DLLEXPORT int Mimpedances_images(
    WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res) {
    int err; // error code
    MTensor elec_tensor = MArgument_getMTensor(Args[0]);
    mcomplex gamma = MArgument_getComplex(Args[1]);
    mcomplex s = MArgument_getComplex(Args[2]);
    mreal mur = MArgument_getReal(Args[3]);
    mcomplex kappa = MArgument_getComplex(Args[4]);
    mint max_eval = MArgument_getInteger(Args[5]);
    mreal req_abs_error = MArgument_getReal(Args[6]);
    mreal req_rel_error = MArgument_getReal(Args[7]);
    mint error_norm = MArgument_getInteger(Args[8]);
    mint integration_type = MArgument_getInteger(Args[9]);
    mcomplex ref_l = MArgument_getComplex(Args[10]);
    mcomplex ref_t = MArgument_getComplex(Args[11]);
    MTensor imag_tensor = MArgument_getMTensor(Args[12]); //TODO error checking
    mint const* dims;
    dims = libData->MTensor_getDimensions(elec_tensor);
    mreal *elec_data, *imag_data;
    elec_data = libData->MTensor_getRealData(elec_tensor);
    imag_data = libData->MTensor_getRealData(imag_tensor);
    size_t num_electrodes = dims[0];
    Electrode* electrodes = (Electrode*) malloc(sizeof(Electrode)*num_electrodes);
    Electrode* images = (Electrode*) malloc(sizeof(Electrode)*num_electrodes);
    double start_point[3], end_point[3];
    double radius;
    _Complex double zi;
    for (size_t i = 0; i < num_electrodes; i++) {
        start_point[0] = elec_data[9*i + 0];
        start_point[1] = elec_data[9*i + 1];
        start_point[2] = elec_data[9*i + 2];
        end_point[0] = elec_data[9*i + 3];
        end_point[1] = elec_data[9*i + 4];
        end_point[2] = elec_data[9*i + 5];
        radius = elec_data[9*i + 6];
        zi = elec_data[9*i + 7] + I*elec_data[9*i + 8];
        populate_electrode(&(electrodes[i]), start_point, end_point, radius, zi);
        start_point[0] = imag_data[9*i + 0];
        start_point[1] = imag_data[9*i + 1];
        start_point[2] = imag_data[9*i + 2];
        end_point[0] = imag_data[9*i + 3];
        end_point[1] = imag_data[9*i + 4];
        end_point[2] = imag_data[9*i + 5];
        radius = imag_data[9*i + 6];
        zi = imag_data[9*i + 7] + I*imag_data[9*i + 8];
        populate_electrode(&(images[i]), start_point, end_point, radius, zi);
    }
    size_t ne2 = num_electrodes*num_electrodes;
    _Complex double* zl = calloc(ne2, sizeof(_Complex double));
    _Complex double* zt = calloc(ne2, sizeof(_Complex double));
    _Complex double gamma1 = gamma.ri[0] + I*gamma.ri[1];
    _Complex double kappa1 = kappa.ri[0] + I*kappa.ri[1];
    _Complex double ref_l1 = ref_l.ri[0] + I*ref_l.ri[1];
    _Complex double ref_t1 = ref_t.ri[0] + I*ref_t.ri[1];
    _Complex double s1 = s.ri[0] + I*s.ri[1];
    impedances_images(electrodes, images, num_electrodes, zl, zt, gamma1, s1,
                      mur, kappa1, ref_l1, ref_t1, max_eval, req_abs_error,
                      req_rel_error, error_norm, integration_type);
    MTensor zlzt;
    mcomplex *data;
    mint out_dim[] = {2*ne2};
    err = libData->MTensor_new(MType_Complex, 1, out_dim, &zlzt);
    data = libData->MTensor_getComplexData(zlzt);
    for (size_t i = 0; i < ne2; i++) {
        data[i].ri[0] = creal(zl[i]);
        data[i].ri[1] = cimag(zl[i]);
    }
    for (size_t i = 0; i < ne2; i++) {
        data[i + ne2].ri[0] = creal(zt[i]);
        data[i + ne2].ri[1] = cimag(zt[i]);
    }
    MArgument_setMTensor(Res, zlzt);
    free(electrodes);
    free(images);
    free(zl);
    free(zt);
    return LIBRARY_NO_ERROR;
}

/** Mharmonic_impedance

*/
DLLEXPORT int Mharmonic_impedance1(
    WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res) {
    int err; // error code
    MTensor elec_tensor = MArgument_getMTensor(Args[0]);
    MTensor imag_tensor = MArgument_getMTensor(Args[1]);
    MTensor s_array = MArgument_getMTensor(Args[2]);
    mreal rsource = MArgument_getReal(Args[3]);
    MTensor mkappa1 = MArgument_getMTensor(Args[4]);
    MTensor mkappa2 = MArgument_getMTensor(Args[5]);
    MTensor mgamma1 = MArgument_getMTensor(Args[6]);
    mint max_eval = MArgument_getInteger(Args[7]);
    mreal req_abs_error = MArgument_getReal(Args[8]);
    mreal req_rel_error = MArgument_getReal(Args[9]);
    mint error_norm = MArgument_getInteger(Args[10]);
    MTensor nodes_tensor = MArgument_getMTensor(Args[11]);
    mint const* dims;
    dims = libData->MTensor_getDimensions(s_array);
    size_t ns = dims[0];
    dims = libData->MTensor_getDimensions(elec_tensor);
    size_t num_electrodes = dims[0];
    mreal *elec_data, *imag_data;
    elec_data = libData->MTensor_getRealData(elec_tensor);
    imag_data = libData->MTensor_getRealData(imag_tensor);
    Electrode* electrodes = (Electrode*) malloc(sizeof(Electrode)*num_electrodes);
    Electrode* images = (Electrode*) malloc(sizeof(Electrode)*num_electrodes);
    double start_point[3], end_point[3];
    double radius;
    for (size_t i = 0; i < num_electrodes; i++) {
        start_point[0] = elec_data[9*i + 0];
        start_point[1] = elec_data[9*i + 1];
        start_point[2] = elec_data[9*i + 2];
        end_point[0] = elec_data[9*i + 3];
        end_point[1] = elec_data[9*i + 4];
        end_point[2] = elec_data[9*i + 5];
        radius = elec_data[9*i + 6];
        populate_electrode(&(electrodes[i]), start_point, end_point, radius, 0.0);
        start_point[0] = imag_data[9*i + 0];
        start_point[1] = imag_data[9*i + 1];
        start_point[2] = imag_data[9*i + 2];
        end_point[0] = imag_data[9*i + 3];
        end_point[1] = imag_data[9*i + 4];
        end_point[2] = imag_data[9*i + 5];
        radius = imag_data[9*i + 6];
        populate_electrode(&(images[i]), start_point, end_point, radius, 0.0);
    }
    mreal* nodes_data;
    nodes_data = libData->MTensor_getRealData(nodes_tensor);
    dims = libData->MTensor_getDimensions(nodes_tensor);
    size_t num_nodes = dims[0];
    //double* nodes = (double*) malloc(sizeof(double)*num_nodes);
    double nodes[num_nodes][3];
    for (size_t i = 0; i < num_nodes; i++) {
        nodes[i][0] = nodes_data[3*i + 0];
        nodes[i][1] = nodes_data[3*i + 1];
        nodes[i][2] = nodes_data[3*i + 2];
    }
    _Complex double* kappa1 = malloc(ns*sizeof(_Complex double));
    _Complex double* kappa2 = malloc(ns*sizeof(_Complex double));
    _Complex double* gamma1 = malloc(ns*sizeof(_Complex double));
    _Complex double* s1 = malloc(ns*sizeof(_Complex double));
    mcomplex *kappa1_data, *kappa2_data, *gamma1_data, *s_data;
    kappa1_data = libData->MTensor_getComplexData(mkappa1);
    kappa2_data = libData->MTensor_getComplexData(mkappa2);
    gamma1_data = libData->MTensor_getComplexData(mgamma1);
    s_data = libData->MTensor_getComplexData(s_array);
    for (size_t i = 0; i < ns; i++) {
        kappa1[i] = kappa1_data[i].ri[0] + I*kappa1_data[i].ri[1];
        kappa2[i] = kappa2_data[i].ri[0] + I*kappa2_data[i].ri[1];
        gamma1[i] = gamma1_data[i].ri[0] + I*gamma1_data[i].ri[1];
        s1[i] = s_data[i].ri[0] + I*s_data[i].ri[1];
    }
    _Complex double* zh = malloc(ns*sizeof(_Complex double));
    harmonic_impedance1(
        ns, s1, kappa1, kappa2, gamma1, electrodes, images, num_electrodes,
        nodes, num_nodes, max_eval, req_abs_error, req_rel_error, error_norm,
        rsource, zh);
    MTensor mzh;
    mcomplex *zh_data;
    mint out_dim[] = {ns};
    err = libData->MTensor_new(MType_Complex, 1, out_dim, &mzh);
    zh_data = libData->MTensor_getComplexData(mzh);
    for (size_t i = 0; i < ns; i++) {
        zh_data[i].ri[0] = creal(zh[i]);
        zh_data[i].ri[1] = cimag(zh[i]);
    }
    MArgument_setMTensor(Res, mzh);
    free(electrodes);
    free(images);
    free(kappa1);
    free(kappa2);
    free(gamma1);
    free(s1);
    free(zh);
    return LIBRARY_NO_ERROR;
}
