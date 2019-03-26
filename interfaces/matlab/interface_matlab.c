#include "interface_matlab.h"
#include "mex.h"
#include "electrode.h"
//#include <math.h>
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
                mxComplexDoubles *pc;
                pc = mxGetComplexDouble(field_array_ptr);
                electrode->zi = (*pc).real + (*pc).imag*I;
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
    /*mexPrintf("start: %f, ", electrode->start_point[0]);
    mexPrintf("%f, ", electrode->start_point[1]);
    mexPrintf("%f\n", electrode->start_point[2]);
    mexPrintf("end: %f, ", electrode->end_point[0]);
    mexPrintf("%f, ", electrode->end_point[1]);
    mexPrintf("%f\n", electrode->end_point[2]);
    mexPrintf("middle: %f, ", electrode->middle_point[0]);
    mexPrintf("%f, ", electrode->middle_point[1]);
    mexPrintf("%f\n", electrode->middle_point[2]);
    mexPrintf("radius: %f\n", electrode->radius);
    mexPrintf("length: %f\n", electrode->length);
    mexPrintf("zi: %f + I*%f\n\n", creal(electrode->zi), cimag(electrode->zi));*/
    return 0;
}
/* FIXME
//auxiliary.c
int equal_points(const double point_1[3], const double point_2[3]) {
    for (size_t i = 0; i < 3; i++) {
        if (fabs(point_1[i] - point_2[i]) > FLT_EPSILON) return 0;
    }
    return 1;
}

double vector_norm(const double start_point[3], const double end_point[3]) {
    double length = 0.0;
    for (int i = 0; i < 3; i++) {
        length += pow(start_point[i] - end_point[i], 2.0);
    }
    length = sqrt(length);
    return length;
}

//Electrode.c
int integrand_double(unsigned ndim, const double *t, void *auxdata,
                     unsigned fdim, double *fval) {
    Integration_data *p = (Integration_data*) auxdata;
    const Electrode* sender = p->sender;
    const Electrode* receiver = p->receiver;
    _Complex double gamma = p->gamma;
    double point_r, point_s, r = 0.0;
    for (size_t i = 0; i < 3; i++) {
        point_r = t[0]*(receiver->end_point[i] - receiver->start_point[i]);
        point_r += receiver->start_point[i];
        point_s = t[1]*(sender->end_point[i] - sender->start_point[i]);
        point_s += sender->start_point[i];
        r += pow(point_r - point_s, 2.0);
    }
    r = sqrt(r);
    _Complex double exp_gr = cexp(-gamma*r);
    fval[0] = creal(exp_gr)/r;
    fval[1] = cimag(exp_gr)/r;
    return 0;
}

int exp_logNf(unsigned ndim, const double *t, void *auxdata, unsigned fdim,
              double *fval) {
    Integration_data *p = (Integration_data*) auxdata;
    const Electrode* sender = p->sender;
    const Electrode* receiver = p->receiver;
    double r1, r2, eta, point_r;
    r1 = 0.0;
    r2 = 0.0;
    eta = 0.0;
    for (size_t i = 0; i < 3; i++) {
        point_r = t[0]*(receiver->end_point[i] - receiver->start_point[i]);
        point_r += receiver->start_point[i];
        r1 += pow(point_r - sender->start_point[i], 2.0);
        r2 += pow(point_r - sender->end_point[i], 2.0);
        eta += pow(point_r - sender->middle_point[i], 2.0);
    }
    r1 = sqrt(r1);
    r2 = sqrt(r2);
    eta = sqrt(eta);
    _Complex double exp_gr;
    if (p->simplified) {
        exp_gr = 1.0;
    } else {
        exp_gr = cexp(- p->gamma*eta);
    }
    double Nf = (r1 + r2 + sender->length)/(r1 + r2 - sender->length);
    exp_gr = exp_gr*log(Nf);
    fval[0] = creal(exp_gr);
    fval[1] = cimag(exp_gr);
    return 0;
}

int integral(const Electrode* sender, const Electrode* receiver,
             _Complex double gamma, size_t max_eval, double req_abs_error,
             double req_rel_error, int error_norm, int integration_type,
             double result[2], double error[2]) {
    Integration_data* auxdata = malloc(sizeof(Integration_data));
    auxdata->sender = sender;
    auxdata->receiver = receiver;
    auxdata->gamma = gamma;
    double tmin[] = {0.0, 0.0};
    double tmax[] = {1.0, 1.0};
    int failure = 1;
    double rbar, lslr;
    _Complex double exp_gr;
    switch (integration_type) {
        case INTG_NONE:
            lslr = (sender->length)*(receiver->length);
            rbar = vector_norm(sender->middle_point, receiver->middle_point);
            exp_gr = cexp(-gamma*rbar)/rbar*lslr;
            result[0] = creal(exp_gr);
            result[1] = cimag(exp_gr);
            failure = 0;
            break;
        case INTG_DOUBLE:
            failure = hcubature(2, integrand_double, auxdata, 2, tmin, tmax,
                                max_eval, req_abs_error, req_rel_error,
                                error_norm, result, error);
            result[0] = result[0] * receiver->length * sender->length;
            result[1] = result[1] * receiver->length * sender->length;
            break;
        case INTG_EXP_LOGNF:
            auxdata->simplified = 0;
            failure = hcubature(2, exp_logNf, auxdata, 1, tmin, tmax, max_eval,
                                req_abs_error, req_rel_error, error_norm,
                                result, error);
            result[0] = result[0] * receiver->length;
            result[1] = result[1] * receiver->length;
            break;
        case INTG_LOGNF:
            auxdata->simplified = 1;
            failure = hcubature(1, exp_logNf, auxdata, 1, tmin, tmax, max_eval,
                                req_abs_error, req_rel_error, error_norm,
                                result, error);
            rbar = vector_norm(sender->middle_point, receiver->middle_point);
            exp_gr = cexp(-gamma*rbar);
            result[0] = creal(exp_gr) * receiver->length;
            result[1] = cimag(exp_gr) * receiver->length;
        break;
    }
    free(auxdata);
    return failure;
}

int calculate_impedances(const Electrode* electrodes, size_t num_electrodes,
                         _Complex double *zl, _Complex double *zt,
                         _Complex double gamma, _Complex double s, double mur,
                         _Complex double kappa, size_t max_eval,
                         double req_abs_error, double req_rel_error,
                         int error_norm, int integration_type) {
    double result[2], error[2], ls, lr, k1[3], k2, cost;
    _Complex double iwu_4pi = s*mur*MU0/(FOUR_PI);
    _Complex double one_4pik = 1.0/(FOUR_PI*kappa);
    _Complex double intg;
    size_t i, k, m;
    int failure;
    // _self and _mutual impedances are not used to reduce the number of
    // operations, as some of them would be done twice or more unnecessarily
    for (i = 0; i < num_electrodes; i++) {
        ls = electrodes[i].length;
        k1[0] = electrodes[i].radius/ls;
        k2 = sqrt(1.0 + k1[0]*k1[0]);
        cost = 2.0*(log( (k2 + 1.)/k1[0] ) - k2 + k1[0]);
        zl[i*num_electrodes + i] = iwu_4pi*ls*cost + electrodes[i].zi;
        zt[i*num_electrodes + i] = one_4pik/ls*cost;
        for (m = 0; m < 3; m++) {
            k1[m] = (electrodes[i].end_point[m] - electrodes[i].start_point[m]);
        }
        for (k = i+1; k < num_electrodes; k++) {
            lr = electrodes[k].length;
            cost = 0.0;
            for (m = 0; m < 3; m++) {
                k2 = (electrodes[k].end_point[m] - electrodes[k].start_point[m]);
                cost += k1[m]*k2;
            }
            cost = abs(cost/(ls*lr));
            failure = integral(&(electrodes[i]), &(electrodes[k]), gamma,
                               max_eval, req_abs_error, req_rel_error,
                               error_norm, integration_type, result, error);
            if (failure) return failure;
            intg = result[0] + I*result[1];
            zl[i*num_electrodes + k] = iwu_4pi*intg*cost;
            zt[i*num_electrodes + k] = one_4pik/(ls*lr)*intg;

            zl[k*num_electrodes + i] = zl[i*num_electrodes + k];
            zt[k*num_electrodes + i] = zt[i*num_electrodes + k];
        }
    }
    return 0;
}

int impedances_images(const Electrode* electrodes, const Electrode* images,
                      size_t num_electrodes, _Complex double *zl,
                      _Complex double *zt, _Complex double gamma,
                      _Complex double s, double mur, _Complex double kappa,
                      _Complex double ref_l, _Complex double ref_t,
                      size_t max_eval, double req_abs_error,
                      double req_rel_error, int error_norm, int integration_type) {
    double result[2], error[2], ls, lr, k1[3], k2, cost;
    _Complex double iwu_4pi = s*mur*MU0/(FOUR_PI);
    _Complex double one_4pik = 1.0/(FOUR_PI*kappa);
    _Complex double intg;
    int failure;
    // _mutual impedances are not used to reduce the number of
    // operations, as some of them would be done twice or more unnecessarily
    for (size_t i = 0; i < num_electrodes; i++) {
        ls = electrodes[i].length;
        for (size_t m = 0; m < 3; m++) {
            k1[m] = (electrodes[i].end_point[m] - electrodes[i].start_point[m]);
        }
        for (size_t k = i; k < num_electrodes; k++) {
            lr = images[k].length;
            cost = 0.0;
            for (size_t m = 0; m < 3; m++) {
                k2 = (images[k].end_point[m] - images[k].start_point[m]);
                cost += k1[m]*k2;
            }
            cost = abs(cost/(ls*lr));
            failure = integral(&(electrodes[i]), &(images[k]), gamma, max_eval,
                               req_abs_error, req_rel_error, error_norm,
                               integration_type, result, error);
            if (failure) return failure;
            intg = result[0] + I*result[1];
            zl[i*num_electrodes + k] += ref_l*iwu_4pi*intg*cost;
            zt[i*num_electrodes + k] += ref_t*one_4pik/(ls*lr)*intg;

            zl[k*num_electrodes + i] = zl[i*num_electrodes + k];
            zt[k*num_electrodes + i] = zt[i*num_electrodes + k];
        }
    }
    return 0;
}
*/
