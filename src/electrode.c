#include <math.h>
#include <float.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "electrode.h"
#include "auxiliary.h"
#include "cubature.h"

int
populate_electrode (Electrode* electrode, const double start_point[3],
                    const double end_point[3], double radius)
{
    if (equal_points(start_point, end_point)) {
        printf("Error: start_point and end_point are equal.\n");
        return 1;
    }
    if (radius <= 0) {
        printf("Error: radius < 0.\n");
        return 2;
    }
    electrode->radius = radius;
    electrode->length = vector_length(start_point, end_point);
    for (size_t i = 0; i < 3; i++) {
        electrode->start_point[i] = start_point[i];
        electrode->end_point[i] = end_point[i];
        electrode->middle_point[i] = (start_point[i] + end_point[i]) / 2.0;
    }
    return 0;
}

bool
equal_electrodes (const Electrode* sender, const Electrode* receiver)
{
    bool c1 = equal_points(sender->start_point, receiver->start_point);
    bool c2 = equal_points(sender->end_point, receiver->end_point);
    bool c3 = equal_points(sender->end_point, receiver->start_point);
    bool c4 = equal_points(sender->start_point, receiver->end_point);
    bool c5 = ( sender->radius - receiver->radius < DBL_EPSILON );
    return (c5 && ( (c1 && c2) || (c3 && c4) ));
}

int
electrodes_file (const char file_name[], Electrode* electrodes,
                 size_t num_electrodes)
{
    FILE* stream = fopen(file_name, "r");
    if (stream == NULL) {
        printf("Cannot open file %s\n", file_name);
        return -10;
    }
    double start_point[3], end_point[3];
    double radius;
    const int nvalues = 7;
    int success = nvalues;
    int pop_error = 0;
    for (size_t i = 0; i < num_electrodes; i++) {
        success = fscanf(stream, "%lf, %lf, %lf, %lf, %lf, %lf, %lf",
                         &start_point[0], &start_point[1], &start_point[2],
                         &end_point[0], &end_point[1], &end_point[2], &radius);
        if (success != nvalues) {
            printf("error reading line %i of file %s\n", (int) (i+1), file_name);
            break;
        }
        success = populate_electrode(&(electrodes[i]), start_point, end_point, radius);
        if (success != 0) {
            printf("Bad input: could not create electrode %i from file %s\n",
                   (int) (i+1), file_name);
            pop_error = success;
            break;
        }
    }
    fclose(stream);
    if (pop_error == 0) {
        return (success - nvalues);
    } else {
        return (pop_error);
    }
}

int
nodes_file (const char file_name[], double* nodes, size_t num_nodes)
{
    FILE* stream = fopen(file_name, "r");
    if (stream == NULL) {
        printf("Cannot open file %s\n", file_name);
        return -10;
    }
    int success = 3;
    for (size_t i = 0; i < num_nodes; i++) {
        success = fscanf(stream, "%lf, %lf, %lf", &nodes[i*3],
                         &nodes[i*3 + 1], &nodes[i*3 + 2]);
        if (success != 3) {
            printf("error reading line %i of file %s\n", (int) (i+1), file_name);
            break;
        }
    }
    fclose(stream);
    return (success - 3);
}

int
segment_electrode (Electrode* electrodes, double* nodes, size_t num_segments,
                   const double* start_point, const double* end_point, double radius)
{
    if (num_segments < 1) {
        printf("Error: number of segments should be greater than 0.\n");
        return 1;
    }
    size_t num_nodes = num_segments + 1;
    double x[num_nodes];
    double y[num_nodes];
    double z[num_nodes];
    linspace(start_point[0], end_point[0], num_nodes, x);
    linspace(start_point[1], end_point[1], num_nodes, y);
    linspace(start_point[2], end_point[2], num_nodes, z);
    for (size_t i = 0; i < num_nodes; i++) {
        nodes[i*3] = x[i];
        nodes[i*3 + 1] = y[i];
        nodes[i*3 + 2] = z[i];
    }
    for (size_t i = 0; i < num_segments; i++) {
        populate_electrode(&(electrodes[i]), (nodes + i*3), (nodes + i*3 + 3), radius);
    }
    return 0;
}

size_t
nodes_from_elecs (double* nodes, Electrode* electrodes, size_t num_electrodes)
{
    // Array of pointers that will have the maximum number of nodes possible.
    // Make them NULL and malloc for each new unique node.
    // Later copy them to nodes
    double* *dummy_nodes = malloc(2*num_electrodes * sizeof(double[3]));
    for (size_t i = 1; i < 2*num_electrodes; i++) {
        dummy_nodes[i] = NULL;
    }
    size_t num_nodes = 0;
    int node_present;
    for (size_t i = 0; i < num_electrodes; i++) {
        // start_point in nodes?
        node_present = 0;
        for (size_t k = 0; k < num_nodes; k++) {
            if (equal_points(electrodes[i].start_point, dummy_nodes[k])) {
                node_present = 1;
                break;
            }
        }
        if (!node_present) {
            dummy_nodes[num_nodes] = electrodes[i].start_point;
            num_nodes++;
        }
        // end_point in nodes?
        node_present = 0;
        for (size_t k = 0; k < num_nodes; k++) {
            if (equal_points(electrodes[i].end_point, dummy_nodes[k])) {
                node_present = 1;
                break;
            }
        }
        if (!node_present) {
            dummy_nodes[num_nodes] = electrodes[i].end_point;
            num_nodes++;
        }
    }
    for (size_t i = 0; i < num_nodes; i++) {
        for (size_t k = 0; k < 3; k++) {
            nodes[i*3 + k] = dummy_nodes[i][k];
        }
    }
    free(dummy_nodes);
    return num_nodes;
}

int
integrand_double (unsigned ndim, const double* t, void *auxdata, unsigned fdim,
                  double* fval)
{
    Integration_data* p = (Integration_data*) auxdata;
    const Electrode* sender = p->sender;
    const Electrode* receiver = p->receiver;
    _Complex double gamma = p->gamma;
    double point_r[3], point_s[3];
    for (size_t i = 0; i < 3; i++) {
        point_r[i] = t[0] * (receiver->end_point[i] - receiver->start_point[i])
                   + receiver->start_point[i];
        point_s[i] = t[1] * (sender->end_point[i] - sender->start_point[i])
                   + sender->start_point[i];
    }
    double r = vector_length(point_s, point_r);
    if (r < DBL_EPSILON) r = receiver->radius;
    _Complex double exp_gr = cexp(-gamma * r);
    fval[0] = creal(exp_gr) / r;
    fval[1] = cimag(exp_gr) / r;
    return 0;
}

int
integrand_single (unsigned ndim, const double* t, void *auxdata, unsigned fdim,
                  double* fval)
{
    Integration_data* p = (Integration_data*) auxdata;
    const Electrode* sender = p->sender;
    const Electrode* receiver = p->receiver;
    _Complex double gamma = p->gamma;
    double point_s[3];
    for (size_t i = 0; i < 3; i++) {
        point_s[i] = t[0] * (sender->end_point[i] - sender->start_point[i])
                   + sender->start_point[i];
    }
    double r = vector_length(point_s, receiver->middle_point);
    if (r < DBL_EPSILON) r = DBL_EPSILON;
    _Complex double exp_gr = cexp(-gamma * r);
    fval[0] = creal(exp_gr) / r;
    fval[1] = cimag(exp_gr) / r;
    return 0;
}

int
logNf (unsigned ndim, const double* t, void *auxdata, unsigned fdim, double* fval)
{
    Integration_data* p = (Integration_data*) auxdata;
    const Electrode* sender = p->sender;
    const Electrode* receiver = p->receiver;
    double point_r[3];
    for (size_t i = 0; i < 3; i++) {
        point_r[i] = t[0]*(receiver->end_point[i] - receiver->start_point[i])
                   + receiver->start_point[i];
    }
    double r1 = vector_length(point_r, sender->start_point);
    double r2 = vector_length(point_r, sender->end_point);
    double Nf = (r1 + r2 + sender->length)/(r1 + r2 - sender->length);
    fval[0] = log(Nf);
    return 0;
}

_Complex double
self_integral (const Electrode* sender)
{
    double ls = sender->length;
    double k1 = sender->radius/ls;
    double k2 = sqrt(1.0 + k1*k1);
    return ( 2.0 * ls * ( log( (k2 + 1.0) / k1 ) - k2 + k1 ) );
}

int
integral (const Electrode* sender, const Electrode* receiver,
          _Complex double gamma, size_t max_eval, double req_abs_error,
          double req_rel_error, int integration_type, double result[2],
          double error[2])
{
    Integration_data* auxdata = malloc(sizeof(Integration_data));
    auxdata->sender = sender;
    auxdata->receiver = receiver;
    auxdata->gamma = gamma;
    double tmin[] = {0.0, 0.0};
    double tmax[] = {1.0, 1.0};
    // Integration is done with a change of variable: from l in (0, L) to
    // t in (0, 1). Hence, dl = L*dt
    int failure = -1;
    result[0] = 0.0;
    result[1] = 0.0;
    error[0] = 0.0;
    error[1] = 0.0;
    switch (integration_type) {
        case INTG_NONE:
            result[0] = vector_length(sender->middle_point, receiver->middle_point);
            failure = 0;
            break;

        case INTG_DOUBLE:
            failure = hcubature(2, integrand_double, auxdata, 2, tmin, tmax,
                                max_eval, req_abs_error, req_rel_error,
                                ERROR_PAIRED, result, error);
            result[0] *= receiver->length * sender->length;
            result[1] *= receiver->length * sender->length;
            break;

        case INTG_SINGLE:
            failure = hcubature(2, integrand_single, auxdata, 1, tmin, tmax,
                                max_eval, req_abs_error, req_rel_error,
                                ERROR_PAIRED, result, error);
            result[0] *= sender->length;
            result[1] *= sender->length;
            break;

        case INTG_MHEM:
            failure = hcubature(1, logNf, auxdata, 1, tmin, tmax, max_eval,
                                req_abs_error, req_rel_error, ERROR_INDIVIDUAL,
                                result, error);
            result[0] *= receiver->length;
            break;

        default:
            failure = -10;
            break;
    }
    free(auxdata);
    if (failure) printf("integration failed with code: %i\n", failure);
    return failure;
}

int
calculate_impedances (_Complex double* zl, _Complex double* zt,
                      const Electrode* electrodes, size_t num_electrodes,
                      _Complex double gamma, _Complex double s, double mur,
                      _Complex double kappa, size_t max_eval, double req_abs_error,
                      double req_rel_error, int integration_type)
{
    double result[2], error[2], ls, lr, k1[3], k2, cost;
    _Complex double iwu_4pi, one_4pik;
    if (integration_type == INTG_MHEM || integration_type == INTG_NONE) {
        iwu_4pi = 1.0;
        one_4pik = 1.0;
    } else {
        iwu_4pi = s * mur * MU0 / (FOUR_PI);
        one_4pik = 1.0 / (FOUR_PI * kappa);
    }
    _Complex double intg = 0;
    int failure = 0;
    //#pragma omp parallel for private(ls, k1, intg, lr, cost, k2, failure)
    for (size_t i = 0; i < num_electrodes; i++) {
        ls = electrodes[i].length;
        /*if (integration_type == INTG_DOUBLE) { // FIXME wrong results
            Electrode sender;
            Electrode receiver;
            double start_point[3];
            double end_point[3];
            lr = electrodes[i].radius;
            start_point[0] = 0.0; start_point[1] = 0.0; start_point[2] = 0.0;
            end_point[0]   =  ls; end_point[1]   = 0.0; end_point[2]   = 0.0;
            populate_electrode(&sender, start_point, end_point, lr);
            start_point[0] = 0.0; start_point[1] = lr; start_point[2] = 0.0;
            end_point[0]   =  ls; end_point[1]   = lr; end_point[2]   = 0.0;
            populate_electrode(&receiver, start_point, end_point, lr);
            failure = integral(&sender, &receiver, gamma,
                               max_eval, req_abs_error, req_rel_error,
                               integration_type, result, error);
            if (failure != 0) return failure;
        } else {
            intg = self_integral(&(electrodes[i]));
        }
        printf("%f, ", cabs(intg));
        printf("%f\n", cabs(self_integral(&(electrodes[i]))));*/
        intg = self_integral(&(electrodes[i]));
        zl[i*num_electrodes + i] = iwu_4pi * intg;
        zt[i*num_electrodes + i] = one_4pik / (ls * ls) * intg;
        for (int m = 0; m < 3; m++) {
            k1[m] = (electrodes[i].end_point[m] - electrodes[i].start_point[m]);
        }
        for (size_t k = (i + 1); k < num_electrodes; k++) {
            lr = electrodes[k].length;
            cost = 0.0;
            for (size_t m = 0; m < 3; m++) {
                k2 = (electrodes[k].end_point[m] - electrodes[k].start_point[m]);
                cost += k1[m] * k2;
            }
            cost = fabs(cost);
            failure = integral(&(electrodes[i]), &(electrodes[k]), gamma,
                               max_eval, req_abs_error, req_rel_error,
                               integration_type, result, error);
            if (failure != 0) return failure; // TODO error check with OpenMP
            intg = result[0] + I*result[1];
            zt[i*num_electrodes + k] = one_4pik / (ls * lr) * intg;
            //zt[k*num_electrodes + i] = zt[i*num_electrodes + k];
            if (cost < DBL_EPSILON) {
                zl[i*num_electrodes + k] = 0.0;
            } else {
                cost = cost / (ls * lr);
                zl[i*num_electrodes + k] = iwu_4pi * intg * cost;
            }
            //zl[k*num_electrodes + i] = zl[i*num_electrodes + k];
        }
    }
    return 0;
}

int
impedances_images (_Complex double* zl, _Complex double* zt,
                   const Electrode* electrodes, const Electrode* images,
                   size_t num_electrodes,   _Complex double gamma,
                   _Complex double s, double mur, _Complex double kappa,
                   _Complex double ref_l, _Complex double ref_t,
                   size_t max_eval, double req_abs_error,
                   double req_rel_error, int integration_type)
{
    double result[2], error[2], ls, lr, k1[3], k2;
    double cost = 0;
    _Complex double iwu_4pi, one_4pik;
    if (integration_type == INTG_MHEM || integration_type == INTG_NONE) {
        iwu_4pi = 1.0;
        one_4pik = 1.0;
        ref_l = 1.0;
        ref_t = 1.0;
    } else {
        iwu_4pi = ref_l * s * mur * MU0 / (FOUR_PI);
        one_4pik = ref_t / (FOUR_PI * kappa);
    }
    _Complex double intg = 0;
    int failure = 0;
    //#pragma omp parallel for private(ls, k1, intg, lr, cost, k2, failure)
    for (size_t i = 0; i < num_electrodes; i++) {
        ls = electrodes[i].length;
        for (int m = 0; m < 3; m++) {
            k1[m] = (electrodes[i].end_point[m] - electrodes[i].start_point[m]);
        }
        for (size_t k = i; k < num_electrodes; k++) {
            lr = images[k].length;
            cost = 0.0;
            for (size_t m = 0; m < 3; m++) {
                k2 = (images[k].end_point[m] - images[k].start_point[m]);
                cost += k1[m] * k2;
            }
            cost = fabs(cost);
            failure = integral(&(electrodes[i]), &(images[k]), gamma, max_eval,
                               req_abs_error, req_rel_error, integration_type,
                               result, error);
            if (failure != 0) return failure; // TODO error check with OpenMP
            intg = result[0] + I*result[1];
            zt[i*num_electrodes + k] += one_4pik / (ls * lr) * intg;
            //zt[k*num_electrodes + i] = zt[i*num_electrodes + k];
            if (cost > DBL_EPSILON) {
                cost = cost / (ls * lr);
                zl[i*num_electrodes + k] += iwu_4pi * intg * cost;
                //zl[k*num_electrodes + i] = zl[i*num_electrodes + k];
            }
        }
    }
    return 0;
}

// TODO integration of fields and potentials using the middle_point
_Complex double
electric_potential (const double* point, const Electrode* electrodes,
                    size_t num_electrodes, const _Complex double* it,
                    _Complex double gamma, _Complex double kappa,
                    size_t max_eval, double req_abs_error, double req_rel_error)
{
    Electrode* point_elec = malloc(sizeof(Electrode));
    double result[2], error[2];
    _Complex double pot = 0.0;
    for (int i = 0; i < 3; i++) {
        point_elec->middle_point[i] = *(point + i);
    }
    int failure;
    for (size_t m = 0; m < num_electrodes; m++) {
        failure = integral(&(electrodes[m]), point_elec, gamma, max_eval,
                           req_abs_error, req_rel_error, INTG_SINGLE,
                           result, error);
        if (failure != 0) return failure; // TODO error check with OpenMP
        pot += (result[0] + I*result[1]) * it[m] / electrodes[m].length;
    }
    free(point_elec);
    return pot/(FOUR_PI*kappa);
}

int
magnetic_potential (const double* point, const Electrode* electrodes,
                    size_t num_electrodes, const _Complex double* il,
                    _Complex double gamma, double mur, size_t max_eval,
                    double req_abs_error, double req_rel_error,
                    _Complex double* va)
{
    int failure = 0;
    Electrode* point_elec = malloc(sizeof(Electrode));
    double result[2], error[2], dx, dy, dz;
    _Complex double pot = 0.0;
    for (int i = 0; i < 3; i++) {
        point_elec->middle_point[i] = *(point + i);
        va[i] = 0.0;
    }
    for (size_t m = 0; m < num_electrodes; m++) {
        if (cabs(il[m]) < DBL_EPSILON) continue;
        dx = electrodes[m].end_point[0] - electrodes[m].start_point[0];
        dy = electrodes[m].end_point[1] - electrodes[m].start_point[1];
        dz = electrodes[m].end_point[2] - electrodes[m].start_point[2];
        failure = integral(&(electrodes[m]), point_elec, gamma, max_eval,
                           req_abs_error, req_rel_error, INTG_SINGLE,
                           result, error);
        if (failure) return failure;
        // divide pot by L as to avoid dividing each dx, dy, dz by L
        pot = (result[0] + I*result[1]) * il[m] / electrodes[m].length;
        va[0] += pot * dx;
        va[1] += pot * dy;
        va[2] += pot * dz;
    }
    free(point_elec);
    _Complex double u_4pi = mur * MU0 / (FOUR_PI);
    for (int i = 0; i < 3; i++) va[i] *= u_4pi;
    return 0;
}

int
elec_field_integrand (unsigned ndim, const double* t, void *auxdata,
                      unsigned fdim, double* fval)
{
    Field_integrand_data* p = (Field_integrand_data*) auxdata;
    const double* point = p->point1;
    const Electrode* electrode = p->electrodes;
    _Complex double gamma = p->gamma;
    double r1, dx[3];
    double r2 = 0.0;
    for (size_t i = 0; i < 3; i++) {
        r1 = t[0] * (electrode->end_point[i] - electrode->start_point[i])
           + electrode->start_point[i];
        dx[i] = point[i] - r1;
        r2 += pow(dx[i], 2.0);
    }
    r1 = sqrt(r2);
    if (r2 < DBL_EPSILON) r2 = DBL_EPSILON;
    _Complex double efield = (1 + gamma * r1) * cexp(-gamma * r1) / r2;
    _Complex double ex;
    for (int i = 0; i < 3; i++) {
        ex = efield * dx[i] / r1;
        fval[2*i] = creal(ex);
        fval[2*i + 1] = cimag(ex);
    }
    return 0;
}

int
electric_field (const double* point, const Electrode* electrodes,
                size_t num_electrodes, const _Complex double* il,
                const _Complex double* it, _Complex double gamma,
                _Complex double s, double mur, _Complex double kappa,
                size_t max_eval, double req_abs_error, double req_rel_error,
                _Complex double* ve)
{
    //for (int i = 0; i < 3; i++) ve[i] = 0.0;
    if (cabs(s) > DBL_EPSILON) {
        magnetic_potential(point, electrodes, num_electrodes, il, gamma, mur,
                           max_eval, req_abs_error, req_rel_error, ve);
        for (int i = 0; i < 3; i++) {
            ve[i] *= -s;
        }
    }
    Field_integrand_data* auxdata = malloc(sizeof(Field_integrand_data));
    auxdata->point1 = point;
    auxdata->gamma = gamma;
    double tmin[] = {0.0};
    double tmax[] = {1.0};
    unsigned fdim = 6;
    double result[fdim], error[fdim];
    int failure = 0;
    _Complex double it_4pik;
    for (size_t m = 0; m < num_electrodes; m++) {
        if (cabs(it[m]) < DBL_EPSILON) continue;
        auxdata->electrodes = (electrodes + m);
        failure = hcubature(fdim, elec_field_integrand, auxdata, 1, tmin, tmax,
                            max_eval, req_abs_error, req_rel_error, ERROR_PAIRED,
                            result, error);
        it_4pik = it[m] / (FOUR_PI * kappa);
        for (int i = 0; i < 3; i++) {
            ve[i] += (result[2*i] + I*result[2*i + 1]) * it_4pik;
        }
    }
    free(auxdata);
    return failure;
}

int
v_mag_pot_integrand (unsigned ndim, const double *t, void *auxdata, unsigned fdim,
                     double *fval)
{
    Field_integrand_data* p = (Field_integrand_data*) auxdata;
    const double* point_1 = p->point1;
    const double* point_2 = p->point2;
    const Electrode* electrodes = p->electrodes;
    size_t num_electrodes = p->num_electrodes;
    const _Complex double* il = p->il;
    _Complex double gamma = p->gamma;
    double mur = p->mur;
    size_t max_eval = p->max_eval;
    double req_abs_error = p->req_abs_error;
    double req_rel_error = p->req_rel_error;
    double point[3];
    double dx[3];
    double length = 0.0;
    for (int i = 0; i < 3; i++) {
        dx[i] = point_2[i] - point_1[i];
        point[i] = t[0] * dx[i] + point_1[i];
        length += pow(dx[i], 2.0);
    }
    if (length < DBL_EPSILON) return 1;
    length = sqrt(length);
    _Complex double va[3];
    magnetic_potential(point, electrodes, num_electrodes, il, gamma, mur,
                       max_eval, req_abs_error, req_rel_error, va);
    _Complex double va_cost = 0.0;
    for (int i = 0; i < 3; i++) va_cost += va[i] * dx[i];
    va_cost = va_cost / length;
    fval[0] = creal(va_cost);
    fval[1] = cimag(va_cost);
    return 0;
}

int
v_elecf_integrand (unsigned ndim, const double *t, void *auxdata, unsigned fdim,
                  double *fval)
{
    Field_integrand_data* p = (Field_integrand_data*) auxdata;
    const double* point_1 = p->point1;
    const double* point_2 = p->point2;
    const Electrode* electrodes = p->electrodes;
    size_t num_electrodes = p->num_electrodes;
    const _Complex double* il = p->il;
    const _Complex double* it = p->it;
    _Complex double gamma = p->gamma;
    _Complex double s = p->s;
    double mur = p->mur;
    _Complex double kappa = p->kappa;
    size_t max_eval = p->max_eval;
    double req_abs_error = p->req_abs_error;
    double req_rel_error = p->req_rel_error;
    double point[3];
    double dx[3];
    double length = 0.0;
    for (int i = 0; i < 3; i++) {
        dx[i] = point_2[i] - point_1[i];
        point[i] = t[0] * dx[i] + point_1[i];
        length += pow(dx[i], 2.0);
    }
    if (length < DBL_EPSILON) return 1;
    length = sqrt(length);
    _Complex double ve[3];
    electric_field(point, electrodes, num_electrodes, il, it, gamma, s, mur,
                   kappa, max_eval, req_abs_error, req_rel_error, ve);
    _Complex double ve_cost = 0.0;
    for (int i = 0; i < 3; i++) ve_cost += ve[i] * dx[i];
    ve_cost /= length;
    fval[0] = creal(ve_cost);
    fval[1] = cimag(ve_cost);
    return 0;
}

_Complex double
voltage (const double* point1, const double* point2,
         const Electrode* electrodes, size_t num_electrodes,
         const _Complex double* il, const _Complex double* it,
         _Complex double gamma, _Complex double s, double mur,
         _Complex double kappa, size_t max_eval, double req_abs_error,
         double req_rel_error)
{
    Field_integrand_data* auxdata = malloc(sizeof(Field_integrand_data));
    auxdata->point1 = point1;
    auxdata->point2 = point2;
    auxdata->electrodes = electrodes;
    auxdata->num_electrodes = num_electrodes;
    auxdata->il = il;
    auxdata->it = it;
    auxdata->gamma = gamma;
    auxdata->s = s;
    auxdata->mur = mur;
    auxdata->kappa = kappa;
    auxdata->max_eval = max_eval;
    auxdata->req_abs_error = req_abs_error;
    auxdata->req_rel_error = req_rel_error;
    double tmin[] = {0.0};
    double tmax[] = {1.0};
    unsigned fdim = 2;
    double result[fdim], error[fdim];
    hcubature(fdim, v_elecf_integrand, auxdata, 1, tmin, tmax, max_eval,
              req_abs_error, req_rel_error, ERROR_PAIRED, result, error);
    free(auxdata);
    return (result[0] + I * result[1]);
    // TODO alternative: integrate electric_field directly
}
