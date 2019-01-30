#include "InterfaceOctave.hpp"

// double array
double* new_doublep(size_t n)
{
    return (new double[n]);
}

int delete_doublep(double* arr)
{
    delete[] arr;
    return 0;
}

int doublep_assign(double* arr, size_t pos, double val)
{
    arr[pos] = val;
    return 0;
}

double doublep_value(double* arr, size_t pos)
{
    return arr[pos];
}

// complex array
std::complex<double>* new_complexp(size_t n)
{
    return (new std::complex<double>[n]);
}

int delete_complexp(std::complex<double>* z)
{
    delete[] z;
    return 0;
}

int complexp_assign(std::complex<double>* z, size_t pos, std::complex<double> val)
{
    z[pos] = val;
    return 0;
}

std::complex<double> complexp_value(std::complex<double>* z, size_t pos)
{
    return z[pos];
}

// Electrode array
Electrode* new_electrodep(size_t n)
{
    return (new Electrode[n]);
}

int delete_electrodep(Electrode* electrodes)
{
    delete[] electrodes;
    return 0;
}

int electrodep_assign(Electrode* electrodes, size_t pos,
                      double start_point[3], double end_point[3],
                      double radius, std::complex<double> zi)
{
    _Complex double* zi_cast = (_Complex double*) (&zi);
    return populate_electrode(&(electrodes[pos]), start_point, end_point,
                              radius, *zi_cast);
}

Electrode electrodep_value(Electrode* electrodes, size_t pos)
{
    return electrodes[pos];
}

int Ocalculate_impedances(Electrode* electrodes, size_t num_electrodes,
                          std::complex<double>* zl, std::complex<double>* zt,
                          std::complex<double> gamma, std::complex<double> s,
                          double mur, std::complex<double> kappa,
                          size_t max_eval, double req_abs_error,
                          double req_rel_error, int error_norm,
                          int integration_type)
{
    _Complex double* zl_cast = (_Complex double*) zl;
    _Complex double* zt_cast = (_Complex double*) zt;
    _Complex double* s_cast = (_Complex double*) &s;
    _Complex double* gamma_cast = (_Complex double*) &gamma;
    _Complex double* kappa_cast = (_Complex double*) &kappa;
    return calculate_impedances(electrodes, num_electrodes, zl_cast, zt_cast,
                                *gamma_cast, *s_cast, mur, *kappa_cast, max_eval,
                                req_abs_error, req_rel_error, error_norm,
                                integration_type);
}

int Oimpedances_images(Electrode* electrodes, Electrode* images,
                       size_t num_electrodes, std::complex<double>* zl,
                       std::complex<double>* zt, std::complex<double> gamma,
                       std::complex<double> s, double mur,
                       std::complex<double> kappa, std::complex<double> ref_l,
                       std::complex<double> ref_t, size_t max_eval,
                       double req_abs_error, double req_rel_error,
                       int error_norm, int integration_type)
{
    _Complex double* zl_cast = (_Complex double*) zl;
    _Complex double* zt_cast = (_Complex double*) zt;
    _Complex double* s_cast = (_Complex double*) &s;
    _Complex double* gamma_cast = (_Complex double*) &gamma;
    _Complex double* kappa_cast = (_Complex double*) &kappa;
    _Complex double* ref_l_cast = (_Complex double*) &ref_l;
    _Complex double* ref_t_cast = (_Complex double*) &ref_t;
    /*impedances_images(electrodes, images, num_electrodes, zl_cast, zt_cast,
                      *gamma_cast, *s_cast, mur, *kappa_cast, *ref_l_cast,
                      *ref_t_cast, max_eval, req_abs_error, req_rel_error,
                      error_norm, integration_type);*/
    return 0;
}

int Oharmonic_impedance1(size_t ns, std::complex<double>* s,
                         std::complex<double>* kappa1,
                         std::complex<double>* kappa2,
                         std::complex<double>* gamma1,
                         Electrode* electrodes, Electrode* images,
                         size_t num_electrodes, double nodes[][3],
                         size_t num_nodes, size_t max_eval, double req_abs_error,
                         double req_rel_error, int error_norm, double rsource,
                         std::complex<double>* zh)
{
    _Complex double* s_cast = (_Complex double*) s;
    _Complex double* gamma1_cast = (_Complex double*) gamma1;
    _Complex double* kappa1_cast = (_Complex double*) kappa1;
    _Complex double* kappa2_cast = (_Complex double*) kappa2;
    _Complex double* zh_cast = (_Complex double*) zh;
    harmonic_impedance1(ns, s_cast, kappa1_cast, kappa2_cast, gamma1_cast,
                        electrodes, images, num_electrodes, nodes,
                        num_nodes, max_eval, req_abs_error, req_rel_error,
                        error_norm, rsource, zh_cast);
    return 0;
}

std::complex<double> sumall(std::complex<double>* z1, size_t n)
{
    /*
    std::complex<double>* z1_cast = reinterpret_cast< std::complex<double>* > (z1);
    std::complex<double> sum = mysumall(z1_cast, n);
    std::complex<double> * sum_cast;
    sum_cast = reinterpret_cast< std::complex<double> *> (&sum);
    return (*sum_cast);*
    return ((std::complex<double>) (foo((std::complex<double>*) (z1), n)));
    */
    std::complex<double> sum = std::complex<double> (0.0, 0.0);
    for (size_t i = 0; i < n; i++)
    {
        sum += z1[i];
    }
    return sum;
}
