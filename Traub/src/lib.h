#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <time.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <complex>
#include <random>
#include <cmath>
#include <chrono>
#include <sstream>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <assert.h>
#include <algorithm>
#include <functional>
#include <iomanip>
#include <omp.h>
#include <fftw3.h>
#include <string.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>

#define RANDOM gsl_rng_uniform(gsl_rng_r)
#define RANDOM_INT(A) gsl_rng_uniform_int(gsl_rng_r, A)
#define RANDOM_GAUSS(S) gsl_ran_gaussian(gsl_rng_r, S)
#define RANDOM_POISSON(M) gsl_ran_poisson(gsl_rng_r, M)
#define INITIALIZE_RANDOM_CLOCK()                                                                        \
    {                                                                                                    \
        gsl_rng_env_setup();                                                                             \
        if (!getenv("GSL_RNG_SEED"))                                                                     \
            gsl_rng_default_seed = std::chrono::high_resolution_clock::now().time_since_epoch().count(); \
        gsl_rng_T = gsl_rng_default;                                                                     \
        gsl_rng_r = gsl_rng_alloc(gsl_rng_T);                                                            \
    }
#define INITIALIZE_RANDOM_F(seed)             \
    {                                         \
        gsl_rng_env_setup();                  \
        if (!getenv("GSL_RNG_SEED"))          \
            gsl_rng_default_seed = seed;      \
        gsl_rng_T = gsl_rng_default;          \
        gsl_rng_r = gsl_rng_alloc(gsl_rng_T); \
    }
#define FREE_RANDOM gsl_rng_free(gsl_rng_r)

static const gsl_rng_type *gsl_rng_T;
static gsl_rng *gsl_rng_r;

#define REAL 0
#define IMAG 1

using dim1I = std::vector<int>;
using dim1 = std::vector<double>;
using dim2 = std::vector<std::vector<double>>;
using dim2I = std::vector<std::vector<int>>;

using std::cout;
using std::endl;
using std::string;
using std::to_string;
using std::vector;

extern unsigned seed;

dim2I readMatrix(std::string, int);
bool fileExists(const std::string &);
void display_timing(double wtime, double cptime);
double get_wall_time();
double get_cpu_time();
dim1 calculate_coupling(const dim2 &, const double);

void calculate_synchrony_measures(const dim2 &,
                                  const dim2 &,
                                  double &,
                                  double &,
                                  double &);
double calculate_voltage_synchrony(const dim2 &);
double calculate_spike_synchrony(dim1 &, const int);
std::vector<double> calculate_isi(const std::vector<double> &v);

template <typename T>
std::vector<T> flatten(const std::vector<std::vector<T>> &v);

void print_spikes_to_file(const std::string, const dim2 &);
double calculate_phase_difference(const dim1 &x);
double calculate_phase_synchrony(const dim2 &voltages);
dim1 angle_arr(const fftw_complex *z, const int n);
void hilbert_transform(const dim1 &in, fftw_complex *out);
dim2I adjMatToAdjList(const dim2I &A);
dim1 unwrap(const dim1 &in);

inline double mean(const vector<double> &v, const int id)
{ /*average the vector from element "id" to end of the vector */
    assert(v.size() > id);
    return accumulate(v.begin() + id, v.end(), 0.0) / (v.size() - id);
}

#endif