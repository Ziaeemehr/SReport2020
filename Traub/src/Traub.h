#include <vector>
#include <iostream>
#include <algorithm>
#include <omp.h>

class Traub;

typedef std::vector<double> dim1;
typedef std::vector<dim1> dim2;
typedef std::vector<std::vector<int>> dim2I;
typedef dim1 (Traub::*odeSys)(const dim1 &);

class Traub
{
private:
    int N;
    double dt;
    double tSimulation;
    double tTransition;
    double coupling;
    double gm;
    double gahp;
    double IApp;
    std::string netLabel;
    dim2I adjList;
    double gSynaptic;
    double vReset = -95.0;
    double spikeThreshold = -55.0;

public:
    void
    setParams(const int N,
              const double dt,
              const double tTransition,
              const double tSimulation,
              const double coupling,
              const double gm,
              const double gahp,
              const double gSynaptic,
              const double IApp,
              const std::string netLable,
              const dim2I &adjList);

    virtual ~Traub() {}

    dim1 setInitialConditions(const std::string randomness);

    void runge_kutta4_integrator(dim1 &, odeSys);
    int integrate(odeSys,
                  const bool,
                  const bool,
                  const bool,
                  const int,
                  const int,
                  const std::string);
    dim1 traubModel(const dim1 &x0);

    // void print_coordinates(std::string, double);
    // void print_times(std::string );
    // dim1 get_times();
    // dim2 get_coordinates();
    // dim1 times;
    // dim2 coordinates;
};
