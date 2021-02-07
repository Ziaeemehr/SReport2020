#include <vector>
#include <iostream>
#include <algorithm>
#include <omp.h>

typedef std::vector<double> dim1;
typedef std::vector<dim1> dim2;
typedef std::vector<std::vector<int>> dim2I;

class WB
{
private:
    const int N;
    const double dt;
    double tfinal;
    int num_steps;

    double phi;
    double couplingOverN;
    double coupling;
    double E_syn;
    double alpha;
    double beta;
    double t_cut;
    double I_mean;
    std::string NET_LABEL;
    dim1 I_app;
    dim1 IC;
    dim2I adjlist;

    int THRESHOLD = -55;
    double Cm = 1.0;
    double E_L = -65.0;
    double E_Na = 55.0;
    double E_K = -90.0;
    double g_L = 0.1;
    double g_Na = 35.0;
    double g_K = 9.0;
    double g_syn = 0.1;

public:
    WB(int iN,
       double itfinal,
       double idt) : N(iN), dt(idt)
    {
        tfinal = itfinal;
        num_steps = int(tfinal / dt);
    }
    void set_params(
        const double,
        const double,
        const double,
        const double,
        const double,
        const double,
        const std::string,
        const dim2I &);

    virtual ~WB() {}

    void runge_kutta4_integrator(dim1 &);
    void set_initial_conditions(const double mu,
                                const double sigma,
                                const std::string);
    int integrate(const bool,
                  const bool,
                  const bool,
                  const bool,
                  const int,
                  const int);
    dim1 dydt(const dim1 &);

    

    // void print_coordinates(std::string, double);
    // void print_times(std::string );
    // dim1 get_times();
    // dim2 get_coordinates();
    // dim1 times;
    // dim2 coordinates;
};
