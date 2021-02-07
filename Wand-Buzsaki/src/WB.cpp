#include "WB.h"
#include "functions.h"
#define M_PI 3.14159265358979323846

/*------------------------------------------------------------*/
void WB::set_initial_conditions(const double i_mean,
                                const double i_sigma,
                                const std::string RANDOMNESS)
{
    I_mean = i_mean;

    if (RANDOMNESS == "clock")
    {
        INITIALIZE_RANDOM_CLOCK();
    }
    else
    {
        INITIALIZE_RANDOM_F(seed);
    }

    IC.resize(4 * N);
    for (int i = 0; i < N; i++)
    {
        IC[i] = RANDOM * 20.0 - 70.0;
        IC[N + i] = RANDOM;
        IC[2 * N + i] = RANDOM;
        IC[3 * N + i] = RANDOM;
    }
    I_app.resize(N);
    for (int i = 0; i < N; i++)
        I_app[i] = RANDOM_GAUSS(i_sigma) + i_mean;

    FREE_RANDOM;
}
/*------------------------------------------------------------*/
void WB::set_params(
    const double tcut,
    const double coupling,
    const double phi,
    const double alpha,
    const double beta,
    const double E_syn,
    const std::string NET_LABEL,
    const dim2I &adjlist)

{
    this->t_cut = tcut;
    this->phi = phi;
    this->alpha = alpha;
    this->beta = beta;
    this->E_syn = E_syn;
    this->NET_LABEL = NET_LABEL;
    this->adjlist = adjlist;
    this->coupling = coupling;
    couplingOverN = coupling / (N + 0.0);
}
/*------------------------------------------------------------*/
int WB::integrate(
    const bool PRINT_SPIKES,
    const bool PRINT_VOLTAGES,
    const bool PRINT_SYNAPSES,
    const bool MEASURE_FREQUENCY,
    const int STEP_PRINT,
    const int ens)
{
    int index_transition = (int)(t_cut / dt);
    int m = int((num_steps - index_transition) / double(STEP_PRINT)) - 1;
    dim2 voltages(N, dim1(m));

    dim1 y = IC;
    dim1 ts(N);
    dim1 ts_new(N);
    dim1I aux(N);
    dim1I s(N);
    int n3 = 3 * N;

    const string subName = NET_LABEL + "-" +
                           to_string(coupling) + "-" +
                           to_string(phi) + "-" +
                           to_string(I_mean) + "-" +
                           to_string(beta) + ".txt";
    std::string spikeFileName = "../data/text/spk-" + subName;
    std::string measurementsFileName = "../data/text/par-" + subName;
    std::string voltageFileName = "../data/text/v-" + subName;
    std::string synapseFileName = "../data/text/syn-" + subName;

    FILE *measurmentsFile;
    FILE *synapseFile;
    measurmentsFile = fopen(measurementsFileName.c_str(), "a");
    if (PRINT_SYNAPSES && ens == 0)
        synapseFile = fopen(synapseFileName.c_str(), "w");

    dim2 spikes;
    spikes.resize(N);
    double t = 0.0;
    unsigned long int counter = 0;

    for (int it = 1; it < num_steps; ++it)
    {
        runge_kutta4_integrator(y);
        t = it * dt;

        for (int i = 0; i < N; i++)
        {
            if ((y[i] > THRESHOLD) && aux[i] == 0)
            {
                aux[i] = 1;
                if (t > t_cut)
                    spikes[i].push_back(t);
            }
            else if ((y[i] < THRESHOLD))
                aux[i] = 0;
        }

        if ((it % STEP_PRINT == 0) && (it > index_transition))
        {
            for (int i = 0; i < N; i++)
                voltages[i][counter] = y[i];
            counter++;
        }

        if (PRINT_SYNAPSES &&
            (it % STEP_PRINT == 0) &&
            (it > index_transition) &&
            (ens == 0))
        {
            for (int i = 0; i < N; i++)
                fprintf(synapseFile, "%18.6f", y[i + n3]);
            fprintf(synapseFile, "\n");
        }
    }

    // cout << "counter = " << counter << endl;
    // cout << "m = " << m << endl;

    double deltaPhase = 0.0;
    double spikeSynchrony = 0.0;
    double voltageSynchrony = 0.0;

    calculate_synchrony_measures(voltages,
                                 spikes,
                                 voltageSynchrony,
                                 spikeSynchrony,
                                 deltaPhase);

    double period = 0.0;
    if (MEASURE_FREQUENCY)
        period = measureAverageISI(spikes);

    fprintf(measurmentsFile, "%10d %18.10f %18.10f %18.10f %18.10f %18.10f \n",
            N, coupling, voltageSynchrony, spikeSynchrony, deltaPhase, 1.0 / period * 1000.0);
    fclose(measurmentsFile);

    // printf("%10d %18.10f %18.10f %18.10f %18.10f \n",
    //        N, coupling, vol_syn, spk_syn, delta_phi);

    if (PRINT_SPIKES && ens == 0)
        print_spikes_to_file(spikeFileName, spikes);

    if (PRINT_VOLTAGES && ens == 0) //wtite 2d vector to file
        print_spikes_to_file(voltageFileName, voltages);

    if (PRINT_SYNAPSES && ens == 0)
        fclose(synapseFile);
    return 0;
}
/*------------------------------------------------------------*/
void WB::runge_kutta4_integrator(dim1 &y)
{
    int n = y.size();
    dim1 k1(n), k2(n), k3(n), k4(n);
    dim1 f(n);

    k1 = dydt(y);

    for (int i = 0; i < n; i++)
        f[i] = y[i] + 0.5 * dt * k1[i];
    k2 = dydt(f);

    for (int i = 0; i < n; i++)
        f[i] = y[i] + 0.5 * dt * k2[i];
    k3 = dydt(f);

    for (int i = 0; i < n; i++)
        f[i] = y[i] + dt * k3[i];
    k4 = dydt(f);

    for (int i = 0; i < n; i++)
        y[i] += (k1[i] + 2.0 * (k2[i] + k3[i]) + k4[i]) * dt / 6.0;
}
/*------------------------------------------------------------*/
dim1 WB::dydt(const dim1 &x0)
{
    const int n2 = 2 * N;
    const int n3 = 3 * N;

    dim1 f(4 * N);

    for (int i = 0; i < N; i++)
    {
        double alpha_m = -0.1 * (x0[i] + 35.0) / (exp(-0.1 * (x0[i] + 35.0)) - 1.0);
        double alpha_h = 0.07 * exp(-(x0[i] + 58.0) / 20.0);
        double alpha_n = -0.01 * (x0[i] + 34.0) / (exp(-0.1 * (x0[i] + 34.0)) - 1.0);

        double beta_m = 4.0 * exp(-(x0[i] + 60.0) / 18.0);
        double beta_h = 1.0 / (exp(-0.1 * (x0[i] + 28.0)) + 1.0);
        double beta_n = 0.125 * exp(-(x0[i] + 44.0) / 80.0);

        double m = alpha_m / (alpha_m + beta_m);
        double F = 1.0 / (1.0 + exp(-0.5 * x0[i]));

        double I_Na = g_Na * m * m * m * x0[i + N] * (x0[i] - E_Na);
        double I_L = g_L * (x0[i] - E_L);
        double I_K = g_K * x0[i + n2] * x0[i + n2] * x0[i + n2] * x0[i + n2] * (x0[i] - E_K);

        double I_syn = 0.0;
        int counter = 0;
        for (int j = 0; j < adjlist[i].size(); j++)
        {
            int k = adjlist[i][j];
            I_syn += (g_syn * x0[k + n3] * (x0[k] - E_syn));
        }

        I_syn = couplingOverN * I_syn;

        f[i] = -I_Na - I_K - I_L - I_syn + I_app[i];                             // dv/dt
        f[i + N] = phi * (alpha_h * (1 - x0[i + N]) - beta_h * x0[i + N]);       // dh/dt
        f[i + 2 * N] = phi * (alpha_n * (1 - x0[i + n2]) - beta_n * x0[i + n2]); // dn/dt
        f[i + 3 * N] = alpha * F * (1 - x0[i + n3]) - beta * x0[i + n3];         // ds/dt
    }

    return f;
}
/*------------------------------------------------------------*/
