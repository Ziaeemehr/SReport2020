#include "Traub.h"
#include "lib.h"

/*------------------------------------------------------------*/
void Traub::setParams(const int N,
                      const double dt,
                      const double tTransition,
                      const double tSimulation,
                      const double coupling,
                      const double gm,
                      const double gahp,
                      const double gSynaptic,
                      const double IApp,
                      const std::string netLabel,
                      const dim2I &adjList)
{
    this->N = N;
    this->dt = dt;
    this->gm = gm;
    this->gahp = gahp;
    this->IApp = IApp;
    this->adjList = adjList;
    this->netLabel = netLabel;
    this->coupling = coupling;
    this->gSynaptic = gSynaptic;
    this->tTransition = tTransition;
    this->tSimulation = tSimulation;
}
/*------------------------------------------------------------*/
dim1 Traub::setInitialConditions(const std::string randomness)
{
    if (randomness == "clock")
        {
            INITIALIZE_RANDOM_CLOCK();
        }
    else
        {
            INITIALIZE_RANDOM_F(seed);
        }
    dim1 initialConditions(7 * N);

    for (int i = 0; i < N; i++)
    {
        initialConditions[i] = RANDOM * 20.0 - 95.0;
        initialConditions[N + i] = RANDOM;
        initialConditions[2 * N + i] = RANDOM;
        initialConditions[3 * N + i] = RANDOM;
        initialConditions[4 * N + i] = RANDOM;
        initialConditions[5 * N + i] = RANDOM;
        initialConditions[6 * N + i] = RANDOM;
    }

    FREE_RANDOM;

    return initialConditions;
}
/*------------------------------------------------------------*/
int Traub::integrate(odeSys ode_system,
                     const bool PRINT_SPIKES,
                     const bool PRINT_VOLTAGE,
                     const bool PRINT_SYNAPSES,
                     const int STEPV,
                     const int ENS,
                     const std::string seedState)
{
    int indexTransition = (int)(tTransition / dt);
    long unsigned numSteps = (int)(tSimulation / dt);

    int m = int((numSteps - indexTransition) / (STEPV + 0.0)) - 1;
    dim2 voltages(N, dim1(m));
    dim2 synapses(N, dim1(m));
    dim1 y = setInitialConditions(seedState);
    dim1I auxilaryVec(N);

    // --------------------------------------------------------------
    std::string spikesFileName = "../data/text/spk-" +
                          netLabel + "-" +
                          to_string(coupling) + "-" +
                          to_string(IApp) + ".txt";
    std::string parname = "../data/text/par-" +
                          netLabel + "-" +
                          to_string(coupling) + "-" +
                          to_string(IApp) + ".txt";
    std::string voltagesFileName = "../data/text/v-" +
                           netLabel + "-" +
                           to_string(coupling) + "-" +
                           to_string(IApp) + ".txt";
    std::string synapsesFileName = "../data/text/s-" +
                           netLabel + "-" +
                           to_string(coupling) + "-" +
                           to_string(IApp) + ".txt";

    FILE *outputFile;
    outputFile = fopen(parname.c_str(), "a");
    // --------------------------------------------------------------

    dim2 spikes;
    spikes.resize(N);
    double t = 0.0;
    unsigned long int counter = 0;

    for (long unsigned it = 1; it < numSteps; ++it)
    {
        t = it * dt;
        runge_kutta4_integrator(y, ode_system);

        if (t > tTransition)
        {
            for (int i = 0; i < N; i++)
            {
                if ((y[i] > spikeThreshold) && auxilaryVec[i] == 0)
                {
                    auxilaryVec[i] = 1;
                    spikes[i].push_back(t);
                }
                else if ((y[i] < spikeThreshold))
                    auxilaryVec[i] = 0;
            }
        }

        if ((it % STEPV == 0) && (it > indexTransition))
        {
            for (int i = 0; i < N; i++)
            {
                voltages[i][counter] = y[i];
                synapses[i][counter] = y[i + (6 * N)];
            }
            counter++;
        }
    }

    double voltageSynchrony = 0.0, spikeSynchrony = 0.0, deltaPhi = 0.0;
    calculate_synchrony_measures(voltages, spikes, voltageSynchrony, spikeSynchrony, deltaPhi);


    fprintf(outputFile, "%10d %18.10f %18.10f %18.10f %18.10f \n",
            N, coupling, voltageSynchrony, spikeSynchrony, deltaPhi);
    fclose(outputFile);

    // printf("%10d %18.10f %18.10f %18.10f %18.10f \n",
    //        N, coupling, voltageSynchrony, spikeSynchrony, deltaPhi);

    if (PRINT_SPIKES && ENS == 0)
        print_spikes_to_file(spikesFileName, spikes);

    if (PRINT_VOLTAGE && ENS == 0) //wtite 2d vector to file
        print_spikes_to_file(voltagesFileName, voltages);
    
    if (PRINT_SYNAPSES && ENS == 0)
        print_spikes_to_file(synapsesFileName, synapses);

    return 0;
}
// /*------------------------------------------------------------*/
void Traub::runge_kutta4_integrator(dim1 &y, odeSys dydt)
{
    int n = y.size();
    dim1 k1(n), k2(n), k3(n), k4(n);
    dim1 f(n);

    k1 = (this->*dydt)(y);

    for (int i = 0; i < n; i++)
        f[i] = y[i] + 0.5 * dt * k1[i];
    k2 = (this->*dydt)(f);

    for (int i = 0; i < n; i++)
        f[i] = y[i] + 0.5 * dt * k2[i];
    k3 = (this->*dydt)(f);

    for (int i = 0; i < n; i++)
        f[i] = y[i] + dt * k3[i];
    k4 = (this->*dydt)(f);

    for (int i = 0; i < n; i++)
        y[i] += (k1[i] + 2.0 * (k2[i] + k3[i]) + k4[i]) * dt / 6.0;
}
/*------------------------------------------------------------*/
dim1 Traub::traubModel(const dim1 &x)
{

    dim1 v(x.begin(), x.begin() + N);
    dim1 m(x.begin() + N, x.begin() + (2 * N));
    dim1 n(x.begin() + (2 * N), x.begin() + (3 * N));
    dim1 h(x.begin() + (3 * N), x.begin() + (4 * N));
    dim1 w(x.begin() + (4 * N), x.begin() + (5 * N));
    dim1 ca(x.begin() + (5 * N), x.begin() + (6 * N));
    dim1 s(x.begin() + (6 * N), x.end());

    dim1 f(7 * N);
    int n2 = 2 * N;
    int n3 = 3 * N;
    int n4 = 4 * N;
    int n5 = 5 * N;
    int n6 = 6 * N;

    for (size_t i = 0; i < N; ++i)
    {
        double mlInf = 1.0 / (1.0 + exp(-(v[i] + 25.0) / 2.5));
        double am = 0.32 * (v[i] + 54.0) / (1.0 - exp(-(v[i] + 54.0) / 4.0));
        double bm = 0.28 * (v[i] + 27.0) / (exp((v[i] + 27.0) / 5.0) - 1.0);
        double ah = 0.128 * exp(-(v[i] + 50.0) / 18.0);
        double bh = 4.0 / (1.0 + exp(-(v[i] + 27.0) / 5.0));
        double an = 0.032 * (v[i] + 52.0) / (1.0 - exp(-(v[i] + 52.0) / 5.0));
        double bn = 0.5 * exp(-(v[i] + 57.0) / 40.0);
        double tw = 100.0 / (3.3 * exp((v[i] + 35.0) / 20.0) + exp(-(v[i] + 35.0) / 20.0));
        double wInf = 1.0 / (1.0 + exp(-(v[i] + 35.0) / 10.0));
        double ica = mlInf * (v[i] - 120.0);
        
        double ISyn = 0.0;
        for (size_t j = 0; j < adjList[i].size(); j++)
        {
            int k = adjList[i][j];
            ISyn += s[k] * vReset;
        }
        ISyn *= (-1.0 * coupling * gSynaptic);
        f[i] = IApp - 100.0 * h[i] * pow(m[i], 3.0) * (v[i] - 50.0) - 
                (80.0 * pow(n[i], 4.0) + gm * w[i] + gahp * (ca[i] / (ca[i] + 1.0))) * 
                (v[i] + 100) - 0.2 * (v[i] + 67.0) - ica + ISyn;
        f[i+N] = am * (1.0 - m[i]) - bm * m[i];
        f[i+n2] = an * (1.0 - n[i]) - bn * n[i];
        f[i+n3] = ah * (1.0 - h[i]) - bh * h[i];
        f[i+n4] = (wInf - w[i]) / tw;
        f[i+n5] = -0.002 * ica - ca[i] / 80.0;
        f[i+n6] = 2.0 * (1.0 - s[i]) / (1.0 + exp(-(v[i] + 10.0) / 10.0)) - 0.1 * s[i];
    }

    return f;
}
/*------------------------------------------------------------*/
