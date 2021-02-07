#include "functions.h"

extern unsigned seed;

/*------------------------------------------------------------*/
double calculate_voltage_synchrony(const dim2 &voltages)
{
    int n = voltages.size();
    int nstep = voltages[0].size();
    dim1 vg(nstep);

    // global volatage at each time step
    for (int j = 0; j < nstep; j++)
    {
        for (int i = 0; i < n; i++)
            vg[j] += voltages[i][j];
        vg[j] /= (n + 0.0);
    }

    double vg_m = mean(vg, 0);

    double O = 0.0;
    for (int j = 0; j < nstep; j++)
    {
        double tmp = (vg[j] - vg_m);
        O += (tmp * tmp);
    }
    O /= (nstep + 0.0);

    dim1 v(nstep);
    double denom = 0.0;

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < nstep; j++)
            v[j] = voltages[i][j];
        double m = mean(v, 0);

        double sigma = 0;
        for (int j = 0; j < nstep; j++)
        {
            double tmp = (v[j] - m);
            sigma += (tmp * tmp);
        }
        sigma /= (nstep + 0.0);
        denom += sigma;
    }
    denom /= (n + 0.0);

    double xi = O / (denom + 0.0);
    return xi;
}
/*------------------------------------------------------------*/
void calculate_synchrony_measures(const dim2 &voltages,
                                  const dim2 &spiketrains,
                                  double &vol_syn,
                                  double &spk_syn,
                                  double &delta_phi)
{
    int n = voltages.size();

    dim1 flatten_vec = flatten(spiketrains);
    delta_phi = calculate_phase_synchrony(voltages);

    if (flatten_vec.size() < 2)
    {
        spk_syn = 0.0;
        vol_syn = 0.0;
    }
    else
    {
        spk_syn = calculate_spike_synchrony(flatten_vec, n);
        vol_syn = calculate_voltage_synchrony(voltages);
    }

    // std::cout << "size " <<flatten_vec.size()<< "\n";
    // for (int i=0; i<n; i++)
    //     cout << spiketrains[i].size()<< " ";
    // cout << "\n";
}
/*------------------------------------------------------------*/
double calculate_spike_synchrony(dim1 &spiketrains, const int nodes)
{
    int n = spiketrains.size();
    std::sort(spiketrains.begin(), spiketrains.end());

    // for(int i=0; i<20; i++)
    //     cout << spiketrains[i] << " ";
    // cout << "\n";

    auto tau = calculate_isi(spiketrains);

    dim1 tau2(tau.size());
    for (int i = 0; i < tau.size(); i++)
        tau2[i] = tau[i] * tau[i];

    double tau_m = mean(tau, 0);
    double tau2_m = mean(tau2, 0);
    double tau_m2 = tau_m * tau_m;
    double burst = ((sqrt(tau2_m - tau_m2)) / (tau_m + 0.0) - 1.0) / (sqrt(nodes) + 0.0);

    return burst;
}
/*------------------------------------------------------------*/
std::vector<double> calculate_isi(const std::vector<double> &v)
{
    int n = v.size();
    std::vector<double> v_isi(n - 1);
    for (int i = 0; i < (n - 1); i++)
        v_isi[i] = v[i + 1] - v[i];
    return v_isi;
}

/*------------------------------------------------------------*/
template <typename T>
std::vector<T> flatten(const std::vector<std::vector<T>> &v)
{
    std::size_t total_size = 0;
    for (const auto &sub : v)
        total_size += sub.size();
    std::vector<T> result;
    result.reserve(total_size);
    for (const auto &sub : v)
        result.insert(result.end(), sub.begin(), sub.end());
    return result;
}
/*------------------------------------------------------------*/
dim2 read_matrix(std::string filename, int Node)
{
    /*get filename and number of row to read a square matrix 
    intput:
        filename: name of text file to read
        Node: number of rows or cols of square matrix
    return:
        matrix as a 2 dimensional double vector

    example: read_matrix("A.txt", 100);
    */
    std::ifstream ifile(filename);

    /*to check if input file exists*/
    if (fileExists(filename))
    {
        dim2 Cij(Node, dim1(Node));

        for (int i = 0; i < Node; i++)
        {
            for (int j = 0; j < Node; j++)
            {
                ifile >> Cij[i][j];
            }
        }
        ifile.close();
        return Cij;
    }
    else
    {
        std::cerr << "\n file : " << filename << " not found \n";
        exit(2);
    }
}
/*------------------------------------------------------------*/
bool fileExists(const std::string &filename)
{
    /*return true if input file exists*/
    struct stat buf;
    if (stat(filename.c_str(), &buf) != -1)
    {
        return true;
    }
    return false;
}
/*------------------------------------------------------------*/
dim1 normal_dist(int N, double mu, double sigma)
{
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(mu, sigma);
    dim1 r(N);
    for (int i = 0; i < N; ++i)
        r[i] = distribution(generator);
    return r;
}
/*------------------------------------------------------------*/
double get_wall_time()
{
    /*measure real passed time */
    struct timeval time;
    if (gettimeofday(&time, NULL))
    {
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
/*------------------------------------------------------------*/
double get_cpu_time()
{
    /*measure cpu passed time*/
    return (double)clock() / CLOCKS_PER_SEC;
}
/*------------------------------------------------------------*/
void display_timing(double wtime, double cptime)
{
    int wh, ch;
    int wmin, cpmin;
    double wsec, csec;
    wh = (int)wtime / 3600;
    ch = (int)cptime / 3600;
    wmin = ((int)wtime % 3600) / 60;
    cpmin = ((int)cptime % 3600) / 60;
    wsec = wtime - (3600. * wh + 60. * wmin);
    csec = cptime - (3600. * ch + 60. * cpmin);
    printf("Wall Time : %d hours and %d minutes and %.4f seconds.\n", wh, wmin, wsec);
    // printf ("CPU  Time : %d hours and %d minutes and %.4f seconds.\n",ch,cpmin,csec);
}
/*------------------------------------------------------------*/
void print_spikes_to_file(const std::string spkname, const dim2 &spikes)
{
    FILE *spkfile;
    spkfile = fopen(spkname.c_str(), "w");

    for (int i = 0; i < spikes.size(); i++)
    {
        for (int j = 0; j < spikes[i].size(); j++)
        {
            fprintf(spkfile, "%18.6f", spikes[i][j]);
        }
        fprintf(spkfile, "\n");
    }

    fclose(spkfile);
}
/*------------------------------------------------------------*/
void hilbert_transform(const dim1 &in, fftw_complex *out)
{
    // copy data into the complex array
    int N = in.size();
    for (int i = 0; i < N; i++)
    {
        out[i][REAL] = in[i];
        out[i][IMAG] = 0;
    }
    // create a DFT plan and execute it
    fftw_plan plan = fftw_plan_dft_1d(N, out, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    int hN = N >> 1; // half of the length(N/2);
    int numRem = hN; // the number of remaining elements

    // multiply the appropriate values by 2
    // those that sould be multiplied by 1 are left intact because they wouldn't change
    for (int i = 1; i < hN; i++)
    {
        out[i][REAL] *= 2;
        out[i][IMAG] *= 2;
    }

    // if the length is even, the number of remaining elements decreases by 1
    if (N % 2 == 0)
        numRem--;

    // if it is odd and greater than 1, the middle value must be multiplied by 2
    else if (N > 1)
    {
        out[hN][REAL] *= 2;
        out[hN][IMAG] *= 2;
    }
    // set the remaining values to 0
    // multiply by 0 gives 0 so we don't care about multiplication
    memset(&out[hN + 1][REAL], 0, numRem * sizeof(fftw_complex));

    // create an IDFT plan and execute it
    plan = fftw_plan_dft_1d(N, out, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    // dom some cleaning
    fftw_destroy_plan(plan);
    fftw_cleanup();

    // scale the IDFT output
    for (int i = 0; i < N; i++)
    {
        out[i][REAL] /= N;
        out[i][IMAG] /= N;
    }
}
/*------------------------------------------------------------*/
dim1 unwrap(const dim1 &in)
{
    int len = in.size();
    dim1 out(len);
    out[0] = in[0];
    for (int i = 1; i < len; i++)
    {
        double d = in[i] - in[i - 1];
        d = d > M_PI ? d - 2 * M_PI : (d < -M_PI ? d + 2 * M_PI : d);
        out[i] = out[i - 1] + d;
    }
    return out;
}
/*------------------------------------------------------------*/
dim1 angle_arr(const fftw_complex *z, const int n)
{
    dim1 angle_z(n);
    for (int i = 0; i < n; i++)
    {
        std::complex<double> z0(z[i][REAL], z[i][IMAG]);
        angle_z[i] = std::arg(z0);
    }

    return angle_z;
}
/*------------------------------------------------------------*/
double calculate_phase_synchrony(const dim2 &voltages)
{
    int n = voltages.size();
    int nstep = voltages[0].size();
    dim1 delta_phi(nstep);
    dim2 inst_phase(n, dim1(nstep));

    for (int i = 0; i < n; i++)
    {
        fftw_complex z[nstep];
        hilbert_transform(voltages[i], z);
        dim1 angle_z = angle_arr(z, nstep);
        // inst_phase[i] = unwrap(angle_z);
        for (int j = 0; j < nstep; j++)
            inst_phase[i][j] = fmod(angle_z[j], (2.0 * M_PI));
    }

    for (int j = 0; j < nstep; j++)
    {
        dim1 x(n);
        for (int i = 0; i < n; i++)
            x[i] = inst_phase[i][j];
        delta_phi[j] = calculate_phase_difference(x);
    }
    double result = mean(delta_phi, 0);

    return result;
}
/*------------------------------------------------------------*/
double calculate_phase_difference(const dim1 &x)
{
    int n = x.size();
    double sum = 0.0;
    for (int i = 1; i < n; i++)
    {
        sum += abs(x[i] - x[0]);
    }
    double result = sum / (n - 1.0);
    return result;
}
/*------------------------------------------------------------*/
dim2I adjmat_to_adjlist(const dim2 &A)
{
    int row = A.size();
    int col = A[0].size();
    dim2I adjlist;
    adjlist.resize(row);

    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            if (A[i][j] > 1.e-8)
                adjlist[i].push_back(j);
        }
    }

    return adjlist;
}
/*------------------------------------------------------------*/

double measureAverageISI(const dim2& spikes)
{
    int counter = 0;
    double period = 0.0;
    int n = spikes.size();
    for (int i = 0; i < n; i++)
    {
        int m = spikes[i].size();
        if (m > 3)
        {
            period += (spikes[i][m - 1] - spikes[i][m - 2]);
            counter++;
        }
    }
    if (counter > 0)
    {
        period /= counter;
        return period;
    }
    else
        return -1.0;
}