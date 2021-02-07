#include "functions.h"
#include "WB.h"

unsigned seed;

int main(int argc, char *argv[])
{
	seed = 1234;
	double wtime = get_wall_time(); //timing

	/*------set simulation constants------*/
	const int N = atoi(argv[1]);
	const double coupling = atof(argv[2]);
	const double tFinal = atof(argv[3]);
	const double tTransision = atof(argv[4]);
	const double I_mean = atof(argv[5]);
	const double phi = atof(argv[6]);
	const double alpha = atof(argv[7]);
	const double beta = atof(argv[8]);
	const double E_syn = atof(argv[9]);
	const std::string NET_LABEL = argv[10];
	const int NUM_ENSEMBLE = atoi(argv[11]);
	const int PRINT_SPIKES = atoi(argv[12]);
	const int PRINT_VOLTAGES = atoi(argv[13]);
	const bool PRINT_SYNAPSE = atoi(argv[14]);
	const bool MEASURE_FREQUENCY = atoi(argv[15]);
	const double dt = 0.05;
	const double I_sigma = 0.0;
	int STEP_PRINT = 5;

	string RANDOMNESS;
	if (NUM_ENSEMBLE > 1)
		RANDOMNESS = "clock";
	else
		RANDOMNESS = "fix";

		/* ----------------------------set parameters--------------------------*/
		std::cout << NET_LABEL << "\t"
				  << coupling << "\t"
				  << I_mean << "\t"
				  << phi << "\n";


	dim2 Cij = read_matrix("dat/" + NET_LABEL + ".dat", N);
	dim2I adjList = adjmat_to_adjlist(Cij);
	
	for (int ens = 0; ens < NUM_ENSEMBLE; ens++)
	{
		WB wb(N, tFinal, dt);
		wb.set_initial_conditions(I_mean, I_sigma, RANDOMNESS);
		wb.set_params(tTransision, coupling,
					  phi, alpha, beta, E_syn,
					  NET_LABEL, adjList);
		wb.integrate(PRINT_SPIKES,
					 PRINT_VOLTAGES,
					 PRINT_SYNAPSE,
					 MEASURE_FREQUENCY,
					 STEP_PRINT,
					 ens);
	}

	wtime = get_wall_time() - wtime;
	// display_timing(wtime, 0); //print elapsed time

	return 0;
}
