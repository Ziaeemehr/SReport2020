#include "lib.h"
#include "Traub.h"

unsigned seed;

int main(int argc, char *argv[])
{
	seed = 1234;
	double wtime = get_wall_time(); //timing

	/*------set simulation constants------*/
	const int N = atoi(argv[1]);
	const double tSimulation = atof(argv[2]);
	const double tTransition = atof(argv[3]);
	const double coupling = atof(argv[4]);
	const double gm = atof(argv[5]);
	const double gahp = atof(argv[6]);
	const double gSynaptic = atof(argv[7]);
	const double IApp = atof(argv[8]);
	const std::string netLabel = argv[9];
	const int numEnsemble = atoi(argv[10]);
	const int PRINT_SPIKES = atoi(argv[11]);
	const int PRINT_VOLTAGE = atoi(argv[12]);
	const int PRINT_SYNAPSE = atoi(argv[13]);
	const double dt = 0.02;

	/* ----------------------------set parameters--------------------------*/
	printf("network = %s, coupling=%10.6f, I=%10.4f\n",
		   netLabel.c_str(),
		   coupling, IApp);

	
	dim2I adjMat = readMatrix("dat/" + netLabel + ".dat", N);
	dim2I adjList = adjMatToAdjList(adjMat);

	for (int ens = 0; ens < numEnsemble; ens++)
	{
		Traub traub;
		traub.setParams(N,
						dt,
						tTransition,
						tSimulation,
						coupling,
						gm,
						gahp,
						gSynaptic,
						IApp,
						netLabel,
						adjList);

		traub.integrate(&Traub::traubModel,
						PRINT_SPIKES,
						PRINT_VOLTAGE,
						PRINT_SYNAPSE, 
						1,  // stepV
						ens, // ensemble
						"fix"
		);
	}

	wtime = get_wall_time() - wtime;
	// display_timing(wtime, 0);

	return 0;
}
