#include "SpectralClustering.h"

int main(int argc, char **argv)
{
	SpectralClustering spectClus(argc, argv);
	spectClus.performClustering();
	return 0;
}
