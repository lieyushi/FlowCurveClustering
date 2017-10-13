#include "DensityClustering.h"

int main(int argc, char **argv)
{
	DensityClustering dclustering(argc, argv);
	dclustering.performClustering();
	return 0;
}