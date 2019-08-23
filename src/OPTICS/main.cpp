/*
 * @brief The main driver function to call the OPTICS clustering algorithm on the data set
 * @author Lieyu Shi
 */


#include "OPTICS.h"

int main(int argc, char **argv)
{
	DensityClustering dclustering(argc, argv);
	dclustering.performClustering();
	return 0;
}
