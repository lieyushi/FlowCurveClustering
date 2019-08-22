/*
 * @brief The main driver function to perform DBSCAN clustering on the data set
 */


#include "DensityClustering.h"


/*
 * @brief The main function
 *
 * @param[in] argc count of argument
 * @param[in] argv char* array of argument
 * @return 0 if successful
 */
int main(int argc, char **argv)
{
	DensityClustering dclustering(argc, argv);
	dclustering.performClustering();
	return 0;
}
