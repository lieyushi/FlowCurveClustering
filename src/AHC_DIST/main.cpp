/*
 * @brief The main driver function to perform the AHC with dist clustering
 * @author Lieyu Shi
 */

#include "AHC.h"


/*
 * @brief The main function to read and perform clustering on the data set from the argument
 *
 * @param[in] argc The count of argument
 * @param[in] argv The argument strings
 * @return 0 if it is successful
 */
int main(int argc, char **argv)
{
	AHC ahc(argc, argv);
	ahc.performClustering();
	return 0;
}
