/*
 * @brief The main driver function to perform AHC clustering on the data set provided from the argument
 * @author Lieyu Shi
 */


#include "AHC.h"


/*
 * @brief The driver main function
 *
 * @param[in] argc count of argument
 * @param[in] argv the char* array of argument
 * @return 0 if successful
 */
int main(int argc, char **argv)
{
	AHC ahc(argc, argv);
	ahc.performClustering();
	return 0;
}
