#include "AHC.h"


/*
 * @brief The driver main function
 * @param argc: count of argument
 * @param argv: the char* array of argument
 * @return 0 if successful
 */
int main(int argc, char **argv)
{
	AHC ahc(argc, argv);
	ahc.performClustering();
	return 0;
}
