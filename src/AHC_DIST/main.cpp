#include "AHC.h"

int main(int argc, char **argv)
{
	AHC ahc(argc, argv);
	ahc.performClustering();
	return 0;
}