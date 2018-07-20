/*
 * main.cpp
 *
 *  Created on: Mar 13, 2018
 *      Author: lieyu
 */

#include "ReadClustering.h"

int main(int argc, char* argv[])
{
	if(argc!=2)
	{
		std::cout << "Error for argument input!" << std::endl;
		exit(1);
	}

	/* create ReadClustering object and get evaluation */

	ReadClustering rc;

	rc.getEvaluation(argv[1]);

	return 0;
}

