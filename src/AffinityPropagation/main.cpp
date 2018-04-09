#include "AffinityPropagation.h"

/* in case of running for many norm, would enable automatic parameter choice */
void setPara(Para& p);


int main(int argc, char **argv)
{
	Para p;

	setPara(p);

	/* enable automatic option */
	bool automatic = true;

	AffinityPropagation ap(argc, argv, p, automatic);

	ap.performClustering();

	return 0;
}



void setPara(Para& p)
{
	/* 1.directly filling with last vertex; 2. uniform sampling, 3. equal-arc sampling */
	p.sampled = 2;

	/* extraction option, 1. centroid, closest and furthest, 2. median, 3. statistical representation */
	p.extractOption = 1;

	/* max iteration for AP clustering */
	p.maxIteration = 20;

}
