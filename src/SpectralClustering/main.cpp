/*
 * @brief The main function to use the spectral clustering on input integral curves
 * @author Lieyu Shi
 */


#include "SpectralClustering.h"


/*
 * @brief in case of running for many norm, would enable automatic parameter choice ahead of time
 *
 * @param[out] The Para class object to be revised
 */
void setPara(Para& p);


int main(int argc, char **argv)
{
	Para p;

	setPara(p);

	/* enable automatic option */
	bool automatic = true;

	SpectralClustering spectClus(argc, argv, p, automatic);

	spectClus.performClustering();

	return 0;
}


/*
 * @brief in case of running for many norm, would enable automatic parameter choice ahead of time
 *
 * @param[out] The Para class object to be revised
 */
void setPara(Para& p)
{
	/* 1.directly filling with last vertex; 2. uniform sampling, 3. equal-arc sampling */
	p.sampled = 2;

	/* Laplacian option: 1.Normalized Laplacian, 2.Unsymmetric Laplacian */
	p.LaplacianOption = 1;

	/* local scaling by sorted distance: true, false */
	p.isDistSorted = true;

	/* post-processing method: 1.k-means, 2.eigenvector rotation*/
	std::cout << "Input the post-processing: 1.k-means, 2.eigenvector rotation: " << std::endl;
	std::cin >> p.postProcessing;
	assert(p.postProcessing==1 || p.postProcessing==2);

	/* derivative method for eigen rotation: 1.numerical derivative, 2.true derivative */
	p.mMethod = 2;

	/* extraction option, 1. centroid, closest and furthest, 2. median, 3. statistical representation */
	p.extractOption = 1;
}
