/*
 * @brief This is the source cpp for the class Initialization.h. It is for the k-means initialization
 * @author Lieyu Shi
 */


#include "Initialization.h"


/*
 * @brief To generate the random coordinates for the k-means initialization
 *
 * @param[out] clusterCenter The random initialization to be updated
 * @param[in] column The column size
 * @param[in] cArray The input matrix coordinates
 * @param[in] Cluster The number of centroids
 */
void Initialization::generateRandomPos(MatrixXf& clusterCenter,
								  	   const int& column,
								       const MatrixXf& cArray,
								       const int& Cluster)
{
	clusterCenter = MatrixXf::Random(Cluster, column);
	MatrixXf range(2, column);
	range.row(0) = cArray.colwise().maxCoeff();  //first row contains max
	range.row(1) = cArray.colwise().minCoeff();  //second row contains min
	VectorXf diffRange = range.row(0)-range.row(1);

	MatrixXf diagonalRange = MatrixXf::Zero(column,column);

#pragma omp parallel for schedule(static) num_threads(8)
	for (int i = 0; i < column; ++i)
	{
		diagonalRange(i,i) = diffRange(i);
	}
	clusterCenter = (clusterCenter+MatrixXf::Constant(Cluster,column,1.0))/2.0;

#pragma omp parallel for schedule(static) num_threads(8)
	for (int i = 0; i < Cluster; ++i)
	{
		clusterCenter.row(i) = clusterCenter.row(i)*diagonalRange+range.row(1);
	}
}


/*
 * @brief To generate the initialization from the samples
 *
 * @param[out] clusterCenter The initialized centroid coordinates
 * @param[in] column The size of column
 * @param[in] cArray The original matrix coordinates as input
 * @paramp[in] Cluster The count of clusters
 */
void Initialization::generateFromSamples(MatrixXf& clusterCenter,
								    	 const int& column,
								    	 const MatrixXf& cArray,
								    	 const int& Cluster)
{
	clusterCenter = MatrixXf(Cluster,column);
	std::vector<int> number(Cluster);
	srand(time(0));

	const int& MaxNum = cArray.rows();

	std::cout << MaxNum << std::endl;

	number[0] = rand()%MaxNum;
	int randNum, chosen = 1;
	bool found;
	for (int i = 1; i < Cluster; ++i)
	{
		do
		{
			randNum = rand()%MaxNum;
			found = false;
			for(int j=0;j<chosen;j++)
			{
				if(randNum==number[j])
				{
					found = true;
					break;
				}
			}
		}while(found!=false);
		number[i] = randNum;
		++chosen;
	}
	assert(chosen==Cluster);
	assert(column==cArray.cols());

#pragma omp parallel for schedule(static) num_threads(8)
	for (int i = 0; i < Cluster; ++i)
	{
		clusterCenter.row(i) = cArray.row(number[i]);
	}
}


/*
 * @brief This is the k-means++ initialization
 *
 * @param[out] clusterCenter The cluster centroid to be initialized
 * @param[in] column The size of column
 * @param[in] cArray The matrix coordinates of the lines
 * @param[in] normOption The norm option
 * @param[in] object The MetricPreparation
 */
void Initialization::generateFarSamples(MatrixXf& clusterCenter,
								   	    const int& column,
								   		const MatrixXf& cArray,
								   		const int& Cluster,
								   		const int& normOption,
								   		const MetricPreparation& object)
{
	assert(column==cArray.cols());
	const int Total = cArray.rows();
	clusterCenter = MatrixXf(Cluster,column);
	int number[Cluster], selection;
	srand(time(0));
	const int& MaxNum = cArray.rows();
	number[0] = rand()%MaxNum;
	int chosen = 1;

	float percentage, nearest, toCentroid;
	VectorXf distance(Total);
	double squredSummation;
	float left, right;
	while(chosen<Cluster)
	{
		percentage = float(rand()/(double)RAND_MAX);
		for (int i = 0; i < Total; ++i)
		{
			nearest = FLT_MAX;
			for (int j = 0; j < chosen; ++j)
			{
				toCentroid = getDisimilarity(cArray, i, number[j], normOption, object);
				if(nearest>toCentroid)
					nearest=toCentroid;
			}
			distance(i)=nearest*nearest;
		}
		squredSummation = distance.sum();
		left = 0.0, right = 0.0;
		for (int i = 0; i < Total; ++i)
		{
			left = right;
			right += float((double)distance(i)/squredSummation);
			if(left < percentage && percentage <= right)
			{
				selection = i;
				break;
			}
		}
		number[chosen] = selection;
		chosen++;
	}

#pragma omp parallel for schedule(static) num_threads(8)
	for (int i = 0; i < Cluster; ++i)
	{
		clusterCenter.row(i) = cArray.row(number[i]);
	}

}
