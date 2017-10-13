#include "Initialization.h"

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

#pragma omp parallel for schedule(dynamic) num_threads(8)
	for (int i = 0; i < column; ++i)
	{
		diagonalRange(i,i) = diffRange(i);
	}
	clusterCenter = (clusterCenter+MatrixXf::Constant(Cluster,column,1.0))/2.0;

#pragma omp parallel for schedule(dynamic) num_threads(8)
	for (int i = 0; i < Cluster; ++i)
	{
		clusterCenter.row(i) = clusterCenter.row(i)*diagonalRange+range.row(1);
	}
}


void Initialization::generateFromSamples(MatrixXf& clusterCenter,
								    	 const int& column,
								    	 const MatrixXf& cArray,
								    	 const int& Cluster)
{
	clusterCenter = MatrixXf(Cluster,column);
	int number[Cluster];
	srand(time(0));
	const int& MaxNum = cArray.rows();
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
		chosen++;
	}
	assert(chosen==Cluster);
	assert(column==cArray.cols());

#pragma omp parallel for schedule(dynamic) num_threads(8)
	for (int i = 0; i < Cluster; ++i)
	{
		clusterCenter.row(i) = cArray.row(number[i]);
	}

}


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
				toCentroid = getDisimilarity(cArray, i, 
							number[j], normOption, object);
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

#pragma omp parallel for schedule(dynamic) num_threads(8)
	for (int i = 0; i < Cluster; ++i)
	{
		clusterCenter.row(i) = cArray.row(number[i]);
	}

}