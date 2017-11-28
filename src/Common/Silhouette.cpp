#include "Silhouette.h"

Silhouette::Silhouette()
{
	sData = std::vector<float>();
	sCluster = std::vector<float>();
}


Silhouette::~Silhouette()
{
	reset();
}


void Silhouette::reset()
{
	sData.clear();
	sCluster.clear();
	sAverage = 0;
}

void Silhouette::computeValue(const int& normOption,
						   	  const MatrixXf& array, 
						   	  const int& Row, 
						      const int& Column,
						   	  const std::vector<int>& group, 
						   	  const MetricPreparation& object,
						   	  const int& groupNumber)
{
	sData.clear();
	sCluster.clear();
	sData = std::vector<float>(Row,0);
	assert(Row==group.size());

	//groupNumber doesn't include noise group
	sCluster = std::vector<float>(groupNumber, 0);
	std::vector<std::vector<int> > storage(groupNumber, std::vector<int>());

	//whether some group marked as -1 noise or not
	int noise = 0;
	for (int i = 0; i < group.size(); ++i)
	{
		if(group[i]<0)
		{
			++noise;
			continue;
		}
		else
			storage[group[i]].push_back(i);
	}

	// compute Silhouette value for each data
	float sSummation = 0;

#pragma omp parallel num_threads(8)
	{
	#pragma omp for nowait
		for (int i = 0; i < group.size(); ++i)
		{	
			if(group[i]<0)
				continue;	
			const float& a_i = getA_i(storage, group, array, i, object, normOption);
			const float& b_i = getB_i(storage, group, array, i, object, normOption);

			float s_i;
			if(abs(a_i-b_i)<1.0e-8)
				s_i = 0;
			else if(a_i<b_i)
				s_i = 1 - a_i/b_i;
			else
				s_i = b_i/a_i - 1;
			if(isnan(s_i))
			{
				std::cout << "Error for nan number!" << std::endl;
				exit(1);
			}
			sData[i] = s_i;

		#pragma omp critical
			sSummation += s_i;
		}
	}
	sAverage = sSummation/(group.size()-noise);

#pragma omp parallel for schedule(dynamic) num_threads(8)
	for (int i = 0; i < sCluster.size(); ++i)
	{
		float& eachCluster = sCluster[i];
		eachCluster = 0;
		const std::vector<int>& clustSet = storage[i];
		for (int j = 0; j < clustSet.size(); ++j)
		{
			eachCluster += sData[clustSet[j]];
		}
		eachCluster/=clustSet.size();
	}
}


const float Silhouette::getA_i(const std::vector<std::vector<int> >& storage,
							   const std::vector<int>& group,
						 	   const MatrixXf& array,
						 	   const int& index,
						 	   const MetricPreparation& object,
						 	   const int& normOption)
{
	const std::vector<int> clusterSet = storage[group[index]];
	float inClusterDist = 0.0;
	for (int j = 0; j < clusterSet.size(); ++j)
	{
		if(clusterSet[j]!=index)
		{
			if(distanceMatrix)
				inClusterDist += distanceMatrix[index][clusterSet[j]];
			else
				inClusterDist += getDist(index, clusterSet[j], object, array, normOption);
		}
	}
	if(isnan(inClusterDist))
	{
		std::cout << "a_i has nan error!" << std::endl;
		exit(1);
	}
	float a_i;
	if(clusterSet.size()==1)
		a_i = 0;
	else
		a_i = inClusterDist/(clusterSet.size()-1);
	return a_i;
}


const float Silhouette::getB_i(const std::vector<std::vector<int> >& storage,
							   const std::vector<int>& group,
							   const MatrixXf& array,
							   const int& index,
							   const MetricPreparation& object,
							   const int& normOption)
{
	float outClusterDist = FLT_MAX, perClusterDist = 0;
	std::vector<int> outClusterSet;
	if(storage.size()==1)
		return 0;
	for (int j = 0; j < storage.size(); ++j) //j is group no.
	{
		if(j!=group[index]) //the other cluster
		{
			outClusterSet = storage[j];//get integer list of this group
			perClusterDist = 0;
			for (int k = 0; k < outClusterSet.size(); ++k)
			{
				if(distanceMatrix)
					perClusterDist+=distanceMatrix[index][outClusterSet[k]];
				else
					perClusterDist += getDist(index, outClusterSet[k], object, array, normOption);
			}
			if(perClusterDist<0)
			{
				std::cout << "Error for negative distance!" << std::endl;
				exit(1);
			}
			perClusterDist/=outClusterSet.size();
			if(outClusterDist>perClusterDist)
				outClusterDist=perClusterDist;
		}
	}
	return outClusterDist;
}


const float Silhouette::getDist(const int& first,
								const int& second,
								const MetricPreparation& object,
								const MatrixXf& array,
								const int& normOption)
{
	float distance = getDisimilarity(array.row(first),array.row(second),
						   			 first,second,normOption,object);
	if(distance<0)
	{
		std::cout << "Error for negative distance!" << std::endl;
		exit(1);
	}
	return distance;
}