/*
 * This class is to calculate the silhouette, Gamma statistics and DB index for clustering evaluation
 */

#include "Silhouette.h"


/*
 * @brief Default constructor to create storage for the silhouette value
 */
Silhouette::Silhouette()
{
	sData = std::vector<float>();
	sCluster = std::vector<float>();
}


/*
 * @brief Destructor to empty the vector memory
 */
Silhouette::~Silhouette()
{
	reset();
}


/*
 * @brief Reset the member variables in the silhouette class
 */
void Silhouette::reset()
{
	sData.clear();
	sCluster.clear();
	sAverage = 0;
}


/*
 * @brief Compute the silhouette value with given input information
 * @param normOption: The norm type
 * @param array: The coordinate matrix of the streamlines
 * @param Row: The size of the row
 * @param Column: The column size
 * @param group: The labels for all the streamlines
 * @param object: A MetricPerparation object for distance matrix computation
 * @param groupNumber: how many clusters are as input
 * @param isPBF: A bool variable to tell whether the coordinates are from PBF fluid simulation or not
 */
void Silhouette::computeValue(const int& normOption,
						   	  const MatrixXf& array, 
						   	  const int& Row, 
						      const int& Column,
						   	  const std::vector<int>& group, 
						   	  const MetricPreparation& object,
						   	  const int& groupNumber,
							  const bool& isPBF)
{

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
	// compute the value
	computeValue(normOption, array, Row, Column, group, object, groupNumber, isPBF, storage);

}


/*
 * @brief Compute the silhouette value for PCA based input, so it will be calculated using Euclidean distance
 * @param array: The matrix coordinates of the dimensionality reduced space
 * @param group: The labels for different streamlines
 * @param groupNumber: The number of clusters that are formed
 * @param isPBF: The bool variable to tell whether it's PBF fluid simulation or not
 */
void Silhouette::computeValue(const Eigen::MatrixXf& array,
							  const std::vector<int>& group,
							  const int& groupNumber,
							  const bool& isPBF)
{
	sData.clear();
	sCluster.clear();

	/* get Row and Column information */
	const int& Row = array.rows();
	const int& Column = array.cols();

	sData = std::vector<float>(Row,0);

	/* assert information */
	assert(Row==group.size());

	std::vector<std::vector<int> > storage(groupNumber, std::vector<int>());
	for (int i = 0; i < group.size(); ++i)
	{
		storage[group[i]].push_back(i);
	}
	/* record labeling information */
	generateGroups(storage);

	//groupNumber doesn't include noise group
	sCluster = std::vector<float>(groupNumber, 0);

	/* if the silhouett computing is not for PBF dataset, then would use distanceMatrix */
	Eigen::MatrixXf distM, idealDistM;
	if(!isPBF)	// not from PBF, so the distance matrix can be assigned
	{
		getMatrixM(array,group,storage,distM,idealDistM);
	}

	std::cout << "Compute silhouette..." << std::endl;
	/* compute silhouette value */
	computeSilhouette(array, group, isPBF, storage, distM);

	std::cout << "silhouette is " << sAverage << std::endl;

	std::cout << "Compute DB index..." << std::endl;
	/* compute DB index */
	computeDBIndex(array, group, storage);
	std::cout << "DB index is " << dbIndex << std::endl;
	/* compute Gamma statistic for distM and idealDistM */
	if(!isPBF)	// only compute Gamma statistics for non-PBF data set
	{
		std::cout << "Compute gamma statistics..." << std::endl;
		computeGammaStatistic(distM,idealDistM);
		std::cout << "Gamma statistics is " << gammaStatistic << std::endl;
		/* garbage collection for eigen::matrix */
		distM.resize(0,0);
		idealDistM.resize(0,0);
	}

}


/*
 * @brief Compute the three clustering evaluation metrics for general norm input
 * @param normOption: The norm option
 * @param array: The matrix coordinates of the streamlines
 * @param Row: The row size
 * @param Column: The column size
 * @param group: The labels for all the streamlines
 * @param object: The MetricPreparation object for distance matrix computation
 * @param groupNumber: number of clusters as input
 * @param isPBF: The bool to say whether it is from PBF or not
 * @param storage: The candidates for each cluster
 */
void Silhouette::computeValue(const int& normOption,
							  const MatrixXf& array,
							  const int& Row,
							  const int& Column,
							  const std::vector<int>& group,
							  const MetricPreparation& object,
							  const int& groupNumber,
							  const bool& isPBF,
							  const std::vector<vector<int> >& storage)
{
	sData.clear();
	sCluster.clear();
	sData = std::vector<float>(Row,0);
	assert(Row==group.size());

	//groupNumber doesn't include noise group
	sCluster = std::vector<float>(groupNumber, 0);
	/* if the silhouett computing is not for PBF dataset, then would use distanceMatrix */
	Eigen::MatrixXf idealDistM;
	if(!isPBF)
	{
		getMatrixM(array,group,storage,idealDistM);
	}

	std::cout << "Compute silhouette..." << std::endl;
	/* compute silhouette value */
	computeSilhouette(array, group, storage, object, normOption);
	std::cout << "Silhouette is " << sAverage << std::endl;

	std::cout << "Compute DB index..." << std::endl;
	/* compute DB index */
	computeDBIndex(array, group, storage, object, normOption);
	std::cout << "DB index is " << dbIndex << std::endl;

	/* compute Gamma statistic for distM and idealDistM */
	if(!isPBF)
	{
		std::cout << "Compute gamma statistics..." << std::endl;
		computeGammaStatistic(idealDistM);
		std::cout << "Gamma statistics is " << gammaStatistic << std::endl;
		/* garbage collection for eigen::matrix */
		idealDistM.resize(0,0);
	}
}


/*
 * @brief Compute the A_i for silhouette value calculation
 * @param storage: The candidates included for all the clusters
 * @param group: The cluster labels of streamlines
 * @param array: The coordinates of the streamlines
 * @param index: The target streamline
 * @param object: The MetricPreparation object
 * @param normOption: The norm option
 */
const float Silhouette::getA_i(const std::vector<std::vector<int> >& storage,
							   const std::vector<int>& group,
						 	   const MatrixXf& array,
						 	   const int& index,
						 	   const MetricPreparation& object,
						 	   const int& normOption)
{
	const std::vector<int>& clusterSet = storage[group[index]];
	float inClusterDist = 0.0, dist;
	for (int j = 0; j < clusterSet.size(); ++j)
	{
		if(clusterSet[j]!=index)
		{
			if(distanceMatrix)
			{
				dist = distanceMatrix[index][clusterSet[j]];
			}
			else
			{
				dist = getDist(index, clusterSet[j], object, array, normOption);
			}
			inClusterDist += dist;
		}
	}
	if(std::isnan(inClusterDist))
	{
		std::cout << "a_i has nan error! " << inClusterDist << std::endl;
		exit(1);
	}
	float a_i;
	if(clusterSet.size()==1)
		a_i = 0;
	else
		a_i = inClusterDist/(clusterSet.size()-1);
	return a_i;
}


/*
 * @brief Compute the B_i for silhouette value calculation
 * @param storage: The candidates included for all the clusters
 * @param group: The cluster labels of streamlines
 * @param array: The coordinates of the streamlines
 * @param index: The target streamline
 * @param object: The MetricPreparation object
 * @param normOption: The norm option
 */
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


/*
 * @brief Compute the similarity distance between two streamlines
 * @param first: The index of first streamline
 * @param second: The index of second streamline
 * @param object: The MetricPreparation class object
 * @param array: The matrix coordinates
 * @param normOption: The norm option
 */
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
	if(isnan(distance) || isinf(distance))
	{
		std::cout << "Error for distance value that is nan or inf!" << std::endl;
		exit(1);
	}
	return distance;
}


/*
 * @brief Compute the distM (distance matrix) and idealDistM (ideal distance matrix) for PCA only
 * @param cArray: The matrix coordinates of the streamlines
 * @param group: The labels for all the streamlines
 * @param storage: The candidates inside all the clusters
 * @param distM: The distance matrix to be assigned the value
 * @param idealDist: The ideal distance matrix with 0 and 1
 */
void Silhouette::getMatrixM(const Eigen::MatrixXf& cArray,
		  	  				const std::vector<int>& group,
							const std::vector<std::vector<int> >& storage,
							Eigen::MatrixXf& distM,
							Eigen::MatrixXf& idealDistM)
{
	const int& Row = cArray.rows();
	const int& Column = cArray.cols();

	/* resize matrix size */
	distM = Eigen::MatrixXf::Zero(Row,Row);
	idealDistM = Eigen::MatrixXf::Constant(Row,Row,1.0);

	/* of course here is not related to distanceMatrix which is a global variable */
#pragma omp parallel for schedule(static) num_threads(8)
	for(int i=0;i<Row;++i)
	{
		for(int j=0;j<Row;++j)
		{
			if(i==j)
				continue;
			/* assign the Euclidean distance of cArray */
			distM(i,j)=(cArray.row(i)-cArray.row(j)).norm();
		}
	}

	const int& groupNumber = storage.size();
#pragma omp parallel for schedule(static) num_threads(8)
	for(int i=0;i<groupNumber;++i)
	{
		/* if i and j in same cluster, then set it to be zero */
		const std::vector<int>& eachVec = storage[i];
		const int& eachSize = eachVec.size();
		for(int j=0;j<eachSize;++j)
		{
			for(int k=0;k<eachSize;++k)
			{
				idealDistM(eachVec[j],eachVec[k]) = 0;
			}
		}
	}
}


/*
 * @brief Compute the distM idealDistM (ideal distance matrix) for non-PBF case and distM already been stored
 * @param cArray: The matrix coordinates of the streamlines
 * @param group: The labels for all the streamlines
 * @param storage: The candidates inside all the clusters
 * @param idealDist: The ideal distance matrix with 0 and 1
 */
void Silhouette::getMatrixM(const Eigen::MatrixXf& cArray,
		  	  				const std::vector<int>& group,
							const std::vector<std::vector<int> >& storage,
							Eigen::MatrixXf& idealDistM)
{
	const int& Row = cArray.rows();
	const int& Column = cArray.cols();

	/* resize matrix size */
	idealDistM = Eigen::MatrixXf::Constant(Row,Row,1.0);

	/* find the ideal matrix inside which the idealDistM(i,j)==0 only if i and j in same cluster */
	const int& groupNumber = storage.size();
#pragma omp parallel for schedule(static) num_threads(8)
	for(int i=0;i<groupNumber;++i)
	{
		/* if i and j in same cluster, then set it to be zero */
		const std::vector<int>& eachVec = storage[i];
		const int& eachSize = eachVec.size();
		for(int j=0;j<eachSize;++j)
		{
			for(int k=0;k<eachSize;++k)
			{
				idealDistM(eachVec[j],eachVec[k]) = 0;
			}
		}
	}
}


/*
 * @brief Compute the A_i for silhouette value calculation for PCA case (use Euclidean distance)
 * @param storage: The candidates included for all the clusters
 * @param group: The cluster labels of streamlines
 * @param array: The coordinates of the streamlines
 * @param index: The target streamline
 * @param isPBF: A bool tag whether the data is from PBF or not
 * @param distM: The calculated distance matrix that has been calculated before
 */
const float Silhouette::getA_i(const std::vector<std::vector<int> >& storage,
							   const std::vector<int>& group,
							   const Eigen::MatrixXf& array,
							   const int& index,
							   const bool& isPBF,
							   const Eigen::MatrixXf& distM)
{
	const std::vector<int>& clusterSet = storage[group[index]];
	float inClusterDist = 0.0;

	int candidate;
	for (int j = 0; j < clusterSet.size(); ++j)
	{
		candidate = clusterSet[j];
		if(candidate!=index)
		{
			if(!isPBF)
				inClusterDist += distM(index,candidate);
			else
				inClusterDist += (array.row(index)-array.row(candidate)).norm();
		}
	}
	if(std::isnan(inClusterDist))
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


/*
 * @brief Compute the B_i for silhouette value calculation for PCA case (use Euclidean distance)
 * @param storage: The candidates included for all the clusters
 * @param group: The cluster labels of streamlines
 * @param array: The coordinates of the streamlines
 * @param index: The target streamline
 * @param isPBF: A bool tag whether the data is from PBF or not
 * @param distM: The calculated distance matrix that has been calculated before
 */
const float Silhouette::getB_i(const std::vector<std::vector<int> >& storage,
							   const std::vector<int>& group,
							   const Eigen::MatrixXf& array,
							   const int& index,
							   const bool& isPBF,
							   const Eigen::MatrixXf& distM)
{
	float outClusterDist = FLT_MAX, perClusterDist = 0;
	std::vector<int> outClusterSet;

	int candidate, outClusterSize;
	for (int j = 0; j < storage.size(); ++j) //j is group no.
	{
		if(j!=group[index]) //the other cluster
		{
			outClusterSet = storage[j];//get integer list of this group
			perClusterDist = 0;

			outClusterSize = outClusterSet.size();

			/* empty cluster which is erroneous */
			if(outClusterSize==0)
			{
				std::cout << "Found empty clusters!" << std::endl;
				exit(1);
			}

			/* get average dist to all elements inside the cluster */
			for (int k = 0; k < outClusterSize; ++k)
			{
				candidate = outClusterSet[k];
				if(!isPBF)
					perClusterDist+=distM(index, candidate);
				else
					perClusterDist += (array.row(index)-array.row(candidate)).norm();
			}
			if(perClusterDist<0)
			{
				std::cout << "Error for negative distance!" << std::endl;
				exit(1);
			}
			perClusterDist/=outClusterSize;
			if(outClusterDist>perClusterDist)
				outClusterDist=perClusterDist;
		}
	}
	return outClusterDist;
}


/*
 * @brief Compute silhouette value for PCA-based only (no norm option)
 * @param array: The matrix coordinates of streamlines
 * @param group: The labels of all the streamlines
 * @param isPBF: The bool tag whether is PBF or not
 * @param storage: The candidates in all the clusters
 * @param distM: The distance matrix as input
 */
void Silhouette::computeSilhouette(const Eigen::MatrixXf& array,
								   const std::vector<int>& group,
								   const bool& isPBF,
								   const std::vector<std::vector<int> >& storage,
								   const Eigen::MatrixXf& distM)
{
	const int& Row = array.rows();

// compute Silhouette value for each data
	float sSummation = 0;

#pragma omp parallel num_threads(8)
	{
	#pragma omp for nowait
		for (int i = 0; i < Row; ++i)
		{
			const float& a_i = getA_i(storage, group, array, i, isPBF, distM);
			const float& b_i = getB_i(storage, group, array, i, isPBF, distM);

			float s_i;
			if(abs(a_i-b_i)<1.0e-8)
				s_i = 0;
			else if(a_i<b_i)
				s_i = 1 - a_i/b_i;
			else
				s_i = b_i/a_i - 1;
			if(std::isnan(s_i))
			{
				std::cout << "Error for nan number!" << std::endl;
				exit(1);
			}
			sData[i] = s_i;

		#pragma omp critical
			sSummation += s_i;
		}
	}
	sAverage = sSummation/group.size();

#pragma omp parallel for schedule(static) num_threads(8)
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


/*
 * @brief Compute silhouette value for general norm input
 * @param array: The matrix coordinates of streamlines
 * @param group: The labels of all the streamlines
 * @param storage: The candidates in all the clusters
 * @param object: The MetricPreparation class object
 * @param normOption: The norm option
 */
void Silhouette::computeSilhouette(const Eigen::MatrixXf& array,
					   	   	       const std::vector<int>& group,
								   const std::vector<std::vector<int> >& storage,
								   const MetricPreparation& object,
								   const int& normOption)
{
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
			if(std::isnan(s_i))
			{
				std::cout << "Error for nan number!" << std::endl;
				exit(1);
			}
			sData[i] = s_i;

		#pragma omp critical
			sSummation += s_i;
		}
	}
	sAverage = sSummation/(group.size());

#pragma omp parallel for schedule(static) num_threads(8)
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


/*
 * @brief Compute DB index for PCA case (only Euclidean distance used)
 * @param array: The matrix coordinates of streamlines
 * @param group: The labels of all the streamlines
 * @param storage: The candidates in all the clusters
 */
void Silhouette::computeDBIndex(const Eigen::MatrixXf& array,
								const std::vector<int>& group,
								const std::vector<std::vector<int> >& storage)
{
	dbIndex = 0.0;

	const int& groupNumber = storage.size();

	const int& Column = array.cols();

	/* calculated the projected-space cenroid */
	Eigen::MatrixXf centroid(groupNumber, Column);

	/* average distance of all elements in cluster to its centroid */
	Eigen::VectorXf averageDist(groupNumber);

#pragma omp parallel for schedule(static) num_threads(8)
	for(int i=0;i<groupNumber;++i)
	{
		Eigen::VectorXf tempCentroid = Eigen::VectorXf::Zero(Column);

		const std::vector<int>& clusterVec = storage[i];
		const int& clusterSize = clusterVec.size();

		for(int j=0;j<clusterSize;++j)
			tempCentroid+=array.row(clusterVec[j]);

		/* get the centroid coordinates */
		centroid.row(i) = tempCentroid/clusterSize;

		float inClusterSum = 0.0, temp_dist;
		for(int j=0;j<clusterSize;++j)
		{
			//inClusterSum+=getDisimilarity(centroid.row(i),array,clusterVec[j],normOption,object);
			inClusterSum+=(array.row(clusterVec[j])-centroid.row(i)).norm();
		}
		averageDist(i) = inClusterSum/clusterSize;
	}

#pragma omp parallel num_threads(8)
	{
	#pragma omp for nowait
		for (int i = 0; i < groupNumber; ++i)
		{
			float maxValue = (float)INT_MIN, ratioDist;
			for (int j=0;j<groupNumber;++j)
			{
				if(i==j)
					continue;
				ratioDist = (averageDist(i)+averageDist(j))/(centroid.row(i)-centroid.row(j)).norm();

				if(maxValue<ratioDist)
					maxValue=ratioDist;
			}

		#pragma omp critical
			dbIndex += maxValue;
		}
	}
	dbIndex/=groupNumber;
}


/*
 * @brief Compute DB index for general input of norm option
 * @param array: The matrix coordinates of streamlines
 * @param group: The labels of all the streamlines
 * @param storage: The candidates in all the clusters
 * @param object: The MetricPreparation class object
 * @param normOption: The norm option
 */
/* compute DB index with normOption as input */
void Silhouette::computeDBIndex(const Eigen::MatrixXf& array,
								const std::vector<int>& group,
								const std::vector<std::vector<int> >& storage,
								const MetricPreparation& object,
								const int& normOption)
{
	dbIndex = 0.0;

	const int& groupNumber = storage.size();

	const int& Column = array.cols();

	/* calculated the projected-space cenroid */
	Eigen::MatrixXf centroid(groupNumber, Column);

	/* average distance of all elements in cluster to its centroid */
	Eigen::VectorXf averageDist(groupNumber);

#pragma omp parallel for schedule(static) num_threads(8)
	for(int i=0;i<groupNumber;++i)
	{
		Eigen::VectorXf tempCentroid = Eigen::VectorXf::Zero(Column);

		const std::vector<int>& clusterVec = storage[i];
		const int& clusterSize = clusterVec.size();

		for(int j=0;j<clusterSize;++j)
			tempCentroid+=array.row(clusterVec[j]);

		/* get the centroid coordinates */
		centroid.row(i) = tempCentroid/clusterSize;

		float inClusterSum = 0.0, temp_dist;
		for(int j=0;j<clusterSize;++j)
		{
			inClusterSum+=getDisimilarity(centroid.row(i),array,clusterVec[j],normOption,object);
		}
		averageDist(i) = inClusterSum/clusterSize;
	}

#pragma omp parallel num_threads(8)
	{
	#pragma omp for nowait
		for (int i = 0; i < groupNumber; ++i)
		{
			float maxValue = (float)INT_MIN, ratioDist, centDist;
			for (int j=0;j<groupNumber;++j)
			{
				if(i==j)
					continue;
				centDist = getDisimilarity(centroid.row(i), centroid.row(j), normOption, object);

				//ratioDist = (averageDist(i)+averageDist(j))/(centroid.row(i)-centroid.row(j)).norm();
				ratioDist = (averageDist(i)+averageDist(j))/centDist;

				if(maxValue<ratioDist)
					maxValue=ratioDist;
			}

		#pragma omp critical
			dbIndex += maxValue;
		}
	}
	dbIndex/=groupNumber;
}


/*
 * @brief Compute the Gamma statistics between two matrices
 * @param distM: The distance matrix
 * @param idealDistM: The ideal distance matrix with only value 1 and 0
 */
void Silhouette::computeGammaStatistic(const Eigen::MatrixXf& distM,
									   const Eigen::MatrixXf& idealDistM)
{
	const int& Row = distM.rows();

	const int& totalNum = Row*(Row-1)/2;
	/* mean of values */
	double u_1 = 0.0, u_2 = 0.0;

	/* E(X*X) */
	double s_1 = 0.0, s_2 = 0.0, numerator = 0.0;

	for(int i=0;i<Row-1;++i)
	{
		for(int j=i+1;j<Row;++j)
		{
			/* update the mean u_1, u_2 */
			u_1+=distM(i,j);
			u_2+=idealDistM(i,j);

			/* update the numerator */
			numerator+=distM(i,j)*idealDistM(i,j);

			/* update the deviation */
			s_1+=distM(i,j)*distM(i,j);
			s_2+=idealDistM(i,j)*idealDistM(i,j);
		}
	}

	/* get mean of distM and idealDistM */
	u_1/=totalNum;
	u_2/=totalNum;

	/* get numerator for the computing */
	numerator-=totalNum*u_1*u_2;

	/* get standard deviation */
	s_1=sqrt(s_1/totalNum-u_1*u_1);
	s_2=sqrt(s_2/totalNum-u_2*u_2);

	if(std::isnan(s_1) || std::isnan(s_2))
	{
		std::cout << "standard deviation has nan error!" << std::endl;
		exit(1);
	}
	gammaStatistic = float(numerator/s_1/s_2/totalNum);
}


/*
 * @brief Compute the Gamma statistics for the general input of norm options
 * @param idealDistM: The ideal distance matrix with 1 and 0 only
 */
void Silhouette::computeGammaStatistic(const Eigen::MatrixXf& idealDistM)
{
	const int& Row = idealDistM.rows();

	const int& totalNum = Row*(Row-1)/2;

	/* mean of values */
	double u_1 = 0.0, u_2 = 0.0;

	/* E(X*X) */
	double s_1 = 0.0, s_2 = 0.0, numerator = 0.0;

	for(int i=0;i<Row-1;++i)
	{
		for(int j=i+1;j<Row;++j)
		{
			u_1+=distanceMatrix[i][j];
			u_2+=idealDistM(i,j);

			numerator+=distanceMatrix[i][j]*idealDistM(i,j);

			s_1+=distanceMatrix[i][j]*distanceMatrix[i][j];
			s_2+=idealDistM(i,j)*idealDistM(i,j);
		}
	}

	/* get mean of distM and idealDistM */
	u_1/=totalNum;
	u_2/=totalNum;

	/* get numerator for the computing */
	numerator-=totalNum*u_1*u_2;

	/* get standard deviation */
	s_1=sqrt(s_1/totalNum-u_1*u_1);
	s_2=sqrt(s_2/totalNum-u_2*u_2);

	if(std::isnan(s_1) || std::isnan(s_2))
	{
		std::cout << "standard deviation has nan error!" << std::endl;
		exit(1);
	}
	gammaStatistic = float(numerator/s_1/s_2/totalNum);
}

