/*
 * DetermClusterNum.cpp, which is the C++ implementation of the hierarchical L-Method
 *
 *  Created on: Aug 22, 2018
 *      Author: lieyu
 */

#include "DetermClusterNum.h"


/*
 * @brief The default constructor
 */
DetermClusterNum::DetermClusterNum() {
	// TODO Auto-generated constructor stub

}


/*
 * @brief Destructor
 */
DetermClusterNum::~DetermClusterNum() {
	// TODO Auto-generated destructor stub
}


/*
 * @brief Use iterative refinement of knee to get optimal number for hierarchical clustering
 * @param eval_graph: The map that contains the cluster number and its merged distance
 */
void DetermClusterNum::iterativeRefinement(std::map<int, float>& eval_graph)
{
	// some necessary pre-processing to remove irregular shapes for the L-method
	removeExtreme(eval_graph);

	// start from the first to search the point with knee
	int cutoff, lastKnee;
	int currentKnee = eval_graph.rbegin()->first;
	cutoff = currentKnee;
	do	// an iterative refinement for the L-method
	{
		lastKnee = currentKnee;
		currentKnee = LMethod(eval_graph, cutoff);
		std::cout << "returned value is " << currentKnee <<", cutoff is " << cutoff << std::endl;
		cutoff = currentKnee*2;
	}while(currentKnee < lastKnee);

	// get the optimal number of clusters
	finalNumOfClusters = currentKnee;

	std::cout << finalNumOfClusters << std::endl;
}


/*
 * @brief Find the knee value by the L-method with a given cutoff value
 * @param eval_graph: The map with cluster numbers and their relative merged distance
 * @param cutoff: The cutoff index point
 * @return An index found by the Lmethod which is related to the knee
 */
const int DetermClusterNum::LMethod(const std::map<int, float>& eval_graph, const int& cutoff)
{
	struct CompObj { float val; int index; };
// #pragma omp declare reduction(minimum : struct CompObj : omp_out = omp_in.val < omp_out.val ? omp_in : omp_out)
	struct CompObj RMSE;
	RMSE.val = FLT_MAX;
	RMSE.index = -1;

	const int& firstIndex = eval_graph.begin()->first;
	/* find the minimal c that minimizes RMSE for the selected cutoff */
#pragma omp parallel num_threads(8)
	{
	#pragma omp nowait
		for(int i=firstIndex;i<=cutoff;++i)
		{
			/* left segment linear least square fitting */
			std::vector<float> index_vec;
			std::vector<float> dist_vec;

			// assign the vector for left segment
			std::map<int, float>::const_iterator iter;
			for(int j=firstIndex;j<=i;++j)
			{
				iter = eval_graph.find(j);
				if(iter!=eval_graph.end())
				{
					index_vec.push_back(iter->first);
					dist_vec.push_back(iter->second);
				}
			}
			Eigen::MatrixXf A_sub(2, index_vec.size());
			A_sub.row(0) = Eigen::VectorXf::Map(&(index_vec[0]), index_vec.size()).transpose();
			A_sub.row(1) = Eigen::VectorXf::Constant(index_vec.size(), 1.0).transpose();
			Eigen::VectorXf b_sub = Eigen::VectorXf::Map(&(dist_vec[0]), index_vec.size());
			A_sub.transposeInPlace();
			int firstRows = A_sub.rows();

			// solve the least-square fitting problems
			Eigen::VectorXf c = A_sub.colPivHouseholderQr().solve(b_sub);
			Eigen::VectorXf error = b_sub-A_sub*c;
			float rmse_l = error.transpose()*error;

			/* right segment linear least square fitting */
			index_vec.clear();
			dist_vec.clear();

			// assignment of the vector
			for(int j=i+1;j<=cutoff;++j)
			{
				iter = eval_graph.find(j);
				if(iter!=eval_graph.end())
				{
					index_vec.push_back(iter->first);
					dist_vec.push_back(iter->second);
				}
			}
			A_sub = Eigen::MatrixXf(2, index_vec.size());
			A_sub.row(0) = Eigen::VectorXf::Map(&(index_vec[0]), index_vec.size()).transpose();
			A_sub.row(1) = Eigen::VectorXf::Constant(index_vec.size(), 1.0).transpose();
			b_sub = Eigen::VectorXf::Map(&(dist_vec[0]), index_vec.size());
			A_sub.transposeInPlace();
			int secondRows = A_sub.rows();

			// least-square fitting problem
			c = A_sub.colPivHouseholderQr().solve(b_sub);
			error = b_sub-A_sub*c;
			float rmse_r = error.transpose()*error;

			/* compute the total weighted error */
			float rmse = float(firstRows)/float(firstRows+secondRows)*rmse_l+
					float(secondRows)/float(firstRows+secondRows)*rmse_r;
			// update the rmse value and index
		#pragma omp critical
			if(RMSE.val>rmse)
			{
				RMSE.val=rmse;
				RMSE.index=i;
			}
		}
	}

	return RMSE.index;
}


/*
 * @brief Record the L-method result in the local file
 * @param normOption: The norm option as input
 */
void DetermClusterNum::recordLMethodResult(const int& normOption)
{
	std::ofstream readme("../dataset/LMethod",ios::out | ios::app);
	if(!readme)
	{
		std::cout << "Error creating readme!" << std::endl;
		exit(1);
	}
	readme << "Optimal cluster number of norm " << normOption << " is " << finalNumOfClusters << std::endl;
	readme << std::endl;
	readme.close();
}


/*
 * @brief Remove extremely dissimilarity mcluster merges
 * @param eval_graph: The map including the cluster numbers and their merged distance
 */
void DetermClusterNum::removeExtreme(std::map<int, float>& eval_graph)
{
	// find the left index with the maximal distance
	float maxDist = -1.0;
	int leftIndex = -1;
	for(auto iter:eval_graph)
	{
		if(maxDist<iter.second)
		{
			maxDist=iter.second;
			leftIndex=iter.first;
		}
	}
	auto iter_index = eval_graph.find(leftIndex);

	// remove some irregular indices
	for(auto iter=eval_graph.begin();iter!=iter_index;)
	{
		if(iter->first<leftIndex&&iter->second<maxDist)
			eval_graph.erase(iter++);
		else
			++iter;
	}
}
