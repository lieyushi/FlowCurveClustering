#ifndef _PCA_CLUSTER_H
#define _PCA_CLUSTER_H

#include "../Common/IOHandler.h"
#include "../Common/Initialization.h"
#include "../Common/Silhouette.h"


struct Ensemble
{
	int number;
	int newIndex;
	int oldIndex;
	Ensemble(const int& number, const int& index):number(number),newIndex(-1), oldIndex(index)
	{}
	bool operator<(const Ensemble& object) const
	{
		return number<object.number;
	}
};


class PCA_Cluster
{
public: 

	static void performPCA_Clustering(const Eigen::MatrixXf& data, 
									  const int& Row, 
									  const int& Column, 
									  std::vector<MeanLine>& massCenter,
		 							  std::vector<int>& group, 
		 							  std::vector<int>& totalNum, 
		 							  std::vector<ExtractedLine>& closest,
		 							  std::vector<ExtractedLine>& furthest,
		 							  float& entropy,
		 							  Silhouette& sil);

	static void performPCA_Clustering(const Eigen::MatrixXf& data, 
									  const int& Row, 
									  const int& Column, 
									  std::vector<MeanLine>& massCenter, 
									  std::vector<int>& group, 
									  std::vector<int>& totalNum, 
									  std::vector<ExtractedLine>& closest, 
									  std::vector<ExtractedLine>& furthest, 
									  const int& Cluster, 
									  float& entropy,
									  Silhouette& sil);

	static void performDirectK_Means(const Eigen::MatrixXf& data, 
									 const int& Row, 
									 const int& Column, 
									 std::vector<MeanLine>& massCenter,
									 std::vector<int>& group, 
									 std::vector<int>& totalNum, 
									 std::vector<ExtractedLine>& closest,
									 std::vector<ExtractedLine>& furthest, 
									 const int& normOption,
									 float& entropy,
									 Silhouette& sil);

	static void performDirectK_Means(const Eigen::MatrixXf& data, 
									 const int& Row, 
									 const int& Column, 
									 std::vector<MeanLine>& massCenter,
									 std::vector<int>& group, 
									 std::vector<int>& totalNum, 
									 std::vector<ExtractedLine>& closest, 
									 std::vector<ExtractedLine>& furthest, 
									 const int& Cluster, 
									 const int& normOption,
									 float& entropy,
									 Silhouette& sil);

private:

	static void performSVD(MatrixXf& cArray, 
						   const Eigen::MatrixXf& data, 
						   const int& Row, 
						   const int& Column, 
						   int& PC_Number, 
						   MatrixXf& SingVec,
	                       VectorXf& meanTrajectory);

	static void performPC_KMeans(const MatrixXf& cArray, 
								 const int& Row,
								 const int& Column, 
								 const int& PC_Number, 
								 const MatrixXf& SingVec, 
								 const VectorXf& meanTrajectory, 
								 std::vector<MeanLine>& massCenter, 
								 const int& Cluster, 
								 std::vector<int>& group,
								 std::vector<int>& totalNum, 
								 std::vector<ExtractedLine>& closest, 
								 std::vector<ExtractedLine>& furthest, 
								 const Eigen::MatrixXf& data,
								 float& entropy,
								 Silhouette& sil);

	static void performAHC(const MatrixXf& cArray, 
						   const int& Row, 
						   const int& Column, 
						   const int& PC_Number, 
						   const MatrixXf& SingVec, 
		                   const VectorXf& meanTrajectory, 
		                   MatrixXf& clusterCenter, 
		                   std::vector<MeanLine>& massCenter);

	static void performFullK_MeansByClusters(const Eigen::MatrixXf& data, 
											 const int& Row, 
											 const int& Column, 
											 std::vector<MeanLine>& massCenter,
						   					 const int& Cluster, 
						   					 std::vector<int>& group, 
						   					 std::vector<int>& totalNum, 
						   					 std::vector<ExtractedLine>& closest, 
						   					 std::vector<ExtractedLine>& furthest, 
						   					 const int& normOption,
						   					 float& entropy,
						   					 Silhouette& sil);

};

#endif

