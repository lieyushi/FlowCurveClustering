#ifndef _PCA_CLUSTER_H
#define _PCA_CLUSTER_H

#include "IOHandler.h"
#include "Initialization.h"
#include "Silhouette.h"
#include "ValidityMeasurement.h"

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


/* record event list and time list */
struct TimeRecorder
{
	std::vector<string> eventList;
	std::vector<string> timeList;
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
									  TimeRecorder& tr,
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
									  TimeRecorder& tr,
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
									 TimeRecorder& tr,
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
									 TimeRecorder& tr,
									 Silhouette& sil);

private:

	static void performSVD(MatrixXf& cArray, 
						   const Eigen::MatrixXf& data, 
						   const int& Row, 
						   const int& Column, 
						   int& PC_Number, 
						   MatrixXf& SingVec,
	                       VectorXf& meanTrajectory,
						   TimeRecorder& tr);

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
								 TimeRecorder& tr,
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
											 TimeRecorder& tr,
											 Silhouette& sil);

};

#endif

