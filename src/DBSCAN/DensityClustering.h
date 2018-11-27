#ifndef _DENSITYCLUSTERING_H
#define _DENSITYCLUSTERING_H

#include "IOHandler.h"
#include "Initialization.h"
#include "Silhouette.h"
#include "ValidityMeasurement.h"
#include <queue>


enum PointType
{
	CORE = 0,
	BORDER,
	NOISE
};

struct PointNode
{
	int type;
	bool visited;
	int group;
	PointNode():type(-1), visited(false), group(-1)
	{}

	~PointNode()
	{}
};

struct DataSet
{
	vector<vector<float> > dataVec;
	Eigen::MatrixXf dataMatrix;
	int maxElements;
	int vertexCount;
	int dimension;

	string strName;
	string fullName;

	DataSet()
	{}

	~DataSet()
	{}
};


class DensityClustering
{
public:
	DensityClustering(const int& argc,
					  char **argv);

	~DensityClustering();

	void performClustering();


private:
	vector<PointNode> nodeVec;
	MetricPreparation object;
	int normOption;
	DataSet ds;

	/* whether it is a PBF dataset or not */
	bool isPBF;

	bool isPathlines;

	void setDataset(const int& argc,
				    char **argv);	//get dataset from file

	void setNormOption();	//get normOption from user input

	void DBSCAN(const float& radius_eps,
				const int& minPts);

	void expandCluster(const int& index,
					   vector<int>& neighbor,
					   const int& cluster_id,
					   const float& radius_eps,
					   const int& minPts);

	const vector<int> regionQuery(const int& index,
								  const float& radius_eps);

	void getDistRange(float& minDist, 
					  float& maxDist);

	const int setMinPts();

	const float setTimesMin(const float& minDist, 
					  		const float& maxDist);

	void extractFeatures(const float& radius_eps,
							   const int& minPts);

	/* two ways for getting eps, one is for user input, the other is average minPts-th dist for all candidates */
	const float getDistThreshold(const int& minPts);

	/* compute minPts-th dist for all candidates */
	const float getAverageDist(const int& minPts);
};

#endif
