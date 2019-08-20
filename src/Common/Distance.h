/*
 * It contains the function to calculate the distance values/similarity measures among integral curves
 */

#ifndef _DISTANCE_H_
#define _DISTANCE_H_

#include "Metric.h"


// the constant for computing the Procrutes distance
extern const int& PROCRUSTES_SIZE;

// the global pointer variable distanceMatrix
extern float** distanceMatrix;


/* ------------------ Compute norm 3 for integral curves ------------------------- */

const float getBMetric_3(const VectorXf& row,
						 const int& size,
						 const int& i,
						 const std::vector<std::vector<float> >& rotationSequence
						);

const float getBMetric_3(const VectorXf& firstRow,
						 const int& size,
						 const VectorXf& secondRow
					    );
const float getBMetric_3(const int& first,
						 const int& second,
						 const std::vector<std::vector<float> >& rotationSequence
					    );
/* ------------------ Finish norm 3 for integral curves ------------------------- */



/* ------------------ Compute norm 6 for integral curves ------------------------- */
const float getBMetric_6(const VectorXf& row,
						 const int& size,
						 const int& i,
						 const std::vector<MultiVariate>& normalMultivariate
						);

const float getBMetric_6(const VectorXf& firstRow,
						 const int& size,
						 const VectorXf& secondRow
						);
const float getBMetric_6(const int& first,
						 const int& second,
						 const std::vector<MultiVariate>& normalMultivariate
						);
/* ------------------ Finish norm 6 for integral curves ------------------------- */


/* ------------------ Compute norm 7 for integral curves ------------------------- */
const float getBMetric_7(const VectorXf& row,
						 const int& size,
						 const int& i,
						 const std::vector<std::vector<float> >& rotationSequence
						);
const float getBMetric_7(const VectorXf& firstRow,
						 const int& size,
						 const VectorXf& secondRow
						);
const float getBMetric_7(const int& first,
						 const int& second,
						 const std::vector<std::vector<float> >& rotationSequence
						);
/* ------------------ Finish norm 7 for integral curves ------------------------- */


/* ------------------ Compute norm 9 for integral curves ------------------------- */
const float getBMetric_9(const VectorXf& row,
						 const int& size,
						 const int& i,
						 const std::vector<MultiVariate>& normalMultivariate
						);
const float getBMetric_9(const VectorXf& firstRow,
						 const int& size,
						 const VectorXf& secondRow
						);
const float getBMetric_9(const int& first,
						 const int& second,
						 const std::vector<MultiVariate>& normalMultivariate
						);
/* ------------------ Finish norm 9 for integral curves ------------------------- */


/* ------------------ Compute B-metric for integral curves ------------------------- */
const float getBMetric(const std::vector<float>& first, 
					   const std::vector<float>& second);

const float getBMetric(const MultiVariate& first, 
					   const MultiVariate& second);
/* ------------------ Finish B-metric for integral curves ------------------------- */


/* ------------------ Compute norm 10 for integral curves ------------------------- */
const float getMetric_10(const VectorXf& centroid,
						 const int& size,
						 const int& index,
						 const std::vector<VectorXf>& unitLength);
const float getMetric_10(const VectorXf& firstRow,
						 const int& size,
						 const VectorXf& secondRow);
const float getMetric_10(const int& first,
						 const int& second,
						 const std::vector<VectorXf>& unitLength);
/* ------------------ Finish norm 10 for integral curves ------------------------- */


// compute MCP distance
const float getMetric_MOP(const VectorXf& first, const VectorXf& second);


// compute Hausdorff distance
const float getMetric_Hausdorff(const VectorXf& first, const VectorXf& second);


/* ------------------ Compute norm 0, 11, 1, 2, 5, 8 for trajectories ------------------------- */
const float getNorm(const Eigen::VectorXf& centroid,
					const Eigen::VectorXf& r2,
					const int& index,
					const int& normOption,
					const std::vector<std::vector<float> >& object,
					const std::vector<std::vector<float> >& objectNorm);
const float getNorm(const VectorXf& centroid,
					const VectorXf& r2,
					const int& firstIndex,
					const int& secondIndex,
					const int& normOption,
					const std::vector<std::vector<float> >& object,
					const std::vector<std::vector<float> >& objectNorm);
const float getNorm(const VectorXf& r1, 
					const VectorXf& r2,
					const int& normOption);
/* ------------------ Finish norm 0, 11, 1, 2, 5, 8 for trajectories ------------------------- */


/* ------------------ Compute norm 14 for trajectories ------------------------- */

/* get signature-based dissimilarity metric given two elements and their histogram*/
const float getSignatureMetric(const Eigen::VectorXf& firstArray,
							   const Eigen::VectorXf& secondArray,
							   const std::vector<float>& firstHist,
							   const std::vector<float>& secondHist);


/* get signature-based dissimilarity metric given centroid */
const float getSignatureMetric(const Eigen::VectorXf& centroid,
							   const Eigen::VectorXf& first,
							   const std::vector<float>& firstHist);


/* get signature-based dissimilarity metric given two centroids */
const float getSignatureMetric(const Eigen::VectorXf& first,
							   const Eigen::VectorXf& second);
/* ------------------ Finish norm 14 for trajectories ------------------------- */


/* ------------------ Compute norm 15 for trajectories ------------------------- */

/* get adapted Procrustes distance. For example, if vec has 100 points, it will calculate mean of 94 point */
const float getProcrustesMetric(const Eigen::VectorXf& first,
								const Eigen::VectorXf& second);


/* get adapted Procrustes distance. For example, if vec has 100 points, it will calculate mean of 14 points */
const float getProcrustesMetricSegment(const Eigen::VectorXf& first,
									   const Eigen::VectorXf& second);
/* ------------------ Finish norm 15 for trajectories ------------------------- */


/* ------------------ Compute norm 16 for trajectories ------------------------- */

/* get illustrative visualization metric for paper An Illustrative Visualization Framework for 3D Vector Fields */
const float getEntropyMetric(const std::vector<float>& firstEntropy,
		                     const std::vector<float>& secondEntropy);

/* get illustrative visualization metric for paper An Illustrative Visualization Framework for 3D Vector Fields */
const float getEntropyMetric(const std::vector<float>& firstEntropy,
		                     const Eigen::VectorXf& array);

/* get illustrative visualization metric for paper An Illustrative Visualization Framework for 3D Vector Fields,
 * given two coordinate vectors */
const float getEntropyMetric(const Eigen::VectorXf& first,
		                     const Eigen::VectorXf& second);
/* ------------------ Finish norm 16 for trajectories ------------------------- */


/* get the revised MCP distance for pathlines from paper
 https://www.sciencedirect.com/science/article/pii/S0097849318300128
 */
const float getPathline_MCP(const Eigen::VectorXf& first,
        					const Eigen::VectorXf& second);


/*
 * The most important API to calculate the similarity measures
 */
const float getDisimilarity(const MatrixXf& data,
							const int& first,
							const int& second,
							const int& normOption,
							const MetricPreparation& object);

/*
 * The most important API to calculate the similarity measures
 */
const float getDisimilarity(const VectorXf& others,
							const MatrixXf& data,
							const int& index,
							const int& normOption,
							const MetricPreparation& object);

/*
 * The most important API to calculate the similarity measures
 */
const float getDisimilarity(const VectorXf& first,
							const VectorXf& second,
							const int& firstIndex,
							const int& secondIndex,
							const int& normOption,
							const MetricPreparation& object);


// compute the distance matrix
void getDistanceMatrix(const MatrixXf& data,
				       const int& normOption,
					   const MetricPreparation& object);


/*
 * The most important API to calculate the similarity measures
 */
const float getDisimilarity(const VectorXf& first,
							const VectorXf& second,
							const int& normOption,
							const MetricPreparation& object);


// clean up the distance matrix
void deleteDistanceMatrix(const int& Row);


/* get rotation for a series of streamlines */
const float getRotation(const std::vector<vector<float> >& streamline, std::vector<float>& rotation);


/* need to store each label for elements for NID computation */
void generateGroups(const std::vector<std::vector<int> >& storage);


#endif
