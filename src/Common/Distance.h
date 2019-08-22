/*
 * @brief It contains the function to calculate the distance values/similarity measures among integral curves
 * @details
 * 	 0: Euclidean Norm
	 1: Fraction Distance Metric
	 2: piece-wise angle average
	 3: Bhattacharyya metric for rotation
	 4: average rotation
	 5: signed-angle intersection
	 6: normal-direction multivariate distribution
	 7: Bhattacharyya metric with angle to a fixed direction
	 8: Piece-wise angle average \times standard deviation
	 9: normal-direction multivariate un-normalized distribution
	 10: x*y/|x||y| borrowed from machine learning
	 11: cosine similarity
	 12: Mean-of-closest point distance (MCP)
	 13: Hausdorff distance min_max(x_i,y_i)
	 14: Signature-based measure from http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6231627
	 15: Procrustes distance take from http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6787131
	 16: entropy-based distance metric taken from http://vis.cs.ucdavis.edu/papers/pg2011paper.pdf
	 17: time-series MCP distance from https://www.sciencedirect.com/science/article/pii/S0097849318300128
			for pathlines only

 * @author Lieyu Shi
 */


#ifndef _DISTANCE_H_
#define _DISTANCE_H_

#include "Metric.h"


/*
 * @brief The length of neighboring vertices for calculating the Procrutes distance
 */
extern const int& PROCRUSTES_SIZE;


/*
 * @brief the global pointer variable for the distance matrix
 */
extern float** distanceMatrix;


/*
 * @brief Calculate the similarity measure 3 for the integral curves
 *
 * @param[in] row Input vector
 * @param[in] size The size of vector
 * @param[in] i The index
 * @param[in] rotationSequence The pre-calculated rotationSequence
 * @return The float value of similarity measure labeled as 3
 */
const float getBMetric_3(const VectorXf& row,
						 const int& size,
						 const int& i,
						 const std::vector<std::vector<float> >& rotationSequence
						);


/*
 * @brief Calculate the distance 3 with given two center trajectories for distance measuring
 *
 * @param[in] firstRow The coordinate of the first line
 * @param[in] size The size of the coordinate
 * @param[in] secondRow The coordinate of the second line
 * @return The float value of the similarity measure labeled as 3
 */
const float getBMetric_3(const VectorXf& firstRow,
						 const int& size,
						 const VectorXf& secondRow
					    );


/*
 * @brief Calculate the distance 3 with given the rotationSequence and two indices
 *
 * @param[in] first The index of the first line
 * @param[in] second The index of the second line
 * @param[in] rotationSeuqnce The rotation sequence vector of the lines as input
 * @return The float value of the similarity measure labeled as 3
 */
const float getBMetric_3(const int& first,
						 const int& second,
						 const std::vector<std::vector<float> >& rotationSequence
					    );


/*
 * @brief Calculate the norm 6 for integral curves given a center trajectory and index of pre-stored vector
 *
 * @param[in] row The coordinate of the first line
 * @param[in] size The size of the coordinate
 * @param[in] i The index of the second line
 * @param[in] normalMultivariate The normalized multivariate vector for all the streamlines
 * @return The float value of norm 6 between two integral curves
 */
const float getBMetric_6(const VectorXf& row,
						 const int& size,
						 const int& i,
						 const std::vector<MultiVariate>& normalMultivariate
						);


/*
 * @brief Calculate the norm 6 for integral curves given two center trajectories for distance measuring
 *
 * @param[in] firstRow The coordinate of the first line
 * @param[in] size The size of the coordinate
 * @param[in] secondRow The coordinate of the second line
 * @return The float value of norm 6 between two integral curves
 */
const float getBMetric_6(const VectorXf& firstRow,
						 const int& size,
						 const VectorXf& secondRow
						);


/*
 * @brief Calculate the norm 6 for integral curves given two indices for distance measuring
 *
 * @param[in] first The index of the first line
 * @param[in] second The index of the second line
 * @param[in] normalMultivariate The normalized multivariate vector for all the streamlines
 * @return The float value of norm 6 between two integral curves
 */
const float getBMetric_6(const int& first,
						 const int& second,
						 const std::vector<MultiVariate>& normalMultivariate
						);


/*
 * @brief Calculate the norm 7 for two integral curves given a center trajectory and index of pre-stored vector
 *
 * @param[in] row The coordinate of the first line
 * @param[in] size The size of the coordinates
 * @param[in] i The index of the second line
 * @param[in] rotationSequence The rotation sequence vector for all the streamlines
 * @return The float value of norm 7 between two integral curves
 */
const float getBMetric_7(const VectorXf& row,
						 const int& size,
						 const int& i,
						 const std::vector<std::vector<float> >& rotationSequence
						);


/*
 * @brief Calculate the norm 7 for two integral curves given two center trajectories for distance measuring
 *
 * @param[in] firstRow The coordinate of the first line
 * @param[in] size The size of the coordinates
 * @param[in] secondRow The coordinate of the second line
 * @return The float value of norm 7 between two integral curves
 */
const float getBMetric_7(const VectorXf& firstRow,
						 const int& size,
						 const VectorXf& secondRow
						);


/*
 * @brief Calculate the norm 7 for two integral curves given two indices of the lines for distance measuring
 *
 * @param[in] first The index of the first line
 * @param[in] second The index of the second line
 * @param[in] rotationSequence The rotation sequence vector for all the streamlines
 * @return The float value of norm 7 between two integral curves
 */
const float getBMetric_7(const int& first,
						 const int& second,
						 const std::vector<std::vector<float> >& rotationSequence
						);


/*
 * @brief Calculate the norm 9 for two integral curves given a center trajectory and index of pre-stored vector
 *
 * @param[in] row The coordinate of the first line
 * @param[in] size The size of the coordinates
 * @param[in] i The index of the second line
 * @param[in] normalMultivariate The normalized multivariate vector for all the streamlines
 * @return The float value of norm 9 between two integral curves
 */
const float getBMetric_9(const VectorXf& row,
						 const int& size,
						 const int& i,
						 const std::vector<MultiVariate>& normalMultivariate
						);


/*
 * @brief Calculate the norm 9 for two integral curves given two center trajectories for distance measuring
 *
 * @param[in] firstRow The coordinate of the first line
 * @param[in] size The size of the coordinate
 * @param[in] secondRow The coordinate of the second line
 * @return The float value of norm 9 between two integral curves
 */
const float getBMetric_9(const VectorXf& firstRow,
						 const int& size,
						 const VectorXf& secondRow
						);


/*
 * @brief Calculate the norm 9 for two integral curves given two indices of streamlines for distance measuring
 *
 * @param[in] first The index of the first line
 * @param[in] second The index of the second line
 * @param[in] normalMultivariate The normalized multivariate vector for all the streamlines
 * @return The float value of norm 9 between two integral curves
 */
const float getBMetric_9(const int& first,
						 const int& second,
						 const std::vector<MultiVariate>& normalMultivariate
						);


/*
 * @brief Calculate the B-metric for two streamlines considered as univariate normal distributions
 * @param[in] firstNorm3 The coordinate of the first line
 * @param[in] secondNorm3 The coordinate of the second line
 * @return The float value
 */
const float getBMetric(const std::vector<float>& first, 
					   const std::vector<float>& second);


/*
 * @brief Calculate the B-metric for two streamlines considered as multivariate normal distributions
 * @param[in] centerNormal The multivariate normal distribution of the first line
 * @param[in] neighNormal The multivriate normal distribution of the second line
 * @return The float value
 */
const float getBMetric(const MultiVariate& first, 
					   const MultiVariate& second);


/*
 * @brief Calculate the norm 10 given two streamlines
 *
 * @param[in] centroid The coordinate of the first line
 * @param[in] size The size of coordinates
 * @param[in] index The index of the second line
 * @param[in] unitLength The unit length vector for all the streamlines
 * @return The float value of norm 10 between two integral curves
 */
const float getMetric_10(const VectorXf& centroid,
						 const int& size,
						 const int& index,
						 const std::vector<VectorXf>& unitLength);


/*
 * @brief Calculate the norm 10 given the coordinates of two streamlines
 *
 * @param[in] firstRow The coordinate of the first line
 * @param[in] size The size of coordinates
 * @param[in] secondRow The coordinate of the second line
 * @return The float value of norm 10 between two integral curves
 */
const float getMetric_10(const VectorXf& firstRow,
						 const int& size,
						 const VectorXf& secondRow);


/*
 * @brief Calculate the norm 10 given the indices of two streamlines
 *
 * @param[in] first The index of the first line
 * @param[in] second The index of the second line
 * @param[in] unitLength The unit length vector for all the streamlines
 * @return The float value of norm 10 between two integral curves
 */
const float getMetric_10(const int& first,
						 const int& second,
						 const std::vector<VectorXf>& unitLength);


/*
 * @brief Get the MCP distance for two lines
 *
 * @param[in] first The first line coordinate
 * @param[in] second The second line coordinate
 * @return A distance value between two lines
 */
const float getMetric_MOP(const VectorXf& first, const VectorXf& second);


/*
 * @brief Get the Hausdorff distance between two lines
 *
 * @param[in] first The first line coordinates
 * @param[in] second The second line coordinates
 * @return The distance value
 */
const float getMetric_Hausdorff(const VectorXf& first, const VectorXf& second);


/*
 * @brief Calculate the norm 0, 1, 2, 5, 8 and 11 for two given streamlines
 *
 * @param[in] centroid The coordinate of the centroid line
 * @param[in] r2 The coordinate of the second line
 * @param[in] index The index of the second line
 * @param[in] normOption The norm option
 * @param[in] pairwise The pairwise vector the first line
 * @param[in] objectNorm The pairwise vector of the second line
 * @return The float value of the similarity measure
 */
const float getNorm(const Eigen::VectorXf& centroid,
					const Eigen::VectorXf& r2,
					const int& index,
					const int& normOption,
					const std::vector<std::vector<float> >& object,
					const std::vector<std::vector<float> >& objectNorm);


/*
 * @brief Calculate the norm 0, 1, 2, 5, 8 and 11 for two given streamlines
 *
 * @param[in] centroid The coordinate of the centroid line
 * @param[in] r2 The coordinate of the second line
 * @param[in] firstIndex The index of the first line
 * @param[in] secondIndex The index of the second line
 * @param[in] normOption The norm option
 * @param[in] pairwise The pairwise vector the first line
 * @param[in] objectNorm The pairwise vector of the second line
 * @return The float value of the similarity measure
 */
const float getNorm(const VectorXf& centroid,
					const VectorXf& r2,
					const int& firstIndex,
					const int& secondIndex,
					const int& normOption,
					const std::vector<std::vector<float> >& object,
					const std::vector<std::vector<float> >& objectNorm);


/*
 * @brief Calculate the norm 0, 1, 2, 5, 8 and 11 for two given streamlines
 *
 * @param[in] r1 The coordinate of the first line
 * @param[in] r2 The coordinate of the second line
 * @param[in] normOption The norm option
 * @return The float value of the similarity measure
 */
const float getNorm(const VectorXf& r1, 
					const VectorXf& r2,
					const int& normOption);


/*
 * @brief Compute the signature-based similarity measure (14) for two integral lines given two elements
 * and their histogram
 *
 * @param[in] firstArray The coordinate of the first line
 * @param[in] secondArray The coordinate of the second line
 * @param[in] firstHist The histogram of signatures for the first line
 * @param[in] secondHist The histogram of signatures for the second line
 * @return The float value for the norm 14
 */
const float getSignatureMetric(const Eigen::VectorXf& firstArray,
							   const Eigen::VectorXf& secondArray,
							   const std::vector<float>& firstHist,
							   const std::vector<float>& secondHist);


/*
 * @brief Calculate the signature-based similarity measure (14) with a line and the centroid line
 *
 * @param[in] centroid The given centroid line coordinate
 * @param[in] first The coordinate of the first line
 * @param[in] firstHist The signature histogram of the first line
 * @return The float value of similarity measure 14
 */
const float getSignatureMetric(const Eigen::VectorXf& centroid,
							   const Eigen::VectorXf& first,
							   const std::vector<float>& firstHist);


/*
 * @brief Calculate the signature-based similarity measure (14) for two streamlines with given coordinates
 *
 * @param[in] first The coordinate of the first line
 * @param[in] second The coordinate of the second line
 * @return The float value of similarity measure 14
 */
const float getSignatureMetric(const Eigen::VectorXf& first,
							   const Eigen::VectorXf& second);


/*
 * @brief Compute the Procrutes distance for two lines with given coordinates
 * @details
 * 	Get adapted Procrustes distance. For example, if vec has 100 points, it will calculate mean of 94 points
 *
 * @param[in] first The coordinate of the first line
 * @param[in] second The coordinate of the second line
 * @return The distance value for Procrutes distance (15)
 */
const float getProcrustesMetric(const Eigen::VectorXf& first,
								const Eigen::VectorXf& second);


/*
 * @brief Compute the Procrutes distance with segments for two lines with given coordinates
 * @details
 * 	Get adapted Procrustes distance. For example, if vec has 100 points, it will calculate mean of 14 points
 *
 * @param[in] first The coordinate of the first line
 * @param[in] second The coordinate of the second line
 * @return The distance value for Procrutes distance (15)
 */
const float getProcrustesMetricSegment(const Eigen::VectorXf& first,
									   const Eigen::VectorXf& second);


/*
 * @brief Calculate the entropy-based metric (16) for two integral curves
 * @details
 * 	The metric formulation can be seen in the paper, An Illustrative Visualization Framework for 3D Vector Fields
 *
 * @param[in] firstEntropy The entropy vector for the first line
 * @param[in] secondEntropy The entropy vector for the second line
 * @return The distance value between two lines
 */
const float getEntropyMetric(const std::vector<float>& firstEntropy,
		                     const std::vector<float>& secondEntropy);


/*
 * @brief Calculate the entropy-based metric (16) for two lines, with one coordinate and one entropy vector
 * @details
 * 	The formula can be seen in the paper, An Illustrative Visualization Framework for 3D Vector Fields, and
 * 	the condition is to be given by one entropy values and another as coordinate vector
 *
 * @param[in] firstEntropy The entropy list of the first line
 * @param[in] array The coordinate of the second line
 * @return The distance value of similarity measure (16)
 */
const float getEntropyMetric(const std::vector<float>& firstEntropy,
		                     const Eigen::VectorXf& array);


/*
 * @brief Calculate the entropy-based metric (16) for two lines, with the coordinates of two lines
 * @details
 * 	The formula can be seen in the paper, An Illustrative Visualization Framework for 3D Vector Fields
 *
 * @param[in] first The coordinate of the first line
 * @param[in] second The coordinate of the second line
 * @return The distance value of similarity measure (16)
 */
const float getEntropyMetric(const Eigen::VectorXf& first,
		                     const Eigen::VectorXf& second);


/*
 * @brief Calculate the time-based MCP (17) for two pathlines
 *
 * @param[in] first The coordinate of the first line
 * @param[in] second The coordinate of the second line
 * @return The distance value between two pathlines
 */
const float getPathline_MCP(const Eigen::VectorXf& first,
        					const Eigen::VectorXf& second);


/*
 * @brief Get the similarity measures given the coordinates and norm option
 * @param[in] data The coordinate matrix
 * @param[in] first The first index i
 * @param[in] second The second index j
 * @param[in] normOption The norm option
 * @param[in] object The MetricPreparation class object for distance computation
 * @return The distance value between line i and j
 */
const float getDisimilarity(const MatrixXf& data,
							const int& first,
							const int& second,
							const int& normOption,
							const MetricPreparation& object);


/*
 * @brief Get the similarity measures given the coordinates and norm option
 *
 * @param[in] others The input line coordinate
 * @param[in] data The coordinate matrix
 * @param[in] index The index of another line i
 * @param[in] normOption The norm option
 * @param[in] object The MetricPreparation class object for distance computation
 * @return The distance value between line i and the input vector others
 */
const float getDisimilarity(const VectorXf& others,
							const MatrixXf& data,
							const int& index,
							const int& normOption,
							const MetricPreparation& object);


/*
 * @brief Get the similarity measures given the coordinates and norm option
 *
 * @param[in] first The first input line coordinate
 * @param[in] second The second input line coordinate
 * @param[in] firstIndex The index of line i
 * @param[in] secondIndex The index of line j
 * @param[in] normOption The norm option
 * @param[in] object The MetricPreparation class object for distance computation
 * @return The distance value between line i and the input vector others
 */
const float getDisimilarity(const VectorXf& first,
							const VectorXf& second,
							const int& firstIndex,
							const int& secondIndex,
							const int& normOption,
							const MetricPreparation& object);


/*
 * @brief Assign values to the distance matrix
 *
 * @param[in] data The coordinate matrix
 * @param[in] normOption The norm option
 * @param[in] object The MetricPreparation object
 */
void getDistanceMatrix(const MatrixXf& data,
				       const int& normOption,
					   const MetricPreparation& object);


/*
 * @brief Get the similarity measures given the coordinates and norm option
 *
 * @param[in] first The first input line coordinate
 * @param[in] second The second input line coordinate
 * @param[in] normOption The norm option
 * @param[in] object The MetricPreparation class object for distance computation
 * @return The distance value between two lines
 */
const float getDisimilarity(const VectorXf& first,
							const VectorXf& second,
							const int& normOption,
							const MetricPreparation& object);


/*
 * @brief Delete the pointer of distance matrix
 *
 * @param[in] Row The row of distance matrix
 */
void deleteDistanceMatrix(const int& Row);


/*
 * @brief Compute and assign the averaged rotation (discrete curvatures) on the cluster representatives
 * @param[in] streamline The cluster representative coordinates
 * @param[out] rotation The rotation vector to be updated
 */
const float getRotation(const std::vector<vector<float> >& streamline, std::vector<float>& rotation);


/*
 * @brief To store the cluster information into the storage file
 *
 * @param[in] storage The vector that contains the candidates for each cluster
 */
void generateGroups(const std::vector<std::vector<int> >& storage);


#endif
