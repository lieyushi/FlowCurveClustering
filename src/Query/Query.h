/*
 * @brief The class to perform the distance-based query for the streamlines/pathlines
 * @author Lieyu Shi
 */


#ifndef _QUERY_H_
#define _QUERY_H_

#include <cmath>
#include <cstring>
#include <cassert>
#include "IOHandler.h"
#include "Distance.h"


#ifndef M_PI
	#define M_PI 3.141592653
#endif

using namespace std;


/*
 * @brief The struct with the distance and the index
 */
struct QueryDistance
{
	float distance;
	int index;
	QueryDistance()
	{}
	QueryDistance(const float& distance, const int& index): 
				  distance(distance), index(index){}
	bool operator<(const QueryDistance& others)
	{
		return distance < others.distance;
	}
};


/*
 * @brief The class to perform streamline query with a given distance option
 * @details
 * 	It firstly supports the search of the streamline candidates with the minimal distance w.r.t. the selected candidate
 * 	It also supports the selection of interesting candidate streamlinea according to rotation value calculation
 */
class Query
{
private:

	/*
	 * @brief The norm option
	 */
	int similarityOption;

	/*
	 * @brief The potentially interesting integral curve index
	 */
	std::vector<int> interestedCurve;

	/*
	 * @brief The rotation values for those interested integral curves
	 */
	std::vector<float> rotation;

	/*
	 * @brief The matrix coordinates of the streamlines
	 */
	Eigen::MatrixXf storage;

	/*
	 * @brief Get the interesting curves with rotation-related calculation
	 *
	 * @param[in] data The matrix coordinates of the streamlines
	 * @param[in] Row The row size
	 * @param[in] Column The column size
	 */
	void getInteresting(const Eigen::MatrixXf& data,
						const int& Row,
						const int& Column);

	/*
	 * @brief Get the closest candidates compared to the given target
	 *
	 * @param[in] target The given target index of streamline
	 * @param[in] closestNumber The number of closest curves for the search
	 * @param[in] normOption The norm option
	 * @param[out] neighbor The calculated indices of the neighboring integral curves
	 */
	void searchClosest(const int& target, 
					   const int& closestNumber, 
					   const int& normOption, 
					   std::vector<int>& neighbor);

public:

	/*
	 * @brief The default constructor for the class
	 */
	Query();

	/*
	 * @brief The constructor with the input parameters and assign the coordinates
	 *
	 * @param[in] data The matrix coordinates as input
	 * @param[in] Row The row size
	 * @param[in] Column The column size
	 */
	Query(const Eigen::MatrixXf& data,
		  const int& Row,
		  const int& Column);

	/*
	 * @brief The destructor
	 */
	~Query();

	/*
	 * @brief Get the results of closest search (distance value and index)
	 *
	 * @param[in] normOption The norm option
	 * @param[out] searchResult The result of the streamline query
	 */
	void getClosestInteresting(const int& normOption,
							   std::vector<StringQuery>& searchResult);

	/*
	 * @brief Calculate the curve indices that are closest to the target candidate streamline
	 *
	 * @param[in] normOption The norm option
	 * @param[out] searchResult The result of the streamline query
	 */
	void getClosestCurve(const int& normOption,
					   	 std::vector<StringQuery>& searchResult);

	/*
	 * @brief Whether the interested candidate vector is empty or not
	 */
	bool interestedEmpty();
};

#endif
