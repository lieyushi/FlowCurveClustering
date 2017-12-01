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


class Query
{
private:
	int similarityOption;
	std::vector<int> interestedCurve;
	std::vector<float> rotation;
	Eigen::MatrixXf storage;

	void getInteresting(const Eigen::MatrixXf& data,
						const int& Row,
						const int& Column);

	void searchClosest(const int& target, 
					   const int& closestNumber, 
					   const int& normOption, 
					   std::vector<int>& neighbor);

public:
	Query();
	Query(const Eigen::MatrixXf& data,
		  const int& Row,
		  const int& Column);
	~Query();

	void getClosestInteresting(const int& normOption,
							   std::vector<StringQuery>& searchResult);

	void getClosestCurve(const int& normOption,
					   	 std::vector<StringQuery>& searchResult);

	bool interestedEmpty();
};

#endif