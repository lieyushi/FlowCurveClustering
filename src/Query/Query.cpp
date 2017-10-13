#include "Query.h"


void Query::getInteresting(const Eigen::MatrixXf& data,
						   const int& Row,
						   const int& Column)
{
	rotation = std::vector<float>(Row);
#pragma omp parallel for schedule(dynamic) num_threads(8)
	for (int i = 0; i < Row; ++i)
	{
		float accumulation = 0.0, leftNorm, rightNorm, dotValue, result;
		const Eigen::VectorXf& array = data.row(i);
		const int& size = Column/3-2;
		Vector3f left, right;
		for (int j = 0; j < size; ++j)
		{
			//std::cout << array[j*3+3] << " " << array[j*3+4] << " " << array[j*3+5] << std::endl;
			left << array(j*3+3)-array(j*3), 
					array(j*3+4)-array(j*3+1), 
					array(j*3+5)-array(j*3+2);
			right << array(j*3+6)-array(j*3+3), 
					 array(j*3+7)-array(j*3+4), 
					 array(j*3+8)-array(j*3+5);
			dotValue = left.dot(right);
			leftNorm = left.norm();
			rightNorm = right.norm();
			if(leftNorm >= 1.0e-8 && rightNorm >=1.0e-8)
			{
				result = dotValue/leftNorm/rightNorm;
				result = min(1.0,(double)result);
				result = max(-1.0,(double)result);
				accumulation += acos(result);
			}
		}
		rotation[i] = accumulation;
	}

	for (int i = 0; i < Row; ++i)
	{
		if(rotation[i]>4.0*M_PI)
			interestedCurve.push_back(i);
	}

	std::cout << "Found " << interestedCurve.size() << " interesting streamlines!" << std::endl;
}


Query::Query()
{
	//default constructor
}
	

Query::Query(const Eigen::MatrixXf& data,
		     const int& Row,
		     const int& Column)
{
	storage = data;
	getInteresting(data, Row, Column);	
}
	
Query::~Query()
{
	interestedCurve.clear();
}


bool Query::interestedEmpty()
{
	return interestedCurve.empty();
}


void Query::getClosestInteresting(const int& normOption,
					   			  std::vector<StringQuery>& searchResult)
{
	assert(!interestedCurve.empty());
	StringQuery tempObject;
	string option;
	int closestNumber;
	int index;
	do
	{
		std::cout << "Input index of interesting curves? " 
				  << "Range from [0, " << (interestedCurve.size()-1) << "]" << std::endl;
		std::cin >> index;
		assert(index >=0 && index < interestedCurve.size());
		std::cin.ignore();

		std::cout << "Input number of closest strings? :" << std::endl;
		std::cin >> closestNumber;
		std::cin.ignore();

		searchClosest(interestedCurve[index], closestNumber, 
					  normOption, tempObject.neighbor);

		tempObject.index = interestedCurve[index];
		searchResult.push_back(tempObject);

		std::cout << "Want to have more string query? Y:Yes, N:No." << std::endl;
		getline(cin, option);
	}while(option.compare(string("Y"))==0 || option.compare(string("y"))==0);
}


void Query::getClosestCurve(const int& normOption,
					   		std::vector<StringQuery>& searchResult)
{
	assert(!interestedCurve.empty());
	StringQuery tempObject;
	string option;
	int closestNumber;
	int index;
	do
	{
		std::cout << "Input index of curves? " 
				  << "Range from [0, " << (storage.rows()-1) << "]" << std::endl;
		std::cin >> index;
		assert(index >=0 && index < storage.rows());
		std::cin.ignore();

		std::cout << "Input number of closest strings? :" << std::endl;
		std::cin >> closestNumber;
		std::cin.ignore();

		searchClosest(index, closestNumber, 
					  normOption, tempObject.neighbor);

		tempObject.index = index;
		searchResult.push_back(tempObject);

		std::cout << "Want to have more string query? Y:Yes, N:No." << std::endl;
		getline(cin, option);
	}while(option.compare(string("Y"))==0 || option.compare(string("y"))==0);
}



void Query::searchClosest(const int& target, 
						  const int& closestNumber, 
					  	  const int& normOption, 
					  	  std::vector<int>& neighbor)
{
	const int& row = storage.rows();
	const int& lineNumber = storage.cols()/3-2;
	std::vector<QueryDistance> distance(row);
	neighbor = std::vector<int>(closestNumber);
#pragma omp parallel for schedule(dynamic) num_threads(8)
	for (int i = 0; i < row; ++i)
	{
		float dist;
		if(i==target)
		{
			distance[i].distance = 0;
			distance[i].index = target;
			continue;
		}
			
		switch(normOption)
		{
		case 2:
		case 5:
		case 8:
		case 11:
			dist = getNorm(storage.row(target), storage.row(i),normOption);
			break;

		case 3:
			dist = getBMetric_3(storage.row(target), lineNumber, storage.row(i));
			break;

		case 4:
			dist = abs(rotation[target]/lineNumber-rotation[i]/lineNumber);
			break;

		case 6:
			dist = getBMetric_6(storage.row(target), lineNumber, storage.row(i));
			break;

		case 7:
			dist = getBMetric_7(storage.row(target), lineNumber, storage.row(i));
			break;

		case 9:
			dist = getBMetric_9(storage.row(target), lineNumber, storage.row(i));
			break;

		case 10:
			dist = getMetric_10(storage.row(target), lineNumber, storage.row(i));
			break;

		case 12:
			dist = getMetric_MOP(storage.row(target), storage.row(i));
		}

		distance[i].distance = -dist;
		distance[i].index = i;
	}

	std::make_heap(distance.begin(), distance.end());
	std::cout << "Target: " << target << ", closest is: "
			  << distance.front().index << std::endl;
	assert(distance.front().index==target);
	std::pop_heap(distance.begin(), distance.end());
	distance.pop_back();
	for (int i = 0; i < closestNumber; ++i)
	{
		neighbor[i] = distance.front().index;
		std::cout << distance.front().index << " " << distance.front().distance << std::endl;
		std::pop_heap(distance.begin(), distance.end());
		distance.pop_back();
	}

	for (int i = 0; i < neighbor.size(); ++i)
	{
		std::cout << neighbor[i] << std::endl;
	}
}