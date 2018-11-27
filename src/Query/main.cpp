#include "Query.h"
#include <sys/time.h>

using namespace std;


void streamlineQuery(const int& argc,
					 char **argv);

int initializationOption;
bool isPathlines;

int main(int argc, char* argv[])
{
	streamlineQuery(argc, argv);
	return 0;
}


void streamlineQuery(const int& argc,
					 char **argv)
{
	while(argc!=3)
	{
		std::cout << "Input argument should have 3!" << endl
		          << "./cluster inputFile data_dimension" << endl;
		exit(1);
	}
	const string& strName = string("../dataset/")+string(argv[1]);
	const int& dimension = atoi(argv[2]);
	std::vector< std::vector<float> > dataVec;
	int vertexCount, maxElements;
	IOHandler::readFile(strName, dataVec, vertexCount, dimension, maxElements);

	/*stringstream ss;
	ss << strName << "_differentNorm_full.vtk";
	const string& fullName = ss.str();
	IOHandler::printVTK(ss.str(), dataVec, vertexCount, dimension);
	ss.str("");*/
	
	Eigen::MatrixXf data;
	int userInput;
	std::cout << "It is pathlines? 1.Yes, 0.No" << std::endl;
	std::cin >> userInput;
	assert(userInput==1||userInput==0);
	isPathlines = (isPathlines==1);

	/* select sampling strategy, and 2 is often for geometric clustering */
	std::cout << "Please choose sampling strategy: " << std::endl
			  << "1.directly filling, 2.uniformly sampling" << std::endl;
	int samplingMethod;
	std::cin >> samplingMethod;
	assert(samplingMethod==1 || samplingMethod==2);
	if(samplingMethod==1)
		IOHandler::expandArray(data, dataVec, dimension, maxElements); //directly filling
	else if(samplingMethod==2)
		IOHandler::sampleArray(data, dataVec, dimension, maxElements); //uniform sampling

	Query q = Query(data, dataVec.size(), maxElements);
	std::vector<StringQuery> searchResult;

	int searchInteresting;
	char isContinued;

	/*  0: Euclidean Norm
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
	*/

	if(isPathlines)
	{
		for (int i = 2; i < 17; ++i)
		{
			if(i!=12)
				continue;

			if(i==17)
			{
				IOHandler::expandArray(data, dataVec, dimension, maxElements); //directly filling
				q = Query(data, dataVec.size(), maxElements);
			}
			//std::cout << "Wanted to continue the query?" << std::endl
			//		  << "Y. Yes; N. No" << std::endl;
			//std::cin >> isContinued;
			//if(isContinued=='N'||isContinued=='n')
			//	break;
			std::cin.ignore();
			std::cout << "-----------------------------------------------" << std::endl
					  << "--------------------- Norm " << i << " string query"
					  << "-----------------" << std::endl;

			std::cout << "Want to search among interesting curves or not?"
					  << " 0.No, 1.Yes!" << std::endl;
			std::cin >> searchInteresting;
			std::cin.ignore();

			if(searchInteresting==0)
				q.getClosestCurve(i,searchResult);
			else if(searchInteresting==1&&!q.interestedEmpty())
				q.getClosestInteresting(i,searchResult);

			for (int j = 0; j < searchResult.size(); ++j)
			{
				IOHandler::printQuery(i,j,searchResult[j], dataVec);
			}
			searchResult.clear();
			std::cout << "--------------------------------------------------" << std::endl;
		}
	}
	else
	{
		for (int i = 2; i < 15; ++i)
		{
			if(i!=12)
				continue;
			//std::cout << "Wanted to continue the query?" << std::endl
			//		  << "Y. Yes; N. No" << std::endl;
			//std::cin >> isContinued;
			//if(isContinued=='N'||isContinued=='n')
			//	break;
			std::cin.ignore();
			std::cout << "-----------------------------------------------" << std::endl
					  << "--------------------- Norm " << i << " string query"
					  << "-----------------" << std::endl;

			std::cout << "Want to search among interesting curves or not?"
					  << " 0.No, 1.Yes!" << std::endl;
			std::cin >> searchInteresting;
			std::cin.ignore();

			if(searchInteresting==0)
				q.getClosestCurve(i,searchResult);
			else if(searchInteresting==1&&!q.interestedEmpty())
				q.getClosestInteresting(i,searchResult);

			for (int j = 0; j < searchResult.size(); ++j)
			{
				IOHandler::printQuery(i,j,searchResult[j], dataVec);
			}
			searchResult.clear();
			std::cout << "--------------------------------------------------" << std::endl;
		}
	}
}
