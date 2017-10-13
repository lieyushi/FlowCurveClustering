#include "Query.h"
#include <sys/time.h>

using namespace std;


void streamlineQuery(const int& argc,
					 char **argv);

int initializationOption;

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
	IOHandler::expandArray(data, dataVec, dimension, maxElements); //directly filling
	//IOHandler::sampleArray(data, dataVec, dimension, maxElements); //uniform sampling

	Query q = Query(data, dataVec.size(), maxElements);
	std::vector<StringQuery> searchResult;

	int searchInteresting;
	char isContinued;
	for (int i = 2; i < 12; ++i)
	{
		std::cout << "Wanted to continue the query?" << std::endl
				  << "Y. Yes; N. No" << std::endl;
		std::cin >> isContinued;
		if(isContinued=='N'||isContinued=='n')
			break;
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
