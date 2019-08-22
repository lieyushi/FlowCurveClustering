/*
 * @brief This class contains the I/O functions for the program
 * @author Lieyu Shi
 */


#include "IOHandler.h"


/*
 * @brief Read the data from the file into the vector<vector<>>
 *
 * @param[in] fileName The file name of the data set
 * @param[out] dataVec The streamline coordinates to be updated
 * @param[out] vertexCount The total vertex count
 * @param[in] dimension 2 or 3 indicates it's 2D or 3D point
 * @param[out] maxElement The max dimension of streamlines
 */
void IOHandler::readFile(const string& fileName, 
						 std::vector< std::vector<float > >& dataVec, 
						 int& vertexCount, 
						 const int& dimension,
						 int& maxElement)
{

	vertexCount = 0;
	std::ifstream fin(fileName.c_str(), ios::in);
	if(!fin)
	{
		std::cout << "Error creating files!" << std::endl;
		exit(1);
	}
	stringstream ss;
	std::vector<float> tempVec;

	string line, part;

	/* read partial number of streamlines */
	//int MAXNUMBER;
	//std::cout << "Input maximal trajectory numbers: " << std::endl;
	//std::cin >> MAXNUMBER; 
	// set currentNumber to record how many streamlines u want to read in
	//int currentNumber = 0;


	/* read partial dimensions of curves */
	//int MAXDIMENSION;
	//std::cout << "Input maximal dimensions: " << std::endl;
	//std::cin >> MAXDIMENSION; 
	// set currentNumber to record how many streamlines u want to read in
	//int currentDimensions;

	std::vector<float> vec(3);
	float temp;
	maxElement = 0;
	while(getline(fin, line) /* && currentNumber < MAXNUMBER*/)
	{
		//currentDimensions = 0;
		int tag = 0, count = 0;
		ss.str(line);
		while(ss>>part /*&& currentDimensions<3*MAXDIMENSION*/)
		{
			/* operations below would remove duplicate vertices because that would damage our computation */
			temp = atof(part.c_str());
			if(tag>=3)
			{
				if(count<3)
				{
					vec[count] = temp;
					++tag;
					++count;
				}
				if(count==3)
				{
					int size = tempVec.size();
					//if(!(abs(vec[0]-tempVec[size-3])<1.0e-5&&abs(vec[1]-tempVec[size-2])<1.0e-5&&abs(vec[2]-tempVec.back())<1.0e-5))
					//{
					tempVec.push_back(vec[0]);
					tempVec.push_back(vec[1]);
					tempVec.push_back(vec[2]);
					//}
					count = 0;
				}
				continue;
			}
			tempVec.push_back(temp);
			++tag;
			//currentDimensions++;
		}
		/* accept only streamlines with at least three vertices */
		if(tempVec.size()/3>2)
		{
			if(maxElement<tempVec.size())
				maxElement = tempVec.size();
			dataVec.push_back(tempVec);
			vertexCount+=tempVec.size();
		}
		tempVec.clear();
		ss.clear();
		ss.str("");
		//currentNumber++;
	}
	fin.close();
	
	
	vertexCount/=dimension;
	std::cout << "File reader has been completed, and it toally has " << dataVec.size() << " trajectories and " 
	          << vertexCount << " vertices!" << std::endl;
	std::cout << "Max dimension is " << maxElement << std::endl;
}


/*
 * @brief Read the particle-based data from the given frame number
 *
 * @param[in] fileName The file of the particle based data
 * @param[out] dataVec The streamline coordinates
 * @param[in] vertexCount The total point count
 * @param[in] dimension 2 (2d) or 3 (3d) points
 * @param[in] trajectoryNum Number of trajectories
 * @param[in] Frame Number of frames
 */
void IOHandler::readFile(const string& fileName, 
						 std::vector< std::vector<float > >& dataVec, 
						 int& vertexCount, 
						 const int& dimension,
						 const int& trajectoryNum, 
						 const int& Frame)
{
	vertexCount = trajectoryNum*(Frame-1);
	dataVec = std::vector< std::vector<float> >(trajectoryNum, std::vector<float> ((Frame-1)*dimension));
#pragma omp parallel for schedule(static) num_threads(8)
	/* from 1 to Frame-1 then pay attention to i index */
	for (int i = 1; i < Frame; ++i)
	{
		stringstream ss;
		ss << fileName << i << ".txt";
		std::ifstream fin(ss.str().c_str(), ios::in);
		if(!fin)
		{
			std::cout << "File doesn't exist for this number!" << std::endl;
			exit(1);
		}
		float firstFloat;
		string line, linePart;

		ss.clear();
		ss.str("");
		for (int j = 0; j < trajectoryNum; ++j)
		{
			getline(fin, line);

			assert(!line.empty());

			ss.str(line);
			ss >> linePart;

			ss >> linePart;
			dataVec[j][(i-1)*dimension] = atof(linePart.c_str());

			ss >> linePart;
			dataVec[j][(i-1)*dimension+1] = atof(linePart.c_str());

			ss >> linePart;
			dataVec[j][(i-1)*dimension+2] = atof(linePart.c_str());
		}

		fin.close();
		std::cout << "File " << i << " has been read in successfully!" << std::endl;
	}


}


/*
 * @brief Print the vtk for streamlines
 *
 * @param[in] fileName The given vtk file name and position
 * @param[in] dataVec The streamline coordinates
 * @param[in] vertexCount the count of vertices
 * @param[in] dimension 2 or 3 indicates 2D or 3D
 * @param[in] clusterNumber the labels for each streamline
 * @param[in] sCluster The float scalar value for each streamline
 */
void IOHandler::printVTK(const string& fileName, 
						 const std::vector< std::vector<float > >& dataVec, 
						 const int& vertexCount, 
						 const int& dimension,
						 const std::vector<int>& clusterNumber,
						 const std::vector<float>& sCluster)
{
	if(clusterNumber.empty()||sCluster.empty())
		return;
	std::ofstream fout(fileName.c_str(), ios::out);
	if(!fout)
	{
		std::cout << "Error creating a new file!" << std::endl;
		exit(1);
	}
	fout << "# vtk DataFile Version 3.0" << std::endl << "Bernard streamline" << std::endl
	     << "ASCII" << std::endl << "DATASET POLYDATA" << std::endl;
	fout << "POINTS " << vertexCount << " float" << std::endl;

	int subSize, arraySize;
	std::vector<float> tempVec;
	for (int i = 0; i < dataVec.size(); ++i)
	{
		tempVec = dataVec[i];
		subSize = tempVec.size()/dimension;
		for (int j = 0; j < subSize; ++j)
		{
			for (int k = 0; k < dimension; ++k)
			{
				fout << tempVec[j*dimension+k] << " ";
			}
			fout << endl;
		}
	}

	fout << "LINES " << dataVec.size() << " " << (vertexCount+dataVec.size()) << std::endl;

	subSize = 0;
	for (int i = 0; i < dataVec.size(); ++i)
	{
		arraySize = dataVec[i].size()/dimension;
		fout << arraySize << " ";
		for (int j = 0; j < arraySize; ++j)
		{
			fout << subSize+j << " ";
		}
		subSize+=arraySize;
		fout << std::endl;
	}
	fout << "POINT_DATA" << " " << vertexCount << std::endl;
	fout << "SCALARS group int 1" << std::endl;
	fout << "LOOKUP_TABLE group_table" << std::endl;

	for (int i = 0; i < dataVec.size(); ++i)
	{
		arraySize = dataVec[i].size()/dimension;
		for (int j = 0; j < arraySize; ++j)
		{
			fout << clusterNumber[i] << std::endl;
		}
	}

	if(!sCluster.empty())
	{
		fout << "SCALARS sCluster float 1" << std::endl;
		fout << "LOOKUP_TABLE sCluster_table" << std::endl;

		for (int i = 0; i < dataVec.size(); ++i)
		{
			arraySize = dataVec[i].size()/dimension;
			for (int j = 0; j < arraySize; ++j)
			{
				fout << sCluster[clusterNumber[i]] << std::endl;
			}
		}
	}

	fout.close();
}


/*
 * @brief Print the vtk of the streamlines
 *
 * @param[in] fileName The name of the data set file
 * @param[in] dataVec The streamline coordinates
 * @param[in] vertexCount total count of vertices
 * @param[in] dimension 2 or 3
 */
void IOHandler::printVTK(const string& fileName, 
						 const std::vector< std::vector<float > >& dataVec, 
						 const int& vertexCount, 
						 const int& dimension)
{
	if(dataVec.empty())
		return;
	std::ifstream fin(fileName.c_str());
	if(fin.good())
		return;
	std::ofstream fout(fileName.c_str(), ios::out);
	if(!fout)
	{
		std::cout << "Error creating a new file!" << std::endl;
		exit(1);
	}
	fout << "# vtk DataFile Version 3.0" << std::endl << "Bernard streamline" << std::endl
	     << "ASCII" << std::endl << "DATASET POLYDATA" << std::endl;
	fout << "POINTS " << vertexCount << " float" << std::endl;

	int subSize, arraySize;
	std::vector<float> tempVec;
	for (int i = 0; i < dataVec.size(); ++i)
	{
		tempVec = dataVec[i];
		subSize = tempVec.size()/dimension;
		for (int j = 0; j < subSize; ++j)
		{
			for (int k = 0; k < dimension; ++k)
			{
				fout << tempVec[j*dimension+k] << " ";
			}
			fout << endl;
		}
	}

	fout << "LINES " << dataVec.size() << " " << (vertexCount+dataVec.size()) << std::endl;

	subSize = 0;
	for (int i = 0; i < dataVec.size(); ++i)
	{
		arraySize = dataVec[i].size()/dimension;
		fout << arraySize << " ";
		for (int j = 0; j < arraySize; ++j)
		{
			fout << subSize+j << " ";
		}
		subSize+=arraySize;
		fout << std::endl;
	}
	fout << "POINT_DATA" << " " << vertexCount << std::endl;
	fout << "SCALARS group int 1" << std::endl;
	fout << "LOOKUP_TABLE group_table" << std::endl;

	for (int i = 0; i < dataVec.size(); ++i)
	{
		arraySize = dataVec[i].size()/dimension;
		for (int j = 0; j < arraySize; ++j)
		{
			fout << i << std::endl;
		}
	}

	fout.close();
}


/*
 * @brief Print the vtk file of the streamlines
 *
 * @param[in] fileName the name of the data set
 * @param[in] dataVec The streamline coordinates
 * @param[in] vertexCount the count of vertices
 * @param[in] dimension 2 or 3
 * @param[in] sCluster the scalar value for streamlines
 */
void IOHandler::printVTK(const string& fileName, 
						 const std::vector<MeanLine>& dataVec, 
						 const int& vertexCount, 
						 const int& dimension,
						 const std::vector<float>& sCluster)
{
	if(dataVec.empty())
		return;
	std::ofstream fout(fileName.c_str(), ios::out);
	if(!fout)
	{
		std::cout << "Error creating a new file!" << std::endl;
		exit(1);
	}
	fout << "# vtk DataFile Version 3.0" << std::endl << "Bernard streamline" << std::endl
	     << "ASCII" << std::endl << "DATASET POLYDATA" << std::endl;
	fout << "POINTS " << vertexCount << " float" << std::endl;

	int subSize, arraySize;
	std::vector<float> tempVec;
	for (int i = 0; i < dataVec.size(); ++i)
	{
		tempVec = dataVec[i].minCenter;
		subSize = tempVec.size()/dimension;
		for (int j = 0; j < subSize; ++j)
		{
			for (int k = 0; k < dimension; ++k)
			{
				fout << tempVec[j*dimension+k] << " ";
			}
			fout << endl;
		}
	}

	fout << "LINES " << dataVec.size() << " " << (vertexCount+dataVec.size()) << std::endl;

	subSize = 0;
	for (int i = 0; i < dataVec.size(); ++i)
	{
		arraySize = dataVec[i].minCenter.size()/dimension;
		fout << arraySize << " ";
		for (int j = 0; j < arraySize; ++j)
		{
			fout << subSize+j << " ";
		}
		subSize+=arraySize;
		fout << std::endl;
	}
	fout << "POINT_DATA" << " " << vertexCount << std::endl;
	fout << "SCALARS group int 1" << std::endl;
	fout << "LOOKUP_TABLE group_table" << std::endl;

	for (int i = 0; i < dataVec.size(); ++i)
	{
		arraySize = dataVec[i].minCenter.size()/dimension;
		for (int j = 0; j < arraySize; ++j)
		{
			fout << dataVec[i].cluster << std::endl;
		}
	}

	if(!sCluster.empty())
	{
		fout << "SCALARS sCluster float 1" << std::endl;
		fout << "LOOKUP_TABLE sCluster_table" << std::endl;

		for (int i = 0; i < dataVec.size(); ++i)
		{
			arraySize = dataVec[i].minCenter.size()/dimension;
			for (int j = 0; j < arraySize; ++j)
			{
				fout << sCluster[dataVec[i].cluster] << std::endl;
			}
		}
	}

	fout.close();
}


/*
 * @brief Direct repeating the last point of the streamlines
 *
 * @param[out] data The matrix coordinates to be updated
 * @param[in] dataVec The line coordinates of streamlines as input
 * @param[in] dimension 2 or 3
 * @param[in] maxElement The max dimension of streamlines
 */
void IOHandler::expandArray(MatrixXf& data, 
							const std::vector< std::vector<float> >& dataVec, 
							const int& dimension, 
							const int& maxElements)
{
	data = Eigen::MatrixXf(dataVec.size(), maxElements);
#pragma omp parallel for schedule(static) num_threads(8)
	for (int i = 0; i < dataVec.size(); ++i)
	{
		const std::vector<float>& eachVec = dataVec[i];
		const int& vecSize = eachVec.size();
		//data.row(i) = Eigen::VectorXf::Map(&(eachVec[0]), vecSize);
		for (int j = 0; j<vecSize; j++)
			data(i,j) = eachVec[j];

		for (int j = vecSize; j < maxElements; j=j+dimension)
		{
			for (int k=0; k<dimension; k++)
				data(i,j+k) = eachVec[vecSize-dimension+k];
		}
	}
}


/*
 * @brief Sample the streamlines with equal arc length
 *
 * @param[out] data The matrix coordinates to be updated
 * @param[in] dataVec The input streamline coordinates
 * @param[in] dimension 2 or 3
 * @param[in] maxElements The max dimension of streamlines
 */
/* sample equal-sized array by equal arcs given a numOfVertices count */
void IOHandler::uniformArcSampling(MatrixXf& data,
								   const std::vector< std::vector<float> >& dataVec,
								   const int& dimension,
								   const int& maxElements)
{
	const int& numOfRows = dataVec.size();

	/* assign memory for required matrix */
	const int& totalSize = 3*maxElements;
	data = Eigen::MatrixXf(numOfRows, totalSize);
#pragma omp parallel for schedule(static) num_threads(8)
	/* compute total length of streamline and record each cumulative length information */
	for(int i=0;i<numOfRows;++i)
	{
		float entireLength = 0.0, lineLength;
		const std::vector<float>& eachVec = dataVec[i];
		const int& vecSize = eachVec.size();
		const int& lineNum = vecSize/3-1;
		Eigen::Vector3f lineSeg;

		vector<float> pairwise(lineNum);
		for(int j=0;j<lineNum;++j)
		{
			lineSeg(0)=eachVec[3*j+3]-eachVec[3*j];
			lineSeg(1)=eachVec[3*j+4]-eachVec[3*j+1];
			lineSeg(2)=eachVec[3*j+5]-eachVec[3*j+2];
			lineLength = lineSeg.norm();
			entireLength+=lineLength;
			pairwise[j]=entireLength;
		}

		float eachLength = entireLength/(maxElements-1);
		Eigen::VectorXf row_i(data.row(i).size());

		/* insert starting vertex */
		row_i(0)=eachVec[0];
		row_i(1)=eachVec[1];
		row_i(2)=eachVec[2];

		/* insert ending vertex */
		row_i(totalSize-3)=eachVec[vecSize-3];
		row_i(totalSize-2)=eachVec[vecSize-2];
		row_i(totalSize-1)=eachVec[vecSize-1];

		float tempLength, tempRatio, ratioComplement;
		int preLine, preVertex, postVertex;
		for(int j=1;j<=maxElements-2;++j)
		{
			/* current length */
			tempLength = j*eachLength;

			/* pre-index of coordinates */
			preLine = std::lower_bound(pairwise.begin(), pairwise.end(), tempLength)-pairwise.begin()-1;

			/* locate vertex index */
			preVertex = preLine+1;
			postVertex = preVertex+1;

			/* use linear interpolation to generate new coordinates */
			tempRatio = (tempLength-pairwise[preLine])/(pairwise[preLine+1]-pairwise[preLine]);
			ratioComplement = 1.0-tempRatio;
			row_i(3*j) = tempRatio*eachVec[3*postVertex]+ratioComplement*eachVec[3*preVertex];
			row_i(3*j+1) = tempRatio*eachVec[3*postVertex+1]+ratioComplement*eachVec[3*preVertex+1];
			row_i(3*j+2) = tempRatio*eachVec[3*postVertex+2]+ratioComplement*eachVec[3*preVertex+2];
		}

		data.row(i) = row_i;
	}
}


/*
 * @brief Sample the streamlines on the intervals to preservce thier geometry shapes
 *
 * @param[out] data The matrix coordinates to be updated
 * @param[in] dataVec The streamline coordinates as input
 * @param[in] dimension 2 or 3
 * @param[in] maxElements The maximal length of streamlines
 */
void IOHandler::sampleArray(MatrixXf& data, 
							const std::vector< std::vector<float> >& dataVec, 
							const int& dimension, 
							const int& maxElements)
{
	/*maxElements = INT_MIN;
	int arraySize;
	for (int i = 0; i < dataVec.size(); ++i)
	{
		arraySize = dataVec[i].size();
		if(maxElements < arraySize)
			maxElements = arraySize;
	}
	std::cout << maxElements << std::endl;*/
	//temp.row(i) = Eigen::VectorXf::Map(&each[0], 10); //must match the column size
	data = Eigen::MatrixXf(dataVec.size(), maxElements);
#pragma omp parallel for schedule(static) num_threads(8)
	for (int i = 0; i < dataVec.size(); ++i)
	{
		const std::vector<float>& eachVec = dataVec[i]; //cached vector<float>
		if(eachVec.size()==maxElements)
		{
			data.row(i) = Eigen::VectorXf::Map(&eachVec[0], maxElements);
		}
		else
		{
			const int& pointNum = eachVec.size()/3;	//current vec point length
			const int& totalNum = maxElements/3; //totally maximal point length
			const int& segNum = pointNum-1;
			const int& averageAdd = (totalNum-pointNum)/segNum; //average point on each segment

		//# of segments with averageAdd+1 sampled points 
			const int& averageRes = (totalNum-pointNum)%segNum; 

			int segmentLength;
			int currentPoint = 0;

			Eigen::Vector3f meanLength, insertedPoint;
			Eigen::Vector3f front, end;

			int j;
			for(j=0; j<segNum;j++)	//traverse all segments
			{
				for(int k=0; k<3;k++)
					data(i,3*currentPoint+k) = eachVec[3*j+k];
				currentPoint++;
				if(j<segNum-averageRes)
					segmentLength = averageAdd;
				else
					segmentLength = averageAdd+1;
				if(segmentLength>=1)
				{
					front << eachVec[3*j], eachVec[3*j+1], eachVec[3*j+2];
					end << eachVec[3*j+3], eachVec[3*j+4], eachVec[3*j+5];			
					meanLength = (end-front)/(segmentLength+1);
					for(int k=1; k<=segmentLength; k++)
					{
						insertedPoint = front+k*meanLength;
						for(int s=0; s<3;s++)
							data(i,3*currentPoint+s) = insertedPoint(s);
						currentPoint++;
					}
				}
			}
			assert(currentPoint==totalNum-1);
			for(int k=0; k<3;k++)
				data(i,3*currentPoint+k) = eachVec[3*j+k];
		}
	}
}


/*
 * @brief Assign vector values to the pointer type
 *
 * @param[out] data The pointer to be updated
 * @param[in] dataVec The streamline coordinates
 * @param[in] dimension 2 or 3
 */
void IOHandler::formArray(float ***data, 
						  const std::vector< std::vector<float> >& dataVec, 
						  const int& dimension)
{
	*data = new float*[dataVec.size()];
#pragma omp parallel for schedule(static) num_threads(8)
	for (int i = 0; i < dataVec.size(); ++i)
	{
		const int& arraySize = dataVec[i].size();
		(*data)[i] = new float[arraySize];
		memcpy(&(*data)[i][0], &(dataVec[i][0]), arraySize*sizeof(float));
	}
}


/*
 * @brief Print the new group labels to the existing vtk file
 *
 * @param[in] dataVec The streamline coordinates as input
 * @param[in] group The label vector for each streamline
 * @param[in] totalNum The size of each cluster just in case to show dominant clusters
 * @param[in] groupName The string text for the labels
 * @param[in] fullName The full name of the vtk to be written in
 * @param[in] dimension 2 or 3
 */
void IOHandler::printToFull(const std::vector< std::vector<float> >& dataVec, 
							const std::vector<int>& group, 
				 			const std::vector<int>& totalNum, 
				 			const string& groupName, 
				 			const string& fullName, 
				 			const int& dimension)
{
	if(group.empty()||totalNum.empty())
		return;
	std::ofstream fout(fullName.c_str(), ios::out | ios::app );
	if(!fout)
	{
		std::cout << "Error opening the file!" << std::endl;
		exit(1);
	}

	fout << "SCALARS " << groupName << " int 1" << std::endl;
	fout << "LOOKUP_TABLE " << groupName+string("_table") << std::endl;

	int arraySize;
	for (int i = 0; i < dataVec.size(); ++i)
	{
		arraySize = dataVec[i].size()/dimension;
		for (int j = 0; j < arraySize; ++j)
		{
			fout << group[i] << std::endl;
		}
	}

	fout << "SCALARS " <<  groupName + "_num" << " int 1" << std::endl;
	fout << "LOOKUP_TABLE " <<  groupName+string("_num_table") << std::endl;

	for (int i = 0; i < dataVec.size(); ++i)
	{
		arraySize = dataVec[i].size()/dimension;
		for (int j = 0; j < arraySize; ++j)
		{
			fout << totalNum[i] << std::endl;
		}
	}
	fout.close(); 
}


/*
 * @brief Print scalar values to the existing vtk file
 *
 * @param[in] dataVec The streamline data set
 * @param[in] sData The scalar values on the streamlines
 * @param[in] groupName The string text for the scalar values
 * @param[in] fullName The full name of the .vtk file
 * @param[in] dimension 2 or 3
 */
void IOHandler::printToFull(const std::vector< std::vector<float> >& dataVec, 
							const std::vector<float>& sData, 
				 			const string& groupName, 
				 			const string& fullName, 
				 			const int& dimension)
{
	if(sData.empty()||dataVec.empty())
		return;
	std::ofstream fout(fullName.c_str(), ios::out | ios::app );
	if(!fout)
	{
		std::cout << "Error opening the file!" << std::endl;
		exit(1);
	}

	if(!sData.empty())	
	{
		fout << "SCALARS " << groupName << " float 1" << std::endl;
		fout << "LOOKUP_TABLE " << groupName+string("_table") << std::endl;

		int arraySize;
		for (int i = 0; i < dataVec.size(); ++i)
		{
			arraySize = dataVec[i].size()/dimension;
			for (int j = 0; j < arraySize; ++j)
			{
				fout << sData[i] << std::endl;
			}
		}
	}

	fout.close(); 
}


/*
 * @brief Print scalar and label values to the existing vtk file
 *
 * @param[in] dataVec The streamline data set
 * @param[in] group The group labels for the streamlines
 * @param[in] sCluster The scalar values on the streamlines
 * @param[in] groupName The string text for the scalar values
 * @param[in] fullName The full name of the .vtk file
 * @param[in] dimension 2 or 3
 */
void IOHandler::printToFull(const std::vector< std::vector<float> >& dataVec, 
							const std::vector<int>& group, 
							const std::vector<float>& sCluster, 
				 			const string& groupName, 
				 			const string& fullName, 
				 			const int& dimension)
{
	if(dataVec.empty()||group.empty()||sCluster.empty())
		return;
	std::ofstream fout(fullName.c_str(), ios::out | ios::app );
	if(!fout)
	{
		std::cout << "Error opening the file!" << std::endl;
		exit(1);
	}

	if(!sCluster.empty())
	{
		fout << "SCALARS " <<  groupName << " float 1" << std::endl;
		fout << "LOOKUP_TABLE " <<  groupName+string("_table") << std::endl;

		int arraySize;
		for (int i = 0; i < dataVec.size(); ++i)
		{
			arraySize = dataVec[i].size()/dimension;
			for (int j = 0; j < arraySize; ++j)
			{
				if(group[i]<0)
					fout << 0 << std::endl;
				else
					fout << sCluster[group[i]] << std::endl;
			}
		}
	}
	fout.close(); 
}


/*
 * @brief delete the array
 *
 * @param[out] data The pointer type
 * @param[in] row The number of rows
 */
void IOHandler::deleteArray(float **data, 
							const int& row)
{
	if(data==NULL)
		return;
#pragma omp parallel for schedule(static) num_threads(8)
	for (int i = 0; i < row; ++i)
	{
		delete[] data[i];
	}
	delete[] data;
}


/*
 * @brief Write two double values in the readme
 *
 * @param[in] PCA_KMeans_delta The time for PCA k-means
 * @param[in] KMeans_delta The time for k-means
 */
void IOHandler::writeReadme(const double& PCA_KMeans_delta, 
							const double& KMeans_delta)
{
	std::ofstream readme("../dataset/README",ios::out | ios::app);
	if(!readme)
	{
		std::cout << "Error creating readme!" << std::endl;
		exit(1);
	}
	readme << "PCA_KMeans time elapse is " << PCA_KMeans_delta << " s." << std::endl
		   << "KMeans time elapse is " << KMeans_delta << " s." << std::endl;
    readme << std::endl;
    readme.close();
}


/*
 * @brief Write string and float vector in the readme
 *
 * @param[in] comment The string to be put in the readme
 * @param[in] sAverage The float vector to be written
 */
void IOHandler::writeReadme(const string& comment,
							const std::vector<float>& sAverage)
{
	if(sAverage.empty())
		return;
	std::ofstream readme("../dataset/README",ios::out | ios::app);
	if(!readme)
	{
		std::cout << "Error creating readme!" << std::endl;
		exit(1);
	}
	readme << comment << std::endl;
	for (int i = 0; i < sAverage.size(); ++i)
	{
		readme << sAverage[i] << std::endl;
	}
    readme << std::endl;
    readme.close();
}


/*
 * @brief Write the double and string vector and cluster in the readme
 *
 * @param[in] timeName The string vector
 * @param[in] timeDiff The double vector
 * @param[in] cluster The number of clusters
 */
void IOHandler::writeReadme(const std::vector<string>& timeName, 
							const std::vector<double>& timeDiff,
							const int& cluster)
{
	if(timeName.empty()||timeDiff.empty())
		return;
	std::ofstream readme("../dataset/README",ios::out | ios::app);
	if(!readme)
	{
		std::cout << "Error creating readme!" << std::endl;
		exit(1);
	}
	assert(timeName.size()==timeDiff.size());
	for (int i = 0; i < timeName.size(); ++i)
	{
		readme << timeName[i] << " is " << timeDiff[i] << " s." << std::endl;
	}
	readme << std::endl;
	readme << "Preset cluster number in K-means is: " << cluster << std::endl;
	readme << std::endl;
    readme.close();
}


/*
 * @brief Write two string vectors in the readme for printing information
 *
 * @param[in] timeName The string vector
 * @param[in] timeDiff The string vector
 * @param[in] cluster The size of cluster
 */
void IOHandler::writeReadme(const std::vector<string>& timeName,
							const std::vector<string>& timeDiff,
							const int& cluster)
{
	if(timeName.empty()||timeDiff.empty())
		return;
	std::ofstream readme("../dataset/README",ios::out | ios::app);
	if(!readme)
	{
		std::cout << "Error creating readme!" << std::endl;
		exit(1);
	}
	assert(timeName.size()==timeDiff.size());
	for (int i = 0; i < timeName.size(); ++i)
	{
		readme << timeName[i] << " " << timeDiff[i] << std::endl;
	}
	readme << "Preset cluster number in K-means is: " << cluster << std::endl;
	readme << std::endl;
    readme.close();
}


/*
 * @brief Write the extracted line information w.r.t. the norm in the readme
 *
 * @param[in] closest The closest extracted lines in the readme
 * @param[in] furthest The furthest extracted lines in the readme
 * @param[in] normOption The norm option
 */
void IOHandler::writeReadme(const std::vector<ExtractedLine>& closest, 
							const std::vector<ExtractedLine>& furthest, 
							const int& normOption)
{
	if(closest.empty()||furthest.empty())
		return;
	std::ofstream readme("../dataset/README",ios::out | ios::app);
	if(!readme)
	{
		std::cout << "Error creating readme!" << std::endl;
		exit(1);
	}
	const string& normStr = "Norm_"+to_string(normOption);
	readme << std::endl;
	readme << normStr+ " closest streamline set has " << closest.size() << " streamlines" << std::endl;
	for (int i = 0; i < closest.size(); ++i)
	{
		readme << closest[i].lineNum << " ";
	}
	readme << std::endl;
	
	readme << std::endl;
	readme << normStr+ " furthest streamline set has " << furthest.size() << " streamlines" << std::endl;
	for (int i = 0; i < furthest.size(); ++i)
	{
		readme << furthest[i].lineNum << " ";
	}
	readme << std::endl;
    readme.close();
}


/*
 * @brief Write closest and furthest extracted lines in the readme
 *
 * @param[in] closest The closest extracted lines
 * @param[in] furthest The furthest extracted lines
 */
void IOHandler::writeReadme(const std::vector<ExtractedLine>& closest, 
							const std::vector<ExtractedLine>& furthest)
{
	if(closest.empty()||furthest.empty())
		return;
	std::ofstream readme("../dataset/README",ios::out | ios::app);
	if(!readme)
	{
		std::cout << "Error creating readme!" << std::endl;
		exit(1);
	}
	readme << std::endl;
	readme << "PCA closest streamline set has " << closest.size() << " streamlines" << std::endl;
	for (int i = 0; i < closest.size(); ++i)
	{
		readme << closest[i].lineNum << " ";
	}
	readme << std::endl;
	
	readme << std::endl;
	readme << "PCA furthest streamline set has " << furthest.size() << " streamlines" << std::endl;
	for (int i = 0; i < furthest.size(); ++i)
	{
		readme << furthest[i].lineNum << " ";
	}
	readme << std::endl;
    readme.close();
}


/*
 * @brief Write the entropy, silhoutte and string in the readme
 *
 * @param[in] entropy The entropy value
 * @param[in] sil The Silhouette class object
 * @param[in] norm_str The string object
 */
void IOHandler::writeReadme(const float& entropy, const Silhouette& sil, const string& norm_str)
{
	std::ofstream readme("../dataset/README",ios::out | ios::app);
	if(!readme)
	{
		std::cout << "Error creating readme!" << std::endl;
		exit(1);
	}
	readme << norm_str << std::endl;
	readme << "The average silhouette: " << sil.sAverage
		   << ", the gamma statistic is: " << sil.gammaStatistic
		   << ", the entropy is: " << entropy
		   << ", the DB index is: " << sil.dbIndex
		   << std::endl;
}


/*
 * @brief Write two float values in the readme
 *
 * @param[in] closestAverage The float to be written
 * @param[in] furthestAverage The float to be written
 */
void IOHandler::writeReadme(const float& closestAverage, const float& furthestAverage)
{
	std::ofstream readme("../dataset/README",ios::out | ios::app);
	if(!readme)
	{
		std::cout << "Error creating readme!" << std::endl;
		exit(1);
	}
	readme << "The average rotation of closest is: " << closestAverage
		   << ", of furthest is: " << furthestAverage << std::endl;
}


/*
 * @brief Assign the ExtractedLine information to the closest streamline coordinates and cluster
 *
 * @param[out] closestStreamline The closest extracted streamlines
 * @param[out] cluster The number of clusters
 * @param[in] closest The closest ExtractedLine object as input
 * @param[out] pointNumber The number of points
 * @param[in] dataVec The input coordinates of the streamlines
 */
void IOHandler::assignVec(std::vector<std::vector<float> >& closestStreamline, 
						  std::vector<int>& cluster,
						  const std::vector<ExtractedLine>& closest,
						  int& pointNumber, 
						  const std::vector< std::vector<float> >& dataVec)
{
	if(closest.empty())
		return;
	closestStreamline = std::vector<std::vector<float> >(closest.size(), std::vector<float>());
	cluster = std::vector<int>(closest.size());
	pointNumber = 0;
	for (int i = 0; i < closestStreamline.size(); ++i)
	{
		closestStreamline[i] = dataVec[closest[i].lineNum];
		pointNumber+=closestStreamline[i].size();
		cluster[i] = closest[i].cluster;
	}	
}


/*
 * @brief Assign the size of clustes to the center
 *
 * @param[out] cluster The size of clusters
 * @param[in] centerMass The extracted mean line
 */
void IOHandler::assignVec(std::vector<int>& cluster,
						  const std::vector<MeanLine>& centerMass)
{
	cluster = std::vector<int>(centerMass.size());
#pragma omp parallel for schedule(static) num_threads(8)
	for (int i = 0; i < cluster.size(); ++i)
	{
		cluster[i] = centerMass[i].cluster;
	}
}


/*
 * @brief Write group information into the local file as backup
 *
 * @param[in] The labels of group for the streamlines
 * @param[in] The streamline coordinates
 */
void IOHandler::writeGroup(const std::vector<int>& group, 
						   const std::vector< std::vector<float> >& dataVec)
{
	if(group.empty()||dataVec.empty())
		return;
	std::ofstream readme("../dataset/group",ios::out);
	if(!readme)
	{
		std::cout << "Error creating readme!" << std::endl;
		exit(1);
	}
	assert(group.size()==dataVec.size());
	for (int i = 0; i < group.size(); ++i)
	{
		readme << group[i] << std::endl;
	}
    readme.close();
}


/*
 * @brief Print the query information from the streamline query
 *
 * @param[in] normOption The norm option
 * @param[in] order The int of order
 * @param[in] queryResult The string query result
 * @param[in] dataVec The streamline coordinates
 */
void IOHandler::printQuery(const int& normOption,
					   	   const int& order,
					   	   const StringQuery& queryResult, 
					  	   const std::vector<std::vector<float> >& dataVec)
{
	stringstream ss;
	ss << "../dataset/norm" << normOption << "_query" << order << "_target.vtk";
	const string& targetStr = ss.str();
	ss.str("");
	ss.clear();

	ss << "../dataset/norm" << normOption << "_query" << order << "_result.vtk";
	const string& resultStr = ss.str();

	/* print out the target streamline vtk file */
	ofstream fTarget(targetStr.c_str(),ios::out);
	if(!fTarget)
	{
		std::cout << "Error creating file!" << std::endl;
		exit(1);
	}
	std::cout << queryResult.index << std::endl;
	std::vector<float> targetVec = dataVec[queryResult.index];
	int pointNumber = targetVec.size()/3;

	fTarget << "# vtk DataFile Version 3.0" << std::endl << "Bernard streamline" << std::endl
	        << "ASCII" << std::endl << "DATASET POLYDATA" << std::endl;
	fTarget << "POINTS " << pointNumber << " float" << std::endl;

	for (int i = 0; i < pointNumber; ++i)
	{
		fTarget << targetVec[i*3] << " " << targetVec[i*3+1]
		        << " " << targetVec[i*3+2] << std::endl;
	}

	fTarget << "LINES " << 1 << " " << (1+pointNumber) << std::endl;
	fTarget << pointNumber << " ";
	for (int i = 0; i < pointNumber; ++i)
	{
		fTarget << i << " ";
	}
	fTarget << std::endl;
	fTarget.close();


	/* print out the streamline query result vtk file */
	ofstream fResult(resultStr.c_str(),ios::out);
	if(!fResult)
	{
		std::cout << "Error creating file!" << std::endl;
		exit(1);
	}
	pointNumber = 0;
	const std::vector<int>& neighbor = queryResult.neighbor;
	for (int i = 0; i < neighbor.size(); ++i)
	{
		pointNumber += dataVec[neighbor[i]].size()/3;
	}

	fResult << "# vtk DataFile Version 3.0" << std::endl << "Bernard streamline" << std::endl
	        << "ASCII" << std::endl << "DATASET POLYDATA" << std::endl;
	fResult << "POINTS " << pointNumber << " float" << std::endl;

	int subArraySize, indexNumber = 0;
	std::vector<float> tempVec;
	for (int i = 0; i < neighbor.size(); ++i)
	{
		tempVec = dataVec[neighbor[i]];
		subArraySize = tempVec.size()/3;
		for (int j = 0; j < subArraySize; ++j)
		{
			fResult << tempVec[3*j] << " " << tempVec[3*j+1] << " "
					<< tempVec[3*j+2] << std::endl;
		}
	}

	fResult << "LINES " << neighbor.size() << " " 
			<< (neighbor.size()+pointNumber) << std::endl;
	for (int i = 0; i < neighbor.size(); ++i)
	{
		tempVec = dataVec[neighbor[i]];
		subArraySize = tempVec.size()/3;
		fResult << subArraySize << " ";
		for (int j = 0; j < subArraySize; ++j)
		{
			fResult << (indexNumber+j) << " ";
		}
		fResult << std::endl;
		indexNumber += subArraySize;
	}
	fResult.close();

}


/*
 * @brief print the coordinates of the streamlines into the local txt file
 *
 * @param[in] The float pointer
 * @param[in] Row The size of rows
 * @param[in] Column The size of columns
 */
void IOHandler::printTXT(float **data,
						 const int& Row,
						 const int& Column)
{
	std::ofstream fout("../dataset/full.txt", ios::out);
	if(!fout)
	{
		std::cout << "Error creating a file!" << std::endl;
		exit(1);
	}
	float *array;
	for (int i = 0; i < Row; ++i)
	{
		array = data[i];
		for (int j = 0; j < Column; ++j)
		{
			fout << array[j] << " ";
		}
		fout << std::endl;
	}
	fout.close();
}


/*
 * @brief Directly align vector object to Eigen::MatrixXf. Map vector to vector as well, or called direct repeating
 *
 * @param[out] equalArray The matrix to be assigned the coordinates
 * @param[in] trajectories The trajectory coordinates
 * @param[in] dimension The dimension
 * @param[in] maxElement The max element
 */
void IOHandler::expandArray(std::vector<std::vector<float> >& equalArray,
							const std::vector<std::vector<float> >& trajectories, 
						 	const int& dimension,
						    const int& maxElement)
{
	equalArray = std::vector<std::vector<float> >(trajectories.size(),
				 std::vector<float>(maxElement));
#pragma omp parallel for schedule(static) num_threads(8)
	for (int i = 0; i < trajectories.size(); ++i)
	{
		std::vector<float>& tempRow = equalArray[i];
		const std::vector<float>& tempTraj = trajectories[i];
		const int& vecSize = tempTraj.size();
		memcpy(&(tempRow[0]), &(tempTraj[0]), vecSize*sizeof(float));
		for (int j = vecSize; j < maxElement; j=j+dimension)
		{
			memcpy(&(tempRow[j]), &(tempTraj[vecSize-dimension]), 
				   dimension*sizeof(float));
		}
	}
}


/*
 * @brief print the group information to the existing full vtk file with label names
 *
 * @param[in] dataVec The streamline coordinates
 * @param[in] group The labels of the streamlines
 * @param[in] fullName The full name of the primary vtk file
 * @param[in] groupName The group name of the labels
 * @param[in] dimension The dimension of the points
 */
void IOHandler::printToFull(const std::vector< std::vector<float> >& dataVec,
				 			const std::vector<int>& group, 
				 			const string& fullName,
				 			const string& groupName,
				 			const int& dimension)
{
	std::ofstream fout(fullName.c_str(), ios::out | ios::app );
	if(!fout)
	{
		std::cout << "Error opening the file!" << std::endl;
		exit(1);
	}

	fout << "SCALARS " << groupName << " int 1" << std::endl;
	fout << "LOOKUP_TABLE " << groupName+string("_table") << std::endl;

	int arraySize;
	for (int i = 0; i < dataVec.size(); ++i)
	{
		arraySize = dataVec[i].size()/dimension;
		for (int j = 0; j < arraySize; ++j)
		{
			fout << group[i] << std::endl;
		}
	}
	fout.close(); 
}


/*
 * @brief print the streamlines with attached scalar values in the vtk format
 *
 * @param[in] fileName The file name of the vtk
 * @param[in] array The streamline coordinates
 * @param[in] sCluster The float vector as input
 * @param[in] dimension The dimension of the points
 */
void IOHandler::printFeature(const string& fileName,
				  			 const std::vector<std::vector<float> >& array,
				  			 const std::vector<float>& sCluster,
				  			 const int& dimension)
{
	if(array.empty() || sCluster.empty())
		return;
	stringstream ss;
	ss << "../dataset/" << fileName;
	ofstream fout(ss.str().c_str(), ios::out);
	if(!fout)
	{
		std::cout << "Error creating file!" << std::endl;
		exit(-1);
	}

	int vertexCount = 0;
	for (int i = 0; i < array.size(); ++i)
	{
		vertexCount += array[i].size();
	}
	vertexCount /= dimension;

	fout << "# vtk DataFile Version 3.0" << std::endl << "Bernard streamline" << std::endl
	     << "ASCII" << std::endl << "DATASET POLYDATA" << std::endl;
	fout << "POINTS " << vertexCount << " float" << std::endl;

	int subSize, arraySize;
	std::vector<float> tempVec;
	for (int i = 0; i < array.size(); ++i)
	{
		tempVec = array[i];
		subSize = tempVec.size()/dimension;
		for (int j = 0; j < subSize; ++j)
		{
			for (int k = 0; k < dimension; ++k)
			{
				fout << tempVec[j*dimension+k] << " ";
			}
			fout << endl;
		}
	}

	fout << "LINES " << array.size() << " " << (vertexCount+array.size()) << std::endl;

	subSize = 0;
	for (int i = 0; i < array.size(); ++i)
	{
		arraySize = array[i].size()/dimension;
		fout << arraySize << " ";
		for (int j = 0; j < arraySize; ++j)
		{
			fout << subSize+j << " ";
		}
		subSize+=arraySize;
		fout << std::endl;
	}
	fout << "POINT_DATA" << " " << vertexCount << std::endl;
	fout << "SCALARS group int 1" << std::endl;
	fout << "LOOKUP_TABLE group_table" << std::endl;

	for (int i = 0; i < array.size(); ++i)
	{
		arraySize = array[i].size()/dimension;
		for (int j = 0; j < arraySize; ++j)
		{
			fout << i << std::endl;
		}
	}

	fout << "SCALARS silhouette float 1" << std::endl;
	fout << "LOOKUP_TABLE silhouette_table" << std::endl;

	for (int i = 0; i < array.size(); ++i)
	{
		arraySize = array[i].size()/dimension;
		for (int j = 0; j < arraySize; ++j)
		{
			fout << sCluster[i] << std::endl;
		}
	}
	fout.close();
}


/*
 * @brief print the streamlines with several scalar values in the vtk file
 *
 * @param[in] fileName The name of the vtk file
 * @param[in] array The streamline coordinates
 * @param[in] sCluster The float vector
 * @param[in] rotation The float vector
 * @param[in] dimension The dimension of the points
 */
void IOHandler::printFeature(const string& fileName,
							 const std::vector<std::vector<float> >& array,
							 const std::vector<float>& sCluster,
							 const std::vector<float>& rotation,
							 const int& dimension)
{
	if(array.empty() || sCluster.empty())
			return;
		stringstream ss;
		ss << "../dataset/" << fileName;
		ofstream fout(ss.str().c_str(), ios::out);
		if(!fout)
		{
			std::cout << "Error creating file!" << std::endl;
			exit(-1);
		}

		int vertexCount = 0;
		for (int i = 0; i < array.size(); ++i)
		{
			vertexCount += array[i].size();
		}
		vertexCount /= dimension;

		fout << "# vtk DataFile Version 3.0" << std::endl << "Bernard streamline" << std::endl
		     << "ASCII" << std::endl << "DATASET POLYDATA" << std::endl;
		fout << "POINTS " << vertexCount << " float" << std::endl;

		int subSize, arraySize;
		std::vector<float> tempVec;
		for (int i = 0; i < array.size(); ++i)
		{
			tempVec = array[i];
			subSize = tempVec.size()/dimension;
			for (int j = 0; j < subSize; ++j)
			{
				for (int k = 0; k < dimension; ++k)
				{
					fout << tempVec[j*dimension+k] << " ";
				}
				fout << endl;
			}
		}

		fout << "LINES " << array.size() << " " << (vertexCount+array.size()) << std::endl;

		subSize = 0;
		for (int i = 0; i < array.size(); ++i)
		{
			arraySize = array[i].size()/dimension;
			fout << arraySize << " ";
			for (int j = 0; j < arraySize; ++j)
			{
				fout << subSize+j << " ";
			}
			subSize+=arraySize;
			fout << std::endl;
		}
		fout << "POINT_DATA" << " " << vertexCount << std::endl;
		fout << "SCALARS group int 1" << std::endl;
		fout << "LOOKUP_TABLE group_table" << std::endl;

		for (int i = 0; i < array.size(); ++i)
		{
			arraySize = array[i].size()/dimension;
			for (int j = 0; j < arraySize; ++j)
			{
				fout << i << std::endl;
			}
		}

		fout << "SCALARS silhouette float 1" << std::endl;
		fout << "LOOKUP_TABLE silhouette_table" << std::endl;

		for (int i = 0; i < array.size(); ++i)
		{
			arraySize = array[i].size()/dimension;
			for (int j = 0; j < arraySize; ++j)
			{
				fout << sCluster[i] << std::endl;
			}
		}

		fout << "SCALARS rotation float 1" << std::endl;
		fout << "LOOKUP_TABLE rotation_table" << std::endl;

		for (int i = 0; i < array.size(); ++i)
		{
			arraySize = array[i].size()/dimension;
			for (int j = 0; j < arraySize; ++j)
			{
				fout << rotation[i] << std::endl;
			}
		}

		fout.close();
}


/*
 * @brief print the label and size of clusters for the streamlines
 *
 * @param[in] dataVec The streamline coordinates
 * @param[in] group The group labels of the streamlines
 * @param[in] storage The size of clusters
 * @param[in] groupName The string of the label name
 * @param[in] fullName The full name of the primary vtk file
 * @param[in] dimension The dimension of points
 */
void IOHandler::printClusters(const std::vector< std::vector<float> >& dataVec, 
							  const std::vector<int>& group, 
				 			  const std::vector<int>& storage, 
				 			  const string& groupName, 
				 			  const string& fullName, 
				 			  const int& dimension)
{
	if(group.empty()||storage.empty())
		return;
	std::ofstream fout(fullName.c_str(), ios::out | ios::app );
	if(!fout)
	{
		std::cout << "Error opening the file!" << std::endl;
		exit(1);
	}

	fout << "SCALARS " << groupName << " int 1" << std::endl;
	fout << "LOOKUP_TABLE " << groupName+string("_table") << std::endl;

	int arraySize;
	for (int i = 0; i < dataVec.size(); ++i)
	{
		arraySize = dataVec[i].size()/dimension;
		for (int j = 0; j < arraySize; ++j)
		{
			fout << group[i] << std::endl;
		}
	}

	fout << "SCALARS " <<  groupName + "_num" << " int 1" << std::endl;
	fout << "LOOKUP_TABLE " <<  groupName+string("_num_table") << std::endl;

	for (int i = 0; i < dataVec.size(); ++i)
	{
		arraySize = dataVec[i].size()/dimension;
		for (int j = 0; j < arraySize; ++j)
		{
			fout << storage[group[i]] << std::endl;
		}
	}
	fout.close(); 
}


/*
 * @brief print the clustered streamlines with noise information generated from DBSCAN
 *
 * @param[in] dataVec The streamline coordinates
 * @param[in] group The labels of streamlines
 * @param[in] storage The size of clusters
 * @param[in] groupName The string of the labels
 * @param[in] fullName The name of the primary vtk
 * @param[in] dimension The dimension of the points
 */
void IOHandler::printClustersNoise(const std::vector< std::vector<float> >& dataVec, 
							  const std::vector<int>& group, 
				 			  const std::vector<int>& storage, 
				 			  const string& groupName, 
				 			  const string& fullName, 
				 			  const int& dimension)
{
	/* in case you've noise, so group_id would be -1 */
	if(group.empty()||storage.empty())
		return;
	std::ofstream fout(fullName.c_str(), ios::out | ios::app );
	if(!fout)
	{
		std::cout << "Error opening the file!" << std::endl;
		exit(1);
	}

	fout << "SCALARS " << groupName << " int 1" << std::endl;
	fout << "LOOKUP_TABLE " << groupName+string("_table") << std::endl;

	int arraySize;
	for (int i = 0; i < dataVec.size(); ++i)
	{
		arraySize = dataVec[i].size()/dimension;
		for (int j = 0; j < arraySize; ++j)
		{
			fout << group[i] << std::endl;
		}
	}

	fout << "SCALARS " <<  groupName + "_num" << " int 1" << std::endl;
	fout << "LOOKUP_TABLE " <<  groupName+string("_num_table") << std::endl;

	for (int i = 0; i < dataVec.size(); ++i)
	{
		arraySize = dataVec[i].size()/dimension;
		for (int j = 0; j < arraySize; ++j)
		{
			fout << storage[group[i]+1] << std::endl;
		}
	}
	fout.close(); 
}


/*
 * @brief generate readme txt with number of clusters, input threshold and activity list, mostly for BIRCH
 *
 * @param[in] activityList The vector to store the activities
 * @param[in] timeList The vector to store the float of time
 * @param[in] normOption The norm option
 * @param[in] numClusters The number of clusters generated
 * @param[in] sValue The silhouette value
 * @param[in] threshold The distance threshold as input for the BIRCH clustering
 */
void IOHandler::generateReadme(const std::vector<string>& activityList,
							   const std::vector<double>& timeList,
							   const int& normOption,
							   const int& numClusters,
							   const float& sValue,
							   const float& threshold)
{
	if(activityList.empty()||timeList.empty())
		return;
	std::ofstream readme("../dataset/README",ios::out | ios::app);
	if(!readme)
	{
		std::cout << "Error creating readme!" << std::endl;
		exit(1);
	}
	readme << "----------------------------------------------" << std::endl;
	readme << "Norm: " << normOption << std::endl;
	readme << "Clusters: " << numClusters << std::endl;
	readme << "Silhouette: " << sValue << std::endl;
	readme << "Input threshold: " << threshold << std::endl;
	for (int i = 0; i < activityList.size(); ++i)
	{
		readme << activityList[i] << timeList[i] << " s." << std::endl;
	}
	readme << std::endl;
    readme.close();
}


/*
 * @brief print the readme information with activity list and time list that records some status during the computation
 *
 * @param[in] activityList The vector that has the relative activities
 * @param[in] timeList The vector of time for respective status
 */
void IOHandler::generateReadme(const std::vector<string>& activityList,
							   const std::vector<string>& timeList)
{
	if(activityList.empty()||timeList.empty())
		return;
	std::ofstream readme("../dataset/README",ios::out | ios::app);
	if(!readme)
	{
		std::cout << "Error creating readme!" << std::endl;
		exit(1);
	}
	readme << "----------------------------------------------" << std::endl;
	for (int i = 0; i < activityList.size(); ++i)
	{
		readme << activityList[i] << timeList[i] << std::endl;
	}
    readme.close();
}


/*
 * @brief write string arrays in the readme
 *
 * @param[in] comments The string to be put in the readme
 */
void IOHandler::writeReadme(const string& comments)
{
	if(comments.empty())
		return;
	std::ofstream readme("../dataset/README",ios::out | ios::app);
	if(!readme)
	{
		std::cout << "Error creating readme!" << std::endl;
		exit(1);
	}
	readme << comments << std::endl;
    readme.close();
}


/*
 * @brief write the finalized group size in the storage
 *
 * @param[in] storage The size of clusters from the clustering results
 */
void IOHandler::writeGroupSize(const std::vector<int>& storage)
{
	if(storage.empty())
		return;
	std::ofstream readme("../dataset/README",ios::out | ios::app);
	if(!readme)
	{
		std::cout << "Error creating readme!" << std::endl;
		exit(1);
	}
	readme << "Final cluster size: " << storage.size() << std::endl;
	for (int i = 0; i < storage.size(); ++i)
	{
		readme << storage[i] << " ";
	}
	readme << std::endl;
    readme.close();
}


/*
 * @brief Write the candidates of each cluster into the local file
 *
 * @param[in] storage The candidates belonging to each cluster
 */
void IOHandler::generateGroups(const std::vector<std::vector<int> >& storage)
{
	if(storage.empty())
		return;
	std::ofstream readme("../dataset/Storage",ios::out|ios::app);
	if(!readme)
	{
		std::cout << "Error creating Storage!" << std::endl;
		exit(1);
	}

	readme << std::endl;
	const int& groupSize = storage.size();
	std::vector<int> element;
	for(int i=0;i<groupSize;++i)
	{
		element = storage[i];
		if(element.empty())
			continue;
		for(int j=0;j<element.size();++j)
			readme << element[j] << " ";
		readme << std::endl;
	}
	std::cout << std::endl;
	readme.close();
}


/*
 * @brief Generate and store the candidates of all the clusters in the local file
 *
 * @param[in] storage The candidates belonging to each cluster
 * @param[in] fileName The name of the file for storing the candidiate information
 */
void IOHandler::generateGroups(const std::vector<std::vector<int> >& storage, const string& fileName)
{
	if(storage.empty())
		return;
	std::ofstream readme(("../dataset/"+fileName).c_str(),ios::out);
	if(!readme)
	{
		std::cout << "Error creating Storage!" << std::endl;
		exit(1);
	}

	readme << std::endl;
	const int& groupSize = storage.size();
	std::vector<int> element;
	for(int i=0;i<groupSize;++i)
	{
		element = storage[i];
		if(element.empty())
			continue;
		for(int j=0;j<element.size();++j)
			readme << element[j] << " ";
		readme << std::endl;
	}
	std::cout << std::endl;
	readme.close();
}


/*
 * @brief Print float and string information into the readme
 *
 * @param[in] value The float value
 * @param[in] dataSet The string type
 * @param[in] clustering The string of the clustering technique
 * @param[in] value_name The string for the value name
 */
void IOHandler::writeReadMe(const float& value, const string& dataSet, const string& clustering,
							const string& value_name)
{
	std::ofstream out_file("../dataset/README", ios::out|ios::app);
	if (!out_file)
	{
		std::cout << "Error for creating README!" << std::endl;
		exit(1);
	}
	out_file << "-----------------------------------------------------------------" << std::endl;
	out_file << "The " << value_name << " of " << clustering << " on dataset "
			 << dataSet << " is " << value << std::endl;
	out_file.close();
}


/*
 * @brief Read the number of clusters as input for different similarity measures from the local file
 *
 * @param[out] clusMap The hash map that records the number of cluster and its norm option
 * @param[in] fileName The name of the file to be read in
 */
void IOHandler::readClusteringNumber(std::unordered_map<int,int>& clusMap, const string& fileName)
{
	std::ifstream readme(("../dataset/"+fileName).c_str(), ios::in);
	if(!readme)
	{
		std::cout << "Error creating Storage!" << std::endl;
		exit(1);
	}
	string line;
	int scopePos;
	while(getline(readme, line))
	{
		scopePos = line.find(":");
		if(scopePos==std::string::npos)
		{
			std::cout << "Error for clustering number reader..." << std::endl;
			exit(1);
		}
		clusMap[std::atoi(line.substr(0,scopePos).c_str())] = std::atoi(line.substr(scopePos+1).c_str());
	}
	readme.close();
}

