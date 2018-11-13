/*
 * ReadClustering.cpp
 *
 *  Created on: Mar 13, 2018
 *      Author: lieyu
 */

#include "ReadClustering.h"

ReadClustering::ReadClustering() {
	// TODO Auto-generated constructor stub

}

ReadClustering::~ReadClustering() {
	// TODO Auto-generated destructor stub
}


/* the public function called by main.cpp */
void ReadClustering::getEvaluation(const char* fileName)
{
	int isPBFInput;
	std::cout << "Is it a PBF dataset? 1.Yes, 0.No." << std::endl;
	std::cin >> isPBFInput;
	assert(isPBFInput==1||isPBFInput==0);
	isPBF = (isPBFInput==1);

	/* read data into ds */
	std::cout << fileName << std::endl;
	readData(fileName);

	/* compute the evaluation */
	computeEvaluation();

	/* output the result to text file */
	writeAnalysis();
}


/* write those analysis framework */
void ReadClustering::writeAnalysis()
{
	/* write information */
	IOHandler::generateReadme(activityList,timeList);

}


/* read data from file and store it into Dataset ds */
void ReadClustering::readData(const char* fileName)
{
	ifstream fin(fileName, ios::in);
	if(!fin)
	{
		std::cout << "Error opening the file!" << std::endl;
		exit(1);
	}

	Eigen::MatrixXf vertexCoordinate;

	/* omit first four lines */
	string line;
	for(int i=0;i<5;++i)
	{
		getline(fin, line);
	}
	/* split the string into three parts */
	stringstream ss(line);
	ss >> line;
	ss >> line;

	ss.str(std::string());

	/* get how many vertex inside */
	const int& vertexCount = atoi(line.c_str());

	vertexCoordinate = Eigen::MatrixXf(vertexCount, 3);

	/* read in vertex coordinates */
	for(int i=0;i<vertexCount;++i)
	{
		/* read one line */
		getline(fin, line);
		/* split and analyze the string */
		ss.str(line);
		for(int j=0;j<3;++j)
		{
			ss >> line;
			vertexCoordinate(i,j) = atof(line.c_str());
		}
		ss.str(std::string());
	}

	/* get how many streamlines you'll have */
	getline(fin, line);

	ss.str(line);
	ss >> line;
	ss >> line;

	ss.str(std::string());

	ds.numOfElements = atoi(line.c_str());

	ds.dataVec = std::vector<std::vector<float> >(ds.numOfElements);

	/* read vertex coordinates into dataVec */

	int maxElements = INT_MIN;

	int vertexNum, index;
	for(int i=0;i<ds.numOfElements;++i)
	{
		getline(fin,line);
		ss.str(line);

		/* explicate vertex count in each line */
		ss>>line;
		vertexNum = atoi(line.c_str());

		/* assign memory */
		std::vector<float>& tempVec = ds.dataVec[i];

		tempVec = std::vector<float>(vertexNum*3);

		maxElements = std::max(maxElements, vertexNum*3);

		for(int j=0;j<vertexNum;++j)
		{
			ss>>line;
			index = atoi(line.c_str());
			for(int k=0;k<3;++k)
				tempVec[3*j+k] = vertexCoordinate(index,k);
		}

		ss.str(std::string());
	}
	for(int i=0;i<3;++i)
	{
		getline(fin,line);
	}
	for(int i=0;i<vertexCount;++i)
	{
		getline(fin,line);
	}
	std::size_t found_int, found, found_scalars;

	int normOption, totalLine, groupIndex;
	while(getline(fin,line))
	{
		found_int = line.find("int");
		found_scalars = line.find("SCALARS");
		/* has int, should be group information */
		if(found_int!=std::string::npos && found_scalars!=std::string::npos)
		{
			ss.str(line);
			ss>>line;
			ss>>line;

			string norm_choice;
			if(strcmp(line.substr(0,3).c_str(), "PCA")==0)
				norm_choice="PCA";
			else
			{
				found = line.find("_");
				found_int = line.find("norm");
				if(found==std::string::npos)
				{
					norm_choice = line;
				}
				else if(found_int>found)
				{
					norm_choice = line.substr(found_int);
				}
				else if(found_int<found)
				{
					norm_choice = line.substr(found_int, found);
				}
			}

			getline(fin,line);

			std::vector<int> tempGroup(ds.numOfElements);

			for(int i=0;i<ds.numOfElements;++i)
			{
				for(int j=0;j<ds.dataVec[i].size()/3;++j)
					getline(fin,line);
				tempGroup[i] = atoi(line.c_str());
			}
			ds.groupAggregate.insert(std::make_pair(norm_choice, tempGroup));

			totalLine = vertexCount+2;
			for(int i=0; i<totalLine; ++i)
			{
				getline(fin,line);
			}
		}
	}

	fin.close();

	IOHandler::sampleArray(ds.array,ds.dataVec,3,maxElements);

	/* compute cluster number */

	std::vector<int> groupArray;
	int max_num;
	std::unordered_map<string, std::vector<int> >::const_iterator iter;
	for(iter=ds.groupAggregate.begin(); iter!=ds.groupAggregate.end(); ++iter)
	{
		groupArray = iter->second;
		max_num = -1;
		if(groupArray.empty())
			continue;

		index = groupArray.size();
		for(int j=0;j<index;++j)
		{
			max_num = std::max(max_num, groupArray[j]);
		}
		max_num+=1;
		ds.maxGroup.insert(make_pair(iter->first, max_num));
	}

}

/* compute four analysis evaluation measures */
void ReadClustering::computeEvaluation()
{
	std::unordered_map<string, std::vector<int> >::const_iterator iter;
	for(iter=ds.groupAggregate.begin(); iter!=ds.groupAggregate.end(); ++iter)
	{
		std::cout << "Processing for " << iter->first << "..." << std::endl;
		computeEvaluation(iter);
	}
}



/* compute evaluation based on norm option */
void ReadClustering::computeEvaluation(std::unordered_map<string, std::vector<int> >::const_iterator& iter)
{
	if(iter->second.empty())
		return;

	ds.neighborVec.clear();
	ds.neighborVec = std::vector<std::vector<int> >(ds.maxGroup[iter->first]);

	const std::vector<int>& groupOfNorm = iter->second;
	const int& groupSize = groupOfNorm.size();

	Silhouette sil;
	ValidityMeasurement vm;

	int totalNum = 0;
	for(int i=0;i<groupSize;++i)
	{
		if(groupOfNorm[i]<0)
			continue;
		ds.neighborVec[groupOfNorm[i]].push_back(i);
		++totalNum;
	}

	/* the 'PCA' option */
	if(strcmp("PCA", iter->first.c_str())==0)
	{
		sil.computeValue(ds.array,groupOfNorm,ds.maxGroup[iter->first],isPBF);
		vm.computeValue(ds.array, groupOfNorm);
	}
	else
	{
		/* count from "norm" for norm option */
		const int& normOption = std::atoi(iter->first.substr(4).c_str());
		std::cout << "This is norm " << normOption << std::endl;
		//if(normOption!=4 && normOption!=15)
		//	return;

		MetricPreparation object(ds.array.rows(), ds.array.cols());
		object.preprocessing(ds.array, ds.array.rows(), ds.array.cols(), normOption);
		/* if the dataset is not PBF, then should record distance matrix for Gamma matrix compution */
		if(!isPBF)
		{
			deleteDistanceMatrix(ds.array.rows());
			getDistanceMatrix(ds.array, normOption, object);
			std::cout << "Distance between 0 and 1 is " << distanceMatrix[0][1] << std::endl;
		}
		sil.computeValue(normOption,ds.array,ds.array.rows(),ds.array.cols(),groupOfNorm,object,
					         ds.maxGroup[iter->first],isPBF, ds.neighborVec);
		vm.computeValue(normOption, ds.array, groupOfNorm, object, isPBF);
	}

	/* compute the entropy */
	float entropy = 0, prob;
	for (int i = 0; i < ds.neighborVec.size(); ++i)
	{
		if(ds.neighborVec[i].size()>0)
		{
			prob = float(ds.neighborVec[i].size())/float(totalNum);
			entropy+=prob*log2f(prob);
		}
	}

	entropy = -entropy/log2f(ds.maxGroup[iter->first]);
	std::cout << "Entropy is " << entropy << std::endl;

	activityList.push_back("Silhouette for "+iter->first+" is: ");
	timeList.push_back(to_string(sil.sAverage));

	activityList.push_back("Gamma statistic for "+iter->first+" is: ");
	timeList.push_back(to_string(sil.gammaStatistic));

	activityList.push_back("Entropy for "+iter->first+" is: ");
	timeList.push_back(to_string(entropy));

	activityList.push_back("DB Index for "+iter->first+" is: ");
	timeList.push_back(to_string(sil.dbIndex));

	activityList.push_back("Validity measurement on "+iter->first+" is: ");
	stringstream fc_ss;
	fc_ss << vm.f_c;
	timeList.push_back(fc_ss.str());

	/* record labeling information */
	// IOHandler::generateGroups(ds.neighborVec, iter->first+"_storage");
}

