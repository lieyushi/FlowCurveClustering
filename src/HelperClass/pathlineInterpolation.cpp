/*
 * This file is for time-based sampling for blood_flow.vtk file. The reason why time-based sampling is needed is
 * because centroid computation between two pathlines are needed in clustering evaluation (clustering algorithms,
 * evaluation metric computation and cluster representative selection). Due to the special properties of pathlines,
 * the point coordinates at the same time frame of different pathlines should be corresponding to same time. Hence
 * time-based sampling strategy is required and needed!
 *
 * This resampling is only suitable for blood_flow.vtk and for other files specific format of functions must be re-
 * designed!!!!!!
 */

#include <fstream>
#include <vector>
#include <algorithm>
#include <iostream>
#include <cstring>
#include <sstream>
#include <stdio.h>
#include <string.h>
#include <map>
#include <climits>
#include <cassert>
#include <tuple>
#include <float.h>
#include <unordered_map>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/SVD>

using namespace std;
using namespace Eigen;

#define MULTIPLIER 8.0

// the pathline point with coordinate and time
struct PathlinePoint
{
	float coordinates[3];
	float time;

	PathlinePoint(): time(-1)
	{
		for(int i=0; i<3; ++i)
			coordinates[i] = -1.0;
	}
};

// read pathlinePoint from .vtk file
void readPathlineRaw(const char* fileName, std::vector<std::vector<PathlinePoint> >& pathlines,
		std::tuple<float,float, float>& timeRange);

// perform interpolation with same time slice for blood flow pathlines
void performInterpolation(const std::vector<std::vector<PathlinePoint> >& pathlines,
		std::vector<Eigen::VectorXf>& interpolatedLine, const std::tuple<float,float,float>& timeRange);

// write into txt file for clustering evaluation of blood flow
void generateLineFile(const std::vector<Eigen::VectorXf>& interpolatedLine, const char* fileName);

int main(int argc, char* argv[])
{
	if(argc!=2)
	{
		std::cout << "Error for argument count. Should be ./executable fileName" << std::endl;
		exit(1);
	}

	// get the proper file name
	auto pos = string(argv[1]).find(".vtk");
	if(pos==std::string::npos)
	{
		std::cout << "Input file is not a .vtk file!" << std::endl;
		exit(1);
	}

	std::vector<std::vector<PathlinePoint> > pathlineRaw;
	std::tuple<float,float,float> timeRange;
	// read raw pathlines with coordinates and time from .vtk file
	readPathlineRaw(argv[1], pathlineRaw, timeRange);

	// interpolate the pathlines so that points of different frames for pathlines reside the exactly same time slide
	std::vector<Eigen::VectorXf> interpolatedLine;
	performInterpolation(pathlineRaw, interpolatedLine, timeRange);

	// write into txt file for clustering evaluation of blood flow
	generateLineFile(interpolatedLine, string(argv[1]).substr(0,pos).c_str());

	return 0;
}


// read pathlinePoint from .vtk file
void readPathlineRaw(const char* fileName, std::vector<std::vector<PathlinePoint> >& pathlines,
		std::tuple<float,float, float>& timeRange)
{	
	std::get<0>(timeRange) = FLT_MAX;
	std::get<1>(timeRange) = -1.0;

	std::vector<PathlinePoint> pointCoordinateVec;

	std::ifstream fin(fileName, ios::in);

	if(fin.fail())
	{
		std::cout << "Error for opening file contents!" << std::endl;
		exit(1);
	}

	string line;

	// ignore the header file of vtk
	for(int i=0; i<5; ++i)
		getline(fin, line);

	stringstream ss(line);
	ss >> line, ss >> line; // get number of points

	ss.clear();
	ss.str("");

	// extract the number of points
	const int& numOfPoints = std::atoi(line.c_str());
	std::cout << "There are " << numOfPoints << " points read from the file!" << std::endl;

	// assign the memory for numOfPoints
	pointCoordinateVec.resize(numOfPoints);

	// start traversal for point coordinates
	int index = 0, separation = 0;
	while(getline(fin, line) && line.size()!=0)
	{
		ss.str(line);
		separation = 0;

		// read the one line of points into the cache
		while(ss>>line)
		{
			pointCoordinateVec[index].coordinates[separation] = std::atof(line.c_str());
			if(separation == 2)
			{
				separation = 0;
				++index;
			}
			else
				++separation;
		}
		ss.clear();
		ss.str("");
	}

	assert(index == numOfPoints);

	// read the number of lines from the txt file
	for(int i=0; i<8; ++i)
		getline(fin, line);

	ss.str(line);
	ss >> line, ss >> line;
	const int& numOfLines = std::atoi(line.c_str());
	std::cout << "There are " << numOfLines << " pathlines read from the file!" << std::endl;
	ss.clear();
	ss.str("");

	// store the line to point arrays, e.g., one line has which indices of points
	std::vector<vector<int> > pointsToLine(numOfLines);

	int numOfPointsForLine, pointIndex;
	for(int i=0; i<numOfLines; ++i)
	{
		getline(fin, line);
		ss.str(line);
		ss >> line;

		// find how many points this line contain
		numOfPointsForLine = std::atoi(line.c_str());

		std::vector<int>& currentLine = pointsToLine[i];
		currentLine.resize(numOfPointsForLine);

		for(int j=0; j<numOfPointsForLine; ++j)
		{
			ss >> line;
			pointIndex = std::atoi(line.c_str());
			currentLine[j] = pointIndex;
		}

		ss.clear();
		ss.str("");
	}

	while(getline(fin, line))
	{
		ss.str(line);
		ss >> line;
		ss.clear();
		ss.str("");
		if(strcmp(line.c_str(), "time")==0)
			break;
	}

	index = 0; 
	while(getline(fin, line) && line.size())
	{
		ss.str(line);
		while(ss >> line)
		{
			pointCoordinateVec[index].time = std::atof(line.c_str());
			++index;
		}
		ss.clear();
		ss.str("");
	}
	assert(index == numOfPoints);
	fin.close();

	std::cout << "File content traversal completed!" << std::endl;

	pathlines.resize(numOfLines);
	std::vector<int> lineIndex;

	float averageSlice = 0.0;
	for(int i=0; i<numOfLines; ++i)
	{
		std::vector<PathlinePoint>& currentLine = pathlines[i];
		lineIndex = pointsToLine[i];
		currentLine.resize(lineIndex.size());

		for(int j=0; j<lineIndex.size(); ++j)
		{
			currentLine[j] = pointCoordinateVec[lineIndex[j]];
		}
		std::get<0>(timeRange) = std::min(std::get<0>(timeRange), currentLine[0].time);
		std::get<1>(timeRange) = std::max(std::get<1>(timeRange), currentLine.back().time);

		averageSlice += (currentLine.back().time-currentLine[0].time)/float(lineIndex.size()-1);
		std::cout << "Pathline " << i << " has starting time " << currentLine[0].time <<
		" and ending time " << currentLine.back().time << std::endl;
	}
	pointCoordinateVec.clear();
	std::get<2>(timeRange) = averageSlice/float(numOfLines);

	std::cout << "Starting time is " << std::get<0>(timeRange) << ", ending time is " << std::get<1>(timeRange) << 
	", and average time slice is " << std::get<2>(timeRange) << std::endl;
}

// perform interpolation with same time slice for blood flow pathlines
void performInterpolation(const std::vector<std::vector<PathlinePoint> >& pathlines,
		std::vector<Eigen::VectorXf>& interpolatedLine, const std::tuple<float,float,float>& timeRange)
{
	// use the time information of timeRange to perform time-based sampling for the pathlines
	const float& starting = std::get<0>(timeRange);
	const float& ending = std::get<1>(timeRange);
	const float& aveSlice = MULTIPLIER*std::get<2>(timeRange);
	interpolatedLine.resize(pathlines.size());

#pragma omp parallel for schedule(static) num_threads(8)
	for(int i=0; i<pathlines.size(); ++i)
	{
		// for each pathlines, will interpolate it from 0, 1, 2, ...., currentEndingTime
		const std::vector<PathlinePoint>& line = pathlines[i];
		std::map<float, Eigen::Vector3f> timeCoordinates;
		for(int j=0; j<line.size(); ++j)
		{
			timeCoordinates[line[j].time] = Eigen::Vector3f(line[j].coordinates[0], line[j].coordinates[1],
					line[j].coordinates[2]);
		}
		int numOfPoints = int((line.back().time-starting)/(aveSlice))+1;

		Eigen::VectorXf& lineCoordinate = interpolatedLine[i];
		lineCoordinate = Eigen::VectorXf(3*numOfPoints);

		float current;
		for(int j=0; j<numOfPoints; ++j)
		{
			current = starting + aveSlice*j;

			// before the occuring time, should be directly repeating the first point
			if(current<=line[0].time)
			{
				for(int k=0; k<3; ++k)
				{
					lineCoordinate(3*j+k) = line[0].coordinates[k];
				}
			}
			// this time has been recorded in the map, directly load the data
			else if(timeCoordinates.find(current)!=timeCoordinates.end())
			{
				for(int k=0; k<3; ++k)
				{
					lineCoordinate(3*j+k) = timeCoordinates[current](k);
				}
			}
			// else, find the left and right time step and perform linear interpolation
			else
			{
				// find the clipping left and right time slices for interpolation
				auto right = timeCoordinates.upper_bound(current);
				auto left = std::prev(right);

				float ratio = (current-left->first)/(right->first-left->first);
				Eigen::Vector3f currentPoint = (1.0-ratio)*left->second + ratio*right->second;
				for(int k=0; k<3; ++k)
				{
					lineCoordinate(3*j+k) = currentPoint(k);
				}
			}
		}
	}
}

/*
 * write into txt file for clustering evaluation of blood flow, each line is point coordinates
 * separated by blankspace, x1 y1 z1 x2 y2 z2
 */
void generateLineFile(const std::vector<Eigen::VectorXf>& interpolatedLine, const char* fileName)
{
	std::ofstream fout(fileName, ios::out);
	if(fout.fail())
	{
		std::cout << "Error for creating ascii file in the folder!" << std::endl;
		exit(1);
	}
	// get the how many pathlines
	const int& lineSize = interpolatedLine.size();
	Eigen::VectorXf pathlineData;
	for(int i=0; i<lineSize; ++i)
	{
		pathlineData = interpolatedLine[i];

		for(int j=0; j<pathlineData.size(); ++j)
		{
			fout << pathlineData(j) << " ";
		}
		fout << std::endl;
	}

	fout.close();
}
