#include <fstream>
#include <vector>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <sstream>

using namespace std;


void readFile(const char* fileName, 
			  std::vector<std::vector<float> >& dataVec, 
			  const int& lineNumber, 
			  const int& vertexNumber);

void writeFile(const char* fileName, 
			   const std::vector<std::vector<float> >& dataVec);


int main()
{
	std::vector<std::vector<float> > dataVec;
	const char* fileName = "hurricane";

	std::cout << "input the number among which one is chosen for streamlines? " << std::endl;
	int lineNumber;
	std::cin >> lineNumber;

	std::cout << "input the number among which one is chosen for vertex? " << std::endl;
	int vertexNumber;
	std::cin >> vertexNumber;

	readFile(fileName, dataVec, lineNumber, vertexNumber);
	writeFile("Hurricane", dataVec);
	return 0;
}


void readFile(const char* fileName, 
			  std::vector<std::vector<float> >& dataVec, 
			  const int& lineNumber, 
			  const int& vertexNumber)
{
	std::ifstream fin(fileName, ios::in);
	if(!fin)
	{
		std::cout << "Error creating files!" << std::endl;
		exit(1);
	}
	stringstream ss;
	std::vector<float> tempVec;

	string line, part;

	std::vector<float> vec(3);
	float temp;

	int lineTag = 0;
	while(getline(fin, line) /* && currentNumber < MAXNUMBER*/)
	{
		//currentDimensions = 0;
		if(lineTag==1)
		{
			lineTag = (lineTag+1)%lineNumber;
			continue;
		}

		int tag = 0, count = 0;
		bool isNext = false;
		ss.str(line);
		while(ss>>part /*&& currentDimensions<3*MAXDIMENSION*/)
		{
			/* operations below would remove duplicate vertices because that would damage our computation */
			if(tag>=3)
			{
				isNext = !isNext;
				tag = (tag+1)%(vertexNumber*3);
				continue;
			}
			temp = atof(part.c_str());
			if(isNext)
			{
				if(count<3)
				{
					vec[count] = temp;
					tag = (tag+1)%(vertexNumber*3);
					++count;
				}
				if(count==3)
				{
					int size = tempVec.size();
					if(!(abs(vec[0]-tempVec[size-3])<1.0e-5&&abs(vec[1]-tempVec[size-2])<1.0e-5&&abs(vec[2]-tempVec.back())<1.0e-5))
					{
						tempVec.push_back(vec[0]);
						tempVec.push_back(vec[1]);
						tempVec.push_back(vec[2]);
					}
					count = 0;
				}
				continue;
			}
			tempVec.push_back(temp);
			tag = (tag+1)%(vertexNumber*3);
			//currentDimensions++;
		}
		/* accept only streamlines with at least three vertices */
		if(tempVec.size()/3>2)
		{
			dataVec.push_back(tempVec);
		}
		tempVec.clear();
		ss.clear();
		ss.str("");
		//currentNumber++;

		lineTag = (lineTag+1)%lineNumber;
	}
	fin.close();

	std::cout << "Finished reading file!" << std::endl;
}

void writeFile(const char* fileName, 
			   const std::vector<std::vector<float> >& dataVec)
{
	std::ofstream fout(fileName, ios::out);
	if(!fout)
	{
		std::cout << "Cannot create a file!" << std::endl;
		exit(1);
	}

	const int& vecSize = dataVec.size();
	std::vector<float> tempVec;
	int tempVecSize;

	for (int i = 0; i < vecSize; ++i)
	{
		tempVec = dataVec[i];
		tempVecSize = tempVec.size();
		for (int j = 0; j < tempVecSize; ++j)
		{
			fout << tempVec[j] << " ";
		}
		fout << std::endl;
	}
	fout.close();

	std::cout << "Finished writing file!" << std::endl;
}