#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <sstream>
#include <cmath>

using namespace std;

void readFile(std::vector<int>& groupSize, const char* fileName);

void computeEntropy(const std::vector<int>& groupSize);

int main(int argc, char* argv[])
{
	if(argc!=2)
	{
		std::cout << "Error for argument input!" << std::endl;
		exit(1);
	}

	std::vector<int> storage;
	readFile(storage, argv[1]);

	computeEntropy(storage);

	return 0;
}


void readFile(std::vector<int>& groupSize, const char* fileName)
{
	ifstream fin(fileName, ios::in);
	if(!fin)
	{
		std::cout << "Error for reading a file!" << std::endl;
		exit(1);
	}

	string line;

	int num;
	while(getline(fin,line))
	{
		num = 0;
		stringstream ss(line);
		while(ss>>line)
			++num;
		if(num > 0)
			groupSize.push_back(num);
	}

	fin.close();
}

void computeEntropy(const std::vector<int>& groupSize)
{
	int total = 0;
	for(int i=0;i<groupSize.size();++i)
		total+=groupSize[i];

	float prob, result = 0.0;
	for(int i=0;i<groupSize.size();++i)
	{
		prob = float(groupSize[i])/float(total);
		result+=prob*log2f(prob);
	}

	result = -result/log2f(groupSize.size());
	std::cout << "Entropy is " << result << std::endl;
}