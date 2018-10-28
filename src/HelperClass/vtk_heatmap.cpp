#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <cstring>
#include <sstream>
#include <float.h>
#include <algorithm>
#include <unordered_map>
#include <iomanip> 

using namespace std;

const float& range_start = 0.1;

/* db index max than this threshold will be taken as maximal */
const float& max_db_index = 5.0;

/*-------------------------------------------------------------
	for silhouette and gamma, just use linear mapping into [range_start, 1.0]
	for db index, for minimal <= db <= max_db_index, use linear mapping, 
				  for max_db_index <= db, directly assigns it as 1.0
    for validity measurement, use log10 for linear mapping 
  -------------------------------------------------------------
*/

float data_range[4][3]=
{
	FLT_MAX, -FLT_MAX, 0,
	FLT_MAX, -FLT_MAX, 0,
	FLT_MAX, -FLT_MAX, 0,
	FLT_MAX, -FLT_MAX, 0
};

struct AverageColumn
{
	float average;
	string name;
	std::vector<float> valueVec;
	std::vector<float> std_vec;
	AverageColumn(const float& average, const string& name): average(average), name(name)
	{}

	AverageColumn(): average(0.0)
	{}
};

struct AverageClustering
{
	float average;
	int originalIndex;
	AverageClustering(const float& average, const int& index): average(average), originalIndex(index)
	{}

	AverageClustering(): average(0.0)
	{}
};


struct BestValue
{
	float value;
	int i, j;

	BestValue(const float& value): value(value), i(-1), j(-1)
	{}

	BestValue(): value(0.0), i(-1), j(-1)
	{}

};

void readData(std::vector<std::vector<float> >& dataVec, const char* fileName);

void createHeatMap(const std::vector<std::vector<float> >& dataVec);

int make_index(const int& i, const int& j, const int& col);

void create_assemble(const std::vector<std::vector<float> >& dataVec);

void create_separate(const std::vector<std::vector<float> >& dataVec);

void create_ranking(const std::vector<std::vector<float> >& dataVec);

void create_latex_table(const std::vector<std::vector<float> >& dataVec);

void get_average_value(std::vector<std::vector<float> >& averageValue, std::vector<std::vector<float> >& standardDeviation,
	              string file_list[], const int& file_size);

void create_std_ranking(const std::vector<std::vector<float> >& averageValue, const std::vector<std::vector<float> >& standardDeviation,
	                    const int& file_size);


int main(int argc, char* argv[])
{
	if(argc!=2)
	{
		std::cout << "Error for argument count!" << std::endl;
		exit(1);
	}

	std::vector<std::vector<float> > dataVec;

	readData(dataVec, argv[1]);

	createHeatMap(dataVec);

	create_assemble(dataVec);

	//create_separate(dataVec);

	create_ranking(dataVec);

	create_latex_table(dataVec);
	

	/*if(argc!=1)
	{
		std::cout << "Get average and std so no need for argument!" << std::endl;
		exit(1);
	}
	std::vector<std::vector<float> > averageValue, standardDeviation;

	string filenames[] = {"bernard_evaluation", "crayfish_evaluation", "cylinder_evaluation", "hurricane_evaluation", "solar_plume_evaluation", "tornado_evaluation"};

	//string filenames[] = {"cylinder_pathlines_evaluation", "tub_pathlines_evaluation"};

	get_average_value(averageValue, standardDeviation, filenames, sizeof(filenames)/sizeof(string));

	create_std_ranking(averageValue, standardDeviation, sizeof(filenames)/sizeof(string));

	create_assemble(averageValue);

	create_latex_table(averageValue);*/

	return 0;
}


void readData(std::vector<std::vector<float> >& dataVec, const char* fileName)
{
	std::ifstream fin(fileName, ios::in);
	if(fin.fail())
	{
		std::cout << "Error for reading data from existing file!" << std::endl;
		exit(1);
	}

	stringstream ss;
	string line;
	std::vector<float> row;
	while(getline(fin, line))
	{
		ss.str(std::string());
		ss.clear();
		ss << line;
		while(ss>>line)
		{
			if(strcmp(line.c_str(), "-")==0)
			{
				row.push_back(-10000.0);
			}
			else
			{
				row.push_back(std::atof(line.c_str()));
			}
		}
		dataVec.push_back(row);
		row.clear();
	}

	fin.close();
}


int make_index(const int& i, const int& j, const int& col)
{
	return j*col+i;
}


void createHeatMap(const std::vector<std::vector<float> >& dataVec)
{
	/* get the limit range of four scalar values */
	const int& rows = dataVec.size();
	const int& cols = dataVec[0].size();

	float value;
	int num;
	int nonZero = 0;
	for (int i = 0; i < rows; ++i)
	{
		for (int j = 0; j < cols; ++j)
		{
			if(dataVec[i][j]<=-9999.0)
				continue;
			num = (i%2)*2+j%2;
			++nonZero;
			data_range[num][0] = std::min(data_range[num][0], dataVec[i][j]);
			data_range[num][1] = std::max(data_range[num][1], dataVec[i][j]);

		}
	}
	for (int i = 0; i < 4; ++i)	
	{
		std::cout << data_range[i][0] << " " << data_range[i][1] << std::endl;
		if(i==3)
		{
			data_range[i][0] = log10(data_range[i][0]);
			data_range[i][1] = log10(data_range[i][1]);
		}
		else if(i==2)
		{
			if(data_range[i][1]>=max_db_index)
				data_range[i][1]=max_db_index;
		}
		data_range[i][2] = data_range[i][1]-data_range[i][0];
	}

	/* generate heatmap values */
	std::ofstream ofs("heatmap.vtk", ios::out);
	if(ofs.fail())
	{
		std::cout << "Error for creating vtk file!" << std::endl;
		exit(1);
	}
	ofs << "# vtk DataFile Version 3.0\n"
        << "matrix_vis" << "\n"
        << "ASCII\n\n"
        << "DATASET POLYDATA\n";
    ofs << "POINTS " << (rows+1)*(cols+1) << " float\n";
    const float& x_step = 0.1;
    const float& y_step = 0.1;
    for (int j = 0; j < rows+1; ++j)
    {
    	for (int i = 0; i < cols+1; ++i)
    	{
    		ofs << i*x_step << " " << j*y_step << " " << 0 << std::endl;
    	}
    }
    ofs << "POLYGONS " << nonZero << " " << 5*nonZero << "\n";
    int x, y;
    for (int j = rows-1; j >=0; --j)
    {
    	for (int i = 0; i<cols; ++i)
    	{
    		if(dataVec[j][i]<=-9999.0)
    			continue;
    		ofs << 4 << " " << make_index(i,rows-1-j,cols+1) << " " << make_index(i+1,rows-1-j,cols+1) << " "
    		    << make_index(i+1,rows-j,cols+1) << " " << make_index(i,rows-j,cols+1) << std::endl;
    	}
    }
    ofs << "CELL_DATA " << nonZero << "\n" << "SCALARS " << "label" << " float 1\n" << "LOOKUP_TABLE default\n";
    for (int j = rows-1; j >=0; --j)
    {
    	for (int i = 0; i < cols; ++i)
    	{
    		if(dataVec[j][i]<=-9999.0)
    			continue;
    		num = (j%2)*2+i%2;

    		if(num<=1)
    		{
    			value = (dataVec[j][i]-data_range[num][0])/data_range[num][2];
    		}
    		else if(num==3)
    		{
    			value = (log10(dataVec[j][i])-data_range[num][0])/data_range[num][2];
    		}
    		else if(num==2)
    		{
    			if(dataVec[j][i]>=max_db_index)
    				value = 1.0;
    			else
    				value = (dataVec[j][i]-data_range[num][0])/data_range[num][2];
    		}

    		if(num==0 || num==1)
    			ofs << value << std::endl;
    		else
    			ofs << 1.0-value << std::endl;

    	}
    }
    ofs.close();

    /* generate boundary grids for 2X2 */
    std::ofstream grid("grid.vtk", ios::out);
    if(grid.fail())
    {
    	std::cout << "Error for creating file!" << std::endl;
    	exit(1);
    }
    grid << "# vtk DataFile Version 3.0\n"
        << "matrix_vis" << "\n"
        << "ASCII\n\n"
        << "DATASET POLYDATA\n";
    grid << "POINTS " << (rows/2+1)*(cols/2+1) << " float\n";
    for (int j = 0; j < rows/2+1; ++j)
    {
    	for (int i = 0; i < cols/2+1; ++i)
    	{
    		grid << i*2.0*x_step << " " << j*2.0*y_step << " " << 0 << std::endl;
    	}
    }
    const int& line_number = rows/2*(cols/2+1)+(rows/2+1)*cols/2;
    grid << "LINES " << line_number << " " << 3*line_number << std::endl;
	for (int j = 0; j < rows/2; ++j)
	{
		for (int i = 0; i < cols/2; ++i)
		{
			grid << 2 << " " << make_index(i,j,cols/2+1) << " " << make_index(i+1,j,cols/2+1) << std::endl;
			grid << 2 << " " << make_index(i,j,cols/2+1) << " " << make_index(i,j+1,cols/2+1) << std::endl;
		}
		grid << 2 << " " << make_index(cols/2,j,cols/2+1) << " " << make_index(cols/2,j+1,cols/2+1) << std::endl;
	}
	for (int i = 0; i < cols/2; ++i)
	{
		grid << 2 << " " << make_index(i,rows/2,cols/2+1) << " " << make_index(i+1,rows/2,cols/2+1) << std::endl;
	}
	grid.close();

}

void create_assemble(const std::vector<std::vector<float> >& dataVec)
{
	/* get the limit range of four scalar values */
	const int& rows = dataVec.size();
	const int& cols = dataVec[0].size();

	float value;
	int num;
	/* generate the assembled four-scalar normalized for R visualization */
    std::ofstream normalized("assembled", ios::out);
    if(normalized.fail())
    {
    	std::cout << "Error for creating wrong files!" << std::endl;
    	exit(1);
    }

    string rownames[] = {"d_E", "d_E_", "d_F", "d_F_", "d_G", "d_G_", "d_R", "d_R_", "d_M", "d_M_", "d_H", "d_H_", "d_S", "d_S_", "d_P", "d_P_"};
    for (int i = 0; i < cols; ++i)
    {
    	normalized << rownames[i] << "\t";
    	for (int j = 0; j < rows;++j)
    	{
    		if(dataVec[j][i]<=-9999.0)
    		{
    			value = 0.0;
    			normalized << value << "\t";
    			continue;
    		}
    		num = (j%2)*2+i%2;
    		if(num<=1)
	    		value = (dataVec[j][i]-data_range[num][0])/data_range[num][2]*(1.0-range_start)+range_start;
	    	else if(num==2)
	    	{
	    		if(dataVec[j][i]>=max_db_index)
	    			value = 1.0;
	    		else
	    			value = (dataVec[j][i]-data_range[num][0])/data_range[num][2]*(1.0-range_start)+range_start;
	    	}
	    	else if(num==3)
	    		value = (log10(dataVec[j][i])-data_range[num][0])/data_range[num][2]*(1.0-range_start)+range_start;

	    	if(num==0||num==1)
	    		normalized << value << "\t";
	    	else if(num==2||num==3)
	    		normalized << (1.0+range_start)-value << "\t";
		}
		normalized << std::endl;
	}

    normalized.close();
}

void create_separate(const std::vector<std::vector<float> >& dataVec)
{
	/* get the limit range of four scalar values */
	const int& rows = dataVec.size();
	const int& cols = dataVec[0].size();

	float value;
	int num;
	/* generate the assembled four-scalar normalized for R visualization */
    std::ofstream silouette("Silhouette", ios::out), similarity("Gamma", ios::out), dbindex("DBindex", ios::out), validity("Validity", ios::out);
    if(silouette.fail() || similarity.fail() || dbindex.fail() || validity.fail())
    {
    	std::cout << "Error for creating wrong files!" << std::endl;
    	exit(1);
    }

    string metric[] = {"d_E", "d_F", "d_G", "d_R", "d_M", "d_H", "d_S", "d_P"};
    for (int i = 0; i < cols; ++i)
    {
    	if(i%2==0)
    	{
	    	silouette << metric[i/2] << "\t";
			dbindex << metric[i/2] << "\t";
		}
		else
		{
			similarity << metric[i/2] << "\t";
			validity << metric[i/2] << "\t";
		}
    	for (int j = 0; j < rows;++j)
    	{
    		num = (j%2)*2+i%2;
    		if(dataVec[j][i]<=-9999.0)
    		{
    			value = 0.0;
    			switch(num)
    			{
    				case 0:
	    				silouette << value << "\t";
	    				break;

	    			case 1:
	    				similarity << value << "\t";
	    				break;

	    			case 2:
	    				dbindex << value << "\t";
	    				break;

	    			case 3:
	    				validity << value << "\t";
	    				break;
	    		}

    			continue;
    		}
    		if(num<=1)
	    		value = (dataVec[j][i]-data_range[num][0])/data_range[num][2]*(1.0-range_start)+range_start;
	    	else if(num==2)
	    	{
	    		if(dataVec[j][i]>=max_db_index)
	    			value = 1.0;
	    		else
	    			value = (dataVec[j][i]-data_range[num][0])/data_range[num][2]*(1.0-range_start)+range_start;
	    	}
	    	else if(num==3)
	    		value = (log10(dataVec[j][i])-data_range[num][0])/data_range[num][2]*(1.0-range_start)+range_start;

	    	if(num==0||num==1)
	    	{
	    		if(num==0)
	    			silouette << value << "\t";
	    		else if(num==1)
	    			similarity << value << "\t";
	    	}
	    	else if(num==2||num==3)
	    	{
	    		if(num==2)
	    			dbindex << (1.0+range_start)-value << "\t";
	    		else if(num==3)
	    			validity << (1.0+range_start)-value << "\t";
	    	}
		}
		if(i%2==0)
		{
			silouette << std::endl;
			dbindex << std::endl;
		}
		else
		{
			similarity << std::endl;
			validity << std::endl;
		}
	}
	silouette.close();
	similarity.close();
	dbindex.close();
	validity << std::endl;
}

void create_ranking(const std::vector<std::vector<float> >& dataVec)
{
	/* get the limit range of four scalar values */
	const int& rows = dataVec.size();
	const int& cols = dataVec[0].size();

	float value;
	int num;

    unordered_map<int, std::vector<AverageColumn> > rankMap;

    for (int i = 0; i < 4; ++i)
    {
    	rankMap.insert(make_pair(i, std::vector<AverageColumn>()));
    }

    string metric[]={"d_E", "d_F", "d_G", "d_R", "d_M", "d_H", "d_S", "d_P"};

    int effective[2];
    for (int i = 0; i < cols; ++i)
    {
    	AverageColumn ac[2];
    	ac[0].name = ac[1].name = metric[i/2];
    	effective[0] = effective[1] = 0;
    	for (int j = 0; j < rows;j+=2)
    	{
    		num = i%2;
    		ac[0].valueVec.push_back(dataVec[j][i]);
    		if (dataVec[j][i]>-9999.0)
    		{
    			++effective[0];
    			ac[0].average+=dataVec[j][i];
    		}

    		num = i%2+2;
    		ac[1].valueVec.push_back(dataVec[j+1][i]);
    		if (dataVec[j+1][i]>-9999.0)
    		{
    			++effective[1];
    			ac[1].average+=dataVec[j+1][i];
    		}
		}
		ac[0].average/=effective[0];
		ac[1].average/=effective[1];

		rankMap[i%2].push_back(ac[0]);
		rankMap[i%2+2].push_back(ac[1]);
	}

	string clustering[]={"kmeans", "kmedoids", "AHC-single", "AHC-average", "BIRCH", "DBSCAN", "OPTICS", "SC-kmeans", "SC-eigen", "AP", "PCA"};
	unordered_map<int, std::vector<AverageClustering> > clusteringMap;
	for (int j = 0; j < rows; ++j)
	{
		AverageClustering ac[2];
		ac[0].originalIndex = ac[1].originalIndex = j/2;
		effective[0] = effective[1] = 0;
		for (int i = 0; i < cols; ++i)
		{
			if (dataVec[j][i]>-9999.0)
    		{
    			++effective[i%2];
    			ac[i%2].average+=dataVec[j][i];
    		}
		}
		ac[0].average/=effective[0];
		ac[1].average/=effective[1];
		clusteringMap[(j%2)*2].push_back(ac[0]);
		clusteringMap[(j%2)*2+1].push_back(ac[1]);
	}


	for (int i = 0; i < 4; ++i)
	{	
		/* ranking silhouette and gamma from largest to smallest */
		if(i<=1)
		{
			std::sort(rankMap[i].begin(), rankMap[i].end(), [](const AverageColumn& a, const AverageColumn& b)
			{return a.average>b.average;});

			std::sort(clusteringMap[i].begin(), clusteringMap[i].end(), [](const AverageClustering& a, const AverageClustering& b)
			{return a.average>b.average;});
		}
		/* ranking db index and validity from smallest to largest */
		else
		{
			std::sort(rankMap[i].begin(), rankMap[i].end(), [](const AverageColumn& a, const AverageColumn& b)
			{return a.average<b.average;});
			std::sort(clusteringMap[i].begin(), clusteringMap[i].end(), [](const AverageClustering& a, const AverageClustering& b)
			{return a.average<b.average;});
		}
	}
	std::unordered_map<int, std::vector<int> > verticalOrder;
	for (int i = 0; i < 4; ++i)
	{
		verticalOrder[i] = std::vector<int>(clusteringMap[i].size());
		for (int j=0; j<clusteringMap[i].size(); ++j)
		{
			verticalOrder[i][j] = clusteringMap[i][j].originalIndex; 
		}
	}

	float realValue;
	/* generate the assembled four-scalar normalized for R visualization */
    string file_names[] = {"silhouette_ranking", "gamma_ranking", "dbindex_ranking", "validity_ranking"};

    string best_names[] = {"silhouette_best", "gamma_best", "dbindex_best", "validity_best"};

    BestValue best_value[4];
    best_value[0].value = best_value[1].value = -FLT_MAX;
    best_value[2].value = best_value[3].value = FLT_MAX;

    std::cout << "Ranking started..." << std::endl;
    for (int i = 0; i < 4; ++i)
    {
    	std::ofstream colnames((to_string(i)+"_colnames").c_str(), ios::out);
    	if(colnames.fail())
    	{
    		std::cout << "Error for creating files!" << std::endl;
    		exit(1);
    	}

    	std::cout << file_names[i] << std::endl;
    	for (int j = 0; j < verticalOrder[i].size(); ++j)
    	{
    		std::cout << "(" << clustering[verticalOrder[i][j]] << "," << clusteringMap[i][j].average << "), ";;
    		colnames << clustering[verticalOrder[i][j]] << " ";
    	}
    	colnames << std::endl;
    	colnames.close();

    	std::ofstream ranked_file(file_names[i].c_str(), ios::out);
    	if (ranked_file.fail())
    	{
    		std::cout << "Error for creating file!" << std::endl;
    		exit(1);
    	}

    	const std::vector<AverageColumn>& element = rankMap[i];
    	for (int j = 0; j < element.size(); ++j)
    	{
    		ranked_file << element[j].name << "\t";
    		for (int k = 0; k<element[j].valueVec.size(); ++k)
    		{		

    			realValue = (element[j].valueVec)[verticalOrder[i][k]];

    			if(realValue<=-9999.0)
    			{
    				ranked_file << 0.0 << "\t";
    				continue;
    			}
    			else
    			{
    				if (i<=1)
	    			{
	    				if(realValue>=best_value[i].value)
	    				{
	    					best_value[i].value=realValue;
	    					best_value[i].i = k;
	    					best_value[i].j = j;
	    				}
	    			}
	    			else if (i>=2)
	    			{
	    				if(realValue<=best_value[i].value)
	    				{
	    					best_value[i].value=realValue;
	    					best_value[i].i = k;
	    					best_value[i].j = j;
	    				}
	    			}
    			}

    			if(i<=1)
    				value = (realValue-data_range[i][0])/data_range[i][2]*(1.0-range_start)+range_start;
    			else if(i==2)
		    	{
		    		if(realValue>=max_db_index)
		    			value = 1.0;
		    		else
		    			value = (realValue-data_range[i][0])/data_range[i][2]*(1.0-range_start)+range_start;
		    	}
    			else
    				value = (log10(realValue)-data_range[i][0])/data_range[i][2]*(1.0-range_start)+range_start;

    			if(i<=1)
    				ranked_file << value << "\t";
    			else
    				ranked_file << (1.0+range_start)-value << "\t";
    		}
    		ranked_file << std::endl;
    	}
    	std::cout << std::endl;
    	ranked_file.close();
    }

    std::cout << std::endl;
    for (int i = 0; i < 4; ++i)
    {
    	std::ofstream best_file(best_names[i].c_str(), ios::out);
    	if(best_file.fail())
    	{
    		std::cout << "Error for creating files!" << std::endl;
    		exit(1);
    	}

    	const std::vector<AverageColumn>& element = rankMap[i];
    	std::cout << best_value[i].value << std::endl;
    	for (int j = 0; j < element.size(); ++j)
    	{
    		best_file << element[j].name << "\t";
    		for (int k = 0; k < element[j].valueVec.size(); ++k)
    		{
    			if (best_value[i].i==k && best_value[i].j==j)
    				best_file << 1.0 << "\t";
    			else
    				best_file << -1.0 << "\t";
    		}
    		best_file << std::endl;
    	}
    	best_file.close();
    }
}


void create_latex_table(const std::vector<std::vector<float> >& dataVec)
{
	const int& rows = dataVec.size();
	const int& cols = dataVec[0].size();
	std::ofstream latex_table("latex_table", ios::out);
	if(latex_table.fail())
	{
		std::cout << "Error for creating latex table w.r.t. data vec!" << std::endl;
		exit(1);
	}

	/* get the best respective value */
	string clustering[] = {"\\textbf{kmeans}", "\\textbf{kmedoids}", "\\textbf{AHC}-single", "\\textbf{AHC}-average", "\\textbf{BIRCH}", 
							"\\textbf{DBSCAN}", "\\textbf{OPTICS}", "\\textbf{SC}-kmeans", "\\textbf{SC}-eigen", "\\textbf{AP}", "\\textbf{PCA}"};
	BestValue ***best_value = new BestValue**[cols/2];
	for (int i = 0; i < cols/2; ++i)
	{
		best_value[i] = new BestValue*[2];
		for(int j=0; j<2; ++j)
		{
			if(j==0)
			{
				best_value[i][j] = new BestValue[2];
				best_value[i][j][0].value = -10000.0;
				best_value[i][j][1].value = -10000.0;
			}
			else if(j==1)
			{
				best_value[i][j] = new BestValue[2];
				best_value[i][j][0].value = FLT_MAX;
				best_value[i][j][1].value = FLT_MAX;
			}
		}
	}

	int col_num, index;
	for (int i = 0; i < cols; ++i)
	{
		col_num = i/2;
		index = i%2;
		for (int j = 0; j < rows; j+=2)
		{
			if (dataVec[j][i]<=-9999.0)
				continue;

			if (dataVec[j][i]>best_value[col_num][0][index].value)
			{
				best_value[col_num][0][index].value = dataVec[j][i];
				best_value[col_num][0][index].i = i;
				best_value[col_num][0][index].j = j;
			}
		}

		for (int j = 1; j < rows; j+=2)
		{
			if (dataVec[j][i]<=-9999.0)
				continue;

			if (dataVec[j][i]<best_value[col_num][1][index].value)
			{
				best_value[col_num][1][index].value = dataVec[j][i];
				best_value[col_num][1][index].i = i;
				best_value[col_num][1][index].j = j;
			}
		}
	}


	string tag[] = {"\\dashuline", "\\underline", "*", "\\textbf"};

	for (int j = 0; j<rows; ++j)
	{
		latex_table << "\\multirow{2}{*}{" << clustering[j/2] << "} ";
		for (int i = 0; i < cols; ++i)
		{
			latex_table << "&";
			if (i==best_value[i/2][0][i%2].i && j==best_value[i/2][0][i%2].j)
				latex_table << tag[i%2] << "{";

			if (dataVec[j][i]<=-9999.0)
				latex_table << "-";
			else
				latex_table << std::fixed << std::setprecision(3) << dataVec[j][i];
			if (i==best_value[i/2][0][i%2].i && j==best_value[i/2][0][i%2].j)
				latex_table << "}";
		}
		latex_table << "\\\\" << std::endl;

		j+=1;
		for (int i = 0; i < cols; ++i)
		{
			if(i%2==0)
			{
				latex_table << "& \\multicolumn{1}{l}{";

				if (dataVec[j][i]<=-9999.0)
					latex_table << "\\hspace{0.23cm}" << "-";
				else if (dataVec[j][i]>=0.01 && dataVec[j][i]<1000)
					latex_table << std::fixed << std::setprecision(3) << dataVec[j][i];
				else
					latex_table << std::scientific << std::setprecision(1) << dataVec[j][i];
				if (i==best_value[i/2][1][0].i && j==best_value[i/2][1][0].j)
					latex_table << tag[2];

				latex_table << "}";
			}
			else
			{
				latex_table << "&";
				if (i==best_value[i/2][1][1].i && j==best_value[i/2][1][1].j)
					latex_table << tag[3] << "{";

				if (dataVec[j][i]<=-9999.0)
					latex_table << "-";
				else if (dataVec[j][i]>=0.01)
					latex_table << std::fixed << std::setprecision(3) << dataVec[j][i];
				else
					latex_table << std::scientific << std::setprecision(1) << dataVec[j][i];
				if (i==best_value[i/2][1][1].i && j==best_value[i/2][1][1].j)
					latex_table << "}";
			}
		}
		latex_table << "\\\\" << std::endl;

		latex_table << "\\hline" << std::endl;
	}

	latex_table.close();

}


void get_average_value(std::vector<std::vector<float> >& averageValue, std::vector<std::vector<float> >& standardDeviation, 
				  string filenames[], const int& file_size)
{
	std::vector<std::vector<float> > totalValue[file_size];

	for (int i = 0; i < file_size; ++i)
	{
		readData(totalValue[i], filenames[i].c_str());
	}

	const int& rows = totalValue[0].size();
	const int& cols = totalValue[0][0].size();

	averageValue = std::vector< std::vector<float> >(rows, std::vector<float>(cols));

	if(file_size>=4)
		standardDeviation = std::vector< std::vector<float> >(rows, std::vector<float>(cols));

	float average, stdeviation, value;
	int effective;
	for (int i = 0; i < cols; ++i)
	{
		for (int j = 0; j < rows; ++j)
		{
			average = stdeviation = 0.0;
			effective = 0;
			for (int k = 0; k < file_size; ++k)
			{
				value = totalValue[k][j][i];
				if (value<=-9999.0)
					continue;
				average+=value;
				if(file_size>=4)
					stdeviation+=value*value;
				++effective;
			}
			if(effective==0)
			{
				averageValue[j][i] = -10000.0;
				if(file_size>=4)
					standardDeviation[j][i] = -10000.0;
			}
			else
			{
				average/=float(effective);
				averageValue[j][i] = average;
				if(file_size>=4)
				{
					stdeviation = stdeviation/float(effective)-average*average;
					if(stdeviation<0)
					{
						std::cout << "Error for one-pass standard deviation computation!" << std::endl;
						exit(1);
					}
					standardDeviation[j][i] = sqrt(stdeviation);
				}
			}
		}
	}

}


void create_std_ranking(const std::vector<std::vector<float> >& dataVec, const std::vector<std::vector<float> >& standardDeviation,
	                    const int& file_size)
{
	/* get the limit range of four scalar values */
	const int& rows = dataVec.size();
	const int& cols = dataVec[0].size();

	float value;
	int num;

	float value_option[] = {FLT_MAX, -FLT_MAX, 0};
	for (int i = 0; i < 3; ++i)
	{
		for (int j=0;j<4;++j)
			data_range[j][i] = value_option[i];
	}

	for (int i = 0; i < rows; ++i)
	{
		for (int j = 0; j < cols; ++j)
		{
			if(dataVec[i][j]<=-9999.0)
				continue;
			num = (i%2)*2+j%2;
			data_range[num][0] = std::min(data_range[num][0], dataVec[i][j]);
			data_range[num][1] = std::max(data_range[num][1], dataVec[i][j]);

		}
	}
	for (int i = 0; i < 4; ++i)	
	{
		std::cout << data_range[i][0] << " " << data_range[i][1] << std::endl;
		if(i==2)
		{
			if(data_range[i][1]>=max_db_index)
				data_range[i][1]=max_db_index;
		}
		else if(i==3)
		{
			data_range[i][0] = log10(data_range[i][0]);
			data_range[i][1] = log10(data_range[i][1]);
		}
		data_range[i][2] = data_range[i][1]-data_range[i][0];
	}


    unordered_map<int, std::vector<AverageColumn> > rankMap;

    for (int i = 0; i < 4; ++i)
    {
    	rankMap.insert(make_pair(i, std::vector<AverageColumn>()));
    }

    string metric[]={"d_E", "d_F", "d_G", "d_R", "d_M", "d_H", "d_S", "d_P"};

    int effective[2];
    for (int i = 0; i < cols; ++i)
    {
    	AverageColumn ac[2];
    	ac[0].name = ac[1].name = metric[i/2];
    	effective[0] = effective[1] = 0;
    	for (int j = 0; j < rows;j+=2)
    	{
    		num = i%2;
    		ac[0].valueVec.push_back(dataVec[j][i]);
    		if(file_size>=4)
    			ac[0].std_vec.push_back(standardDeviation[j][i]);
    		if (dataVec[j][i]>-9999.0)
    		{
    			++effective[0];
    			ac[0].average+=dataVec[j][i];
    		}

    		num = i%2+2;
    		ac[1].valueVec.push_back(dataVec[j+1][i]);
    		if(file_size>=4)
    			ac[1].std_vec.push_back(standardDeviation[j+1][i]);
    		if (dataVec[j+1][i]>-9999.0)
    		{
    			++effective[1];
    			ac[1].average+=dataVec[j+1][i];
    		}
		}
		ac[0].average/=effective[0];
		ac[1].average/=effective[1];

		rankMap[i%2].push_back(ac[0]);
		rankMap[i%2+2].push_back(ac[1]);
	}

	string clustering[]={"kmeans", "kmedoids", "AHC-single", "AHC-average", "BIRCH", "DBSCAN", "OPTICS", "SC-kmeans", "SC-eigen", "AP", "PCA"};
	unordered_map<int, std::vector<AverageClustering> > clusteringMap;
	for (int j = 0; j < rows; ++j)
	{
		AverageClustering ac[2];
		ac[0].originalIndex = ac[1].originalIndex = j/2;
		effective[0] = effective[1] = 0;
		for (int i = 0; i < cols; ++i)
		{
			if (dataVec[j][i]>-9999.0)
    		{
    			++effective[i%2];
    			ac[i%2].average+=dataVec[j][i];
    		}
		}
		ac[0].average/=effective[0];
		ac[1].average/=effective[1];
		clusteringMap[(j%2)*2].push_back(ac[0]);
		clusteringMap[(j%2)*2+1].push_back(ac[1]);
	}

	for (int i = 0; i < 4; ++i)
	{	
		/* ranking silhouette and gamma from largest to smallest */
		if(i<=1)
		{
			std::sort(rankMap[i].begin(), rankMap[i].end(), [](const AverageColumn& a, const AverageColumn& b)
			{return a.average>b.average;});

			std::sort(clusteringMap[i].begin(), clusteringMap[i].end(), [](const AverageClustering& a, const AverageClustering& b)
			{return a.average>b.average;});
		}
		/* ranking db index and validity from smallest to largest */
		else
		{
			std::sort(rankMap[i].begin(), rankMap[i].end(), [](const AverageColumn& a, const AverageColumn& b)
			{return a.average<b.average;});
			std::sort(clusteringMap[i].begin(), clusteringMap[i].end(), [](const AverageClustering& a, const AverageClustering& b)
			{return a.average<b.average;});
		}
	}
	std::unordered_map<int, std::vector<int> > verticalOrder;
	for (int i = 0; i < 4; ++i)
	{
		verticalOrder[i] = std::vector<int>(clusteringMap[i].size());
		for (int j=0; j<clusteringMap[i].size(); ++j)
		{
			verticalOrder[i][j] = clusteringMap[i][j].originalIndex; 
		}
	}

	float realValue;
	/* generate the assembled four-scalar normalized for R visualization */
    string file_names[] = {"silhouette_ranking", "gamma_ranking", "dbindex_ranking", "validity_ranking"};

    string best_names[] = {"silhouette_best", "gamma_best", "dbindex_best", "validity_best"};

    BestValue best_value[4];
    best_value[0].value = best_value[1].value = -FLT_MAX;
    best_value[2].value = best_value[3].value = FLT_MAX;

    for (int i = 0; i < 4; ++i)
    {
    	std::ofstream colnames((to_string(i)+"_colnames").c_str(), ios::out);
    	if(colnames.fail())
    	{
    		std::cout << "Error for creating files!" << std::endl;
    		exit(1);
    	}
    	for (int j = 0; j < verticalOrder[i].size(); ++j)
    	{
    		colnames << clustering[verticalOrder[i][j]] << " ";
    	}
    	colnames << std::endl;
    	colnames.close();

    	std::ofstream ranked_file(file_names[i].c_str(), ios::out);
    	if (ranked_file.fail())
    	{
    		std::cout << "Error for creating file!" << std::endl;
    		exit(1);
    	}

    	const std::vector<AverageColumn>& element = rankMap[i];
    	for (int j = 0; j < element.size(); ++j)
    	{
    		ranked_file << element[j].name << "\t";
    		for (int k = 0; k<element[j].valueVec.size(); ++k)
    		{
    			realValue = (element[j].valueVec)[verticalOrder[i][k]];
    			if(realValue<=-9999.0)
    			{
    				ranked_file << 0.0 << "\t";
    				continue;
    			}
    			else
    			{
    				if (i<=1)
	    			{
	    				if(realValue>=best_value[i].value)
	    				{
	    					best_value[i].value=realValue;
	    					best_value[i].i = k;
	    					best_value[i].j = j;
	    				}
	    			}
	    			else if (i>=2)
	    			{
	    				if(realValue<=best_value[i].value)
	    				{
	    					best_value[i].value=realValue;
	    					best_value[i].i = k;
	    					best_value[i].j = j;
	    				}
	    			}
    			}

    			if(i<=1)
    				value = (realValue-data_range[i][0])/data_range[i][2]*(1.0-range_start)+range_start;
    			if(i==2)
				{
					if(realValue>=max_db_index)
						value = 1.0;
					else
						value = (realValue-data_range[i][0])/data_range[i][2]*(1.0-range_start)+range_start;
				}
    			else if(i==3)
    				value = (log10(realValue)-data_range[i][0])/data_range[i][2]*(1.0-range_start)+range_start;

    			if(i<=1)
    				ranked_file << value << "\t";
    			else
    				ranked_file << (1.0+range_start)-value << "\t";
    		}
    		ranked_file << std::endl;
    	}
    	ranked_file.close();
    }


    std::cout << std::endl;
    for (int i = 0; i < 4; ++i)
    {
    	std::ofstream best_file(best_names[i].c_str(), ios::out);
    	if(best_file.fail())
    	{
    		std::cout << "Error for creating files!" << std::endl;
    		exit(1);
    	}

    	const std::vector<AverageColumn>& element = rankMap[i];
    	for (int j = 0; j < element.size(); ++j)
    	{
    		best_file << element[j].name << "\t";
    		for (int k = 0; k < element[j].valueVec.size(); ++k)
    		{
    			if (best_value[i].i==k && best_value[i].j==j)
    				best_file << 1.0 << "\t";
    			else
    				best_file << -1.0 << "\t";
    		}
    		best_file << std::endl;
    	}
    	best_file.close();
    }

    if(file_size>=4)
    {

	    float std_range[4][3]=
		{
			FLT_MAX, -FLT_MAX, 0,
			FLT_MAX, -FLT_MAX, 0,
			FLT_MAX, -FLT_MAX, 0,
			FLT_MAX, -FLT_MAX, 0
		};

		for (int i = 0; i < rows; ++i)
		{
			for (int j = 0; j < cols; ++j)
			{
				if(standardDeviation[i][j]<=-9999.0)
					continue;
				num = (i%2)*2+j%2;
				std_range[num][0] = std::min(std_range[num][0], standardDeviation[i][j]);
				std_range[num][1] = std::max(std_range[num][1], standardDeviation[i][j]);

			}
		}
		for (int i = 0; i < 4; ++i)	
		{
			std_range[i][2] = std_range[i][1]-std_range[i][0];
			std::cout << "[" << std_range[i][0] << "," << std_range[i][1] << "]: " << std_range[i][2] << std::endl;
		}


	    string std_names[] = {"silhouette_std", "gamma_std", "dbindex_std", "validity_std"};

	    for (int i = 0; i < 4; ++i)
	    {
	    	std::ofstream ranked_file(std_names[i].c_str(), ios::out);
	    	if (ranked_file.fail())
	    	{
	    		std::cout << "Error for creating file!" << std::endl;
	    		exit(1);
	    	}

	    	const std::vector<AverageColumn>& element = rankMap[i];
	    	for (int j = 0; j < element.size(); ++j)
	    	{
	    		ranked_file << element[j].name << "\t";
	    		for (int k = 0; k<element[j].std_vec.size(); ++k)
	    		{
	    			realValue = (element[j].std_vec)[verticalOrder[i][k]];
	    			if(realValue<=-9999.0)
	    			{
	    				ranked_file << 0.0 << "\t";
	    				continue;
	    			}
	    			value = (realValue-std_range[i][0])/std_range[i][2]*(1.0-range_start)+range_start;

	    			ranked_file << value << "\t";
	    		}
	    		ranked_file << std::endl;
	    	}
	    	ranked_file.close();
	    }
	}

}