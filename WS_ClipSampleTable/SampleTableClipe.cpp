#include "SampleTableClipe.h"
#include<fstream>
#include<cstdlib>
using namespace std;

SampleTableClipe::SampleTableClipe(string sampleTableAndExtent,string clippedSampleTable)
{
	this->clippedSampleTableFn = clippedSampleTable;
	this->noExtentConstraintFlag = -9999;
	//分割字符串
	/*const char *separator = ";" ;
	int length = sampleTableAndExtent.length();
	char *m = new char[length + 1];
	strcpy(m, sampleTableAndExtent.c_str());
	char *p = new char[length + 1];
	p = strtok(p, separator);
	vector<string> tmp;
	while(p)
	{
		tmp.push_back(p);
		p = strtok(NULL,separator);
	}*/
	vector<string> tmp;
	this->sampleTable->parseStr(sampleTableAndExtent, ';', tmp); 



	this->sampleTableFn = tmp[0];
	this->sampleTable = new Table(sampleTableFn);

	if(tmp.size() == 2)
	{
		vector<string> extent;
		//分割字符串，与上面是两种方法
		this->sampleTable->parseStr(tmp[1], ' ', extent);
		//sampleTable.parseStr(tmp[1], ' ', extent);
		vector<double> clipExtent;
		for (int i = 0; i < 4; i++)
		{
			this->clipExtent.push_back(atof(extent[i].c_str()));
		}
	}
	else
	{
		
		for (int i = 0; i < 4; i++)
		{
			this->clipExtent.push_back(noExtentConstraintFlag);
		}
	}
}


void SampleTableClipe::ClipSampleTableByExtent()
{
	vector<string> header = this->sampleTable->GetHeader();
	vector<vector<string> > values = this->sampleTable->GetValues();
	vector<vector<string> > clippedValues;
	if (clipExtent[0] != noExtentConstraintFlag && clipExtent[1] != noExtentConstraintFlag && clipExtent[2] != noExtentConstraintFlag && clipExtent[3] != noExtentConstraintFlag)
	{
		double left = clipExtent[0];
		double bottom = clipExtent[1];
		double right = clipExtent[2];
		double top = clipExtent[3];
		int xIndex = this->sampleTable->GetColumnIndexByName("X");
		int yIndex = this->sampleTable->GetColumnIndexByName("Y");
		for (unsigned int i = 0; i < values.size(); i++)
		{
			double record_x = atof(this->sampleTable->GetValue(xIndex, i).c_str());
			double record_y = atof(this->sampleTable->GetValue(yIndex, i).c_str());
			if (record_x >= left && record_x <= right && record_y >= bottom && record_y <= top)
			{
				clippedValues.push_back(values[i]);
			}
		}
	}
	else
	{
		clippedValues = values;
	}
	Table *clippedTable = new Table(header, clippedValues);
	clippedTable->WriteAsCSVFile(clippedSampleTableFn);
}