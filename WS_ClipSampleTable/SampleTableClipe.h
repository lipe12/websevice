#ifndef SampleTableClipe_H
#define SampleTableClipe_H
#include "Table.h"
#include<iostream>
#include<string>
#include<fstream>
#include<list>
#include<vector>
using namespace std;

class SampleTableClipe
{
public:
	string sampleTableFn;
	Table *sampleTable;
	vector<double> clipExtent;
	double noExtentConstraintFlag;
	string clippedSampleTableFn;
	//init
	SampleTableClipe(string sampleTableAndExtent,string clippedSampleTable);
	//method
	void ClipSampleTableByExtent();
};
#endif