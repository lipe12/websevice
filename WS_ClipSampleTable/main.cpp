#include<iostream>
#include<vector>
#include "Table.h"
#include "SampleTableClipe.h"
using namespace std;


int main(int argc, char *argv[])
{
	
	char *sampleTableAndExtent = argv[1];
	char *clippedSampleTable = argv[2];
	string st(sampleTableAndExtent);
	string cs(clippedSampleTable);
	//string st = "E:\\C++pro\\xuancheng_Layer_B.CSV;118.717 30.762 119.358 31.126";
	//string cs = "E:\\C++pro\\xuancheng.CSV";
	SampleTableClipe tablecliper(st, cs);
	tablecliper.ClipSampleTableByExtent();
	return 0;
}