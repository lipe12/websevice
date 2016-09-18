#ifndef SAMPLEINTEGRATION_H
#define SAMPLEINTEGRATION_H

#define VERY_SMALL 0.000001
#include <string>
#include <vector>
#include <math.h>
#include <limits>
using namespace std;

class SampleIntegration {
public:
	SampleIntegration();
	SampleIntegration(string s);
	void SetIntegrationMethod(string integrationMethod);
	string GetIntegrationMethod();
	double GetSampleSimilarity(double* SimilaritiesOfCategories, int count, int* LayersInCategories);
	void GetSampleSimi(double** & result, double** climatesimi, double** geogloySimi, 
		double** terrainSimi, double** vegeSimi, double** otherSimi, int row, int col, 
		int climateLyrCnt, int geologyLyrCnt, int terrainLyrCnt, int VegeLyrCnt, int otherLyrCnt);
private:
	string integrationMethod;
};


#endif
