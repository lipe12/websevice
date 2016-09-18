#ifndef CATEGORYINTEGRATION_H
#define CATEGORYINTEGRATION_H

#include <string>
#include <vector>
using namespace std;

class CategoryIntegration
{
public:
	CategoryIntegration();
	CategoryIntegration(string s);
	void setIntegrationMethod(string IntegrationMethod);
	string getIntegrationMethod();
	void getCategorySimilarity(vector<double**>& simiArray, double** &result, int row_size, int col_size);
	double GetCategorySimilarity(double* SimilaritiesOfAttributes, int n);
private:
	string IntegrationMethod;
};

#endif