#ifndef ATTRIBUTERULE_H
#define ATTRIBUTERULE_H

#define VERY_SMALL 0.000001

#include <string>
#include <cmath>
using namespace std;

class AttributeRule
{
public:
	AttributeRule();
	void setCategory(std::string Category);
	void setDistCalculationMethod(std::string DistCalculationMethod);
	std::string getCategory();
	std::string getDistCalculationMethod();
	double getAttributeSimilarity(double AttriVal, double sampleVal, double range, double sum, int counter, double dataStd);//, double sum, int counter, double dataStd
private:
	std::string Category;
	std::string DistCalculationMethod;
};


#endif

