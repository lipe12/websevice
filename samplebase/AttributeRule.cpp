//Calculate similarity based on each kind of environmental variables
//climate similarity
//terrain similarity
//geology similarity
//vegetation similarity

#include <iostream>
#include "AttributeRule.h"
#include <stdio.h>
AttributeRule::AttributeRule(){

}

void AttributeRule::setCategory(std::string Category){
	this->Category = Category;
}

void AttributeRule::setDistCalculationMethod(std::string method){
	this->DistCalculationMethod = method;
}

std::string AttributeRule::getCategory(){
	return this->Category;
}

std::string AttributeRule::getDistCalculationMethod(){
	return this->DistCalculationMethod;
}

double AttributeRule::getAttributeSimilarity(double AttriVal, double sampleVal, double range, double sum, int counter, double dataStd){
	if (DistCalculationMethod == "Gower"){
		double i = (double) 1 - fabs(AttriVal - sampleVal) / range;
        return i;

    }else if (DistCalculationMethod=="Boolean"){
		if (fabs(sampleVal - AttriVal) <= VERY_SMALL)
			return 1;
		else
            return 0;
	}else if(DistCalculationMethod == "Gaussian")
	{
        double i = exp(-pow(AttriVal - sampleVal, 2)/(2 * pow(dataStd, 4) / (sum/counter)));
        return i;
	}else{
		return 0;
	}
}
