#include "CategoryIntegration.h"
#include <math.h>
#include <iostream>
#include <stdio.h>
CategoryIntegration::CategoryIntegration(){
}

CategoryIntegration::CategoryIntegration(string s){
	this->IntegrationMethod = s;
}
void CategoryIntegration::setIntegrationMethod(string method){
	this->IntegrationMethod = method;
}

string CategoryIntegration::getIntegrationMethod(){
	return this->IntegrationMethod;
}

double CategoryIntegration::GetCategorySimilarity(double* SimilaritiesOfAttributes, int n)
{
    if (n == 0)//n是每个类别下环境因子图层的个数
        return 0;

    if (IntegrationMethod == "Average")
    {
        double Sum = 0;
        for (int i = 0; i < n; i++)
        {
            Sum += SimilaritiesOfAttributes[i];

        }
        return Sum / n;
    }
    else
    {

        double Min = 99999;
        for (int i = 0; i < n; i++)
        {
            if (SimilaritiesOfAttributes[i] <= Min)
            //if (SimilaritiesOfAttributes[i] < Min)// by jjc
            {
                Min = SimilaritiesOfAttributes[i];
            }
        }
        return Min;
    }
}

void CategoryIntegration::getCategorySimilarity(vector<double**>& simiArray, double** & result, int row_size, int col_size){

	if (simiArray.size()==0){
		for(int rowIdx = 0; rowIdx < row_size; rowIdx ++){
			for(int colIdx = 0; colIdx < col_size; colIdx ++){
				result[rowIdx][colIdx] = 0;
			}
		}
	}else if(this->IntegrationMethod.compare("Average") == 0){

		for(int rowIdx = 0; rowIdx < row_size; rowIdx ++){
			for(int colIdx = 0; colIdx < col_size; colIdx ++){
				double sum = 0;
				bool noData = false;
				for (int sampleIdx = 0; !noData && sampleIdx < simiArray.size(); sampleIdx ++){
					double v = simiArray[sampleIdx][rowIdx][colIdx];
					//if(v == -9999){
                    if(v < -100){
						noData = true;
					}else{
						sum += v;
					}
				}
				if(noData)
					result[rowIdx][colIdx]=-9999;
				else
					result[rowIdx][colIdx] = sum / simiArray.size();
			}
		}
	}else{
		for(int rowIdx = 0; rowIdx < row_size; rowIdx ++){
			for(int colIdx = 0; colIdx < col_size; colIdx ++){
				double min = 9999;
				bool noData = false;
				for (int sampleIdx = 0; !noData && sampleIdx < simiArray.size(); sampleIdx ++){
					double v = (simiArray[sampleIdx])[rowIdx][colIdx];
					//if(v == -9999){
                    if(v< -100){
						noData = true;
					}else{
						if(v<min)
							min = v;
					}
				}
				if(noData)
					result[rowIdx][colIdx]=-9999;
				else
					result[rowIdx][colIdx] = min;
			}
		}
	}
}
