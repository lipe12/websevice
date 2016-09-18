#include "SampleIntegration.h"
#include <cmath>

SampleIntegration::SampleIntegration(){
}

SampleIntegration::SampleIntegration(string s):integrationMethod(s)
{
}


void SampleIntegration::SetIntegrationMethod(string integrationMethod){
	this->integrationMethod = integrationMethod;
}

string SampleIntegration::GetIntegrationMethod(){
	return this->integrationMethod;
}

//此函数用来综合某个样点与某个待推测点的多个类别的相似度（比如，综合气候，母质，地形等相似度得到一个样点与一个待推测点的相似度）
double SampleIntegration::GetSampleSimilarity(double* SimilaritiesOfCategories, int count, int* LayersInCategories)
{
    if (integrationMethod == "Average")
    {
		//LayersInCategories[1]表示输入的geology类别下图层的个数
        if (LayersInCategories[1] > 0) //geology will be used
        {
                    
            if (abs(SimilaritiesOfCategories[1]) > -99999)//SimilaritiesOfCategories[1]是母质类别的相似度,为何不直接=1呢？
            {
                double Sum = 0;
                int EffectiveCategory = 0;
                for (int i = 0; i < count; i++)
                {
                    if (i != 1 && LayersInCategories[i] > 0)
                    {
                        Sum += SimilaritiesOfCategories[i];
                        EffectiveCategory++;
                    }
                }
                return Sum / EffectiveCategory;
            }
            else
                return 0;
                    
        }
        else
        {
            double Sum = 0;
            int EffectiveCategory = 0;
            for (int i = 0; i < count; i++)
            {
                if (i != 1 && LayersInCategories[i] > 0)
                {
                    Sum += SimilaritiesOfCategories[i];
                    EffectiveCategory++;
                }
            }
            return Sum / EffectiveCategory;
        }
    }
    else
    {
               
        double Min = 99999;
        for (int i = 0; i < count; i++)
        {
            //if (SimilaritiesOfCategories[i] <= Min && LayersInCategories[i] > 0)
            if (SimilaritiesOfCategories[i] < Min && LayersInCategories[i] > 0) // by jjc
            {
                Min = SimilaritiesOfCategories[i];
            }
        }

        return Min;
    }
       
}

void SampleIntegration::GetSampleSimi(double** & result, double** climateSimi, double** geologySimi, double** terrainSimi, 
	double** vegeSimi, double** otherSimi, int row_size, int col_size, 
	int climateLyrCnt, int geologyLyrCnt, int terrainLyrCnt, int vegeLyrCnt, int otherLyrCnt){
	if(integrationMethod == "Average"){
		if (geologyLyrCnt > 0){
			for(int colIdx = 0; colIdx < col_size; colIdx ++){
				for(int rowIdx = 0; rowIdx < row_size; rowIdx ++){
					if(geologySimi[rowIdx][colIdx]==1){
						double sumSimi = 0;
						int effectiveCategory = 0;
						if(climateLyrCnt > 0){
							sumSimi += climateSimi[rowIdx][colIdx];
							effectiveCategory++;
						}
						if(terrainLyrCnt > 0){
							sumSimi += terrainSimi[rowIdx][colIdx];
							effectiveCategory++;
						}
						if(vegeLyrCnt > 0){
							sumSimi += vegeSimi[rowIdx][colIdx];
							effectiveCategory++;
						}
						if(otherLyrCnt > 0){
							sumSimi += otherSimi[rowIdx][colIdx];
							effectiveCategory++;
						}
						result[rowIdx][colIdx] = sumSimi / effectiveCategory;
					}else{
						result[rowIdx][colIdx] = 0;
					}
				}
			}
		}else{
			int effectiveCategory = 0;
			if(climateLyrCnt > 0)
				effectiveCategory ++;
			if(terrainLyrCnt > 0)
				effectiveCategory ++;
			if(vegeLyrCnt > 0)
				effectiveCategory ++;
			if(otherLyrCnt > 0)
				effectiveCategory ++;

			for(int colIdx = 0; colIdx < col_size; colIdx ++){
				for(int rowIdx = 0; rowIdx < row_size; rowIdx ++){
					double sumSimi = 0;
					if(climateLyrCnt > 0){
						sumSimi += climateSimi[rowIdx][colIdx];
					}
					if(terrainLyrCnt > 0){
						sumSimi += terrainSimi[rowIdx][colIdx];
					}
					if(vegeLyrCnt > 0){
						sumSimi += vegeSimi[rowIdx][colIdx];
					}
					if(otherLyrCnt > 0){
						sumSimi += otherSimi[rowIdx][colIdx];
					}
					result[rowIdx][colIdx] = sumSimi / effectiveCategory;
				}
			}
		}
	}else{
		for(int colIdx = 0; colIdx < col_size; colIdx ++){
			for(int rowIdx = 0; rowIdx < row_size; rowIdx ++){
				double minSimi = 9999;
				if(climateLyrCnt > 0)
					minSimi = min(minSimi,climateSimi[rowIdx][colIdx]);
				if(geologyLyrCnt > 0)
					minSimi = min(minSimi,geologySimi[rowIdx][colIdx]);
				if(terrainLyrCnt > 0)
					minSimi = min(minSimi,terrainSimi[rowIdx][colIdx]);
				if(vegeLyrCnt > 0)
					minSimi = min(minSimi,vegeSimi[rowIdx][colIdx]);
				if(otherLyrCnt > 0)
					minSimi = min(minSimi,otherSimi[rowIdx][colIdx]);
				result[rowIdx][colIdx] = minSimi;
			}
		}
	}
}
	
