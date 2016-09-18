#include <sstream>
#include "ParallelSBMp.h"
#include <iostream>
#include <math.h>

#include <fstream>
#include <iomanip>

using namespace std;

double ParallelSBMp::string_to_double( const std::string& s )
{
   std::istringstream i(s);
   double x;
   if (!(i >> x))
     return 0;
   return x;
 }

ParallelSBMp::ParallelSBMp(){

}

ParallelSBMp::ParallelSBMp(double * sampleValue){

   

}

ParallelSBMp::ParallelSBMp(double * sampleValue,double * sampleClimateV, double * sampleGeologyV, double * sampleTerrainV,
             double * sampleVegetationV, double * sampleOtherV,
             vector<string> &AttriRules,
             double * climateVRange,double * geologyVRange,double * terrainVRange,double * vegeVRange,double * otherVRange,
             string catIntegrationMethod, string sampleIntegrationMethod,
			 string uncertaintyThreshold
			 ):attriRules(AttriRules.size()),catIntegration(catIntegrationMethod),sampleIntegration(sampleIntegrationMethod)
{

	 this->sampleValue = sampleValue;

     this->sampleClimateV = sampleClimateV;
     this->sampleGeologyV = sampleGeologyV;
     this->sampleTerrainV = sampleTerrainV;
	 this->sampleVegetationV = sampleVegetationV;
	 this->sampleOtherV = sampleOtherV;

	 this->climateVRange = climateVRange;
	 this->geologyVRange = geologyVRange;
	 this->terrainVRange = terrainVRange;
	 this->vegeVRange = vegeVRange;
	 this->otherVRange = otherVRange;

	 for (int i = 0; i < AttriRules.size(); i ++){
			vector<string> ruleParameterArray;
			split(AttriRules[i],"?",ruleParameterArray);
			this->attriRules[i].setCategory(ruleParameterArray[0]);
			this->attriRules[i].setDistCalculationMethod(ruleParameterArray[1]);
	}
	this->uncertaintyThreshold = atof(uncertaintyThreshold.c_str());


}

void ParallelSBMp::getPropertyMap(double *climateStd,double *geoStd,double *terrainStd,double *vegetationStd,double *otherStd,
					double ** climate2Dvalues, double ** geo2Dvalues,double ** terrain2Dvalues, double ** vege2Dvalues, double ** other2Dvalues,
				   double ** propertyVals, double ** uncertaintyVals,int numOfSamples,int ClimateLyrCnt, int GeogLyrCnt,
				    int TerrainLyrCnt, int VegeLyrCnt, int OtherLyrCnt,int Rows, int Cols, double noData,double *climateAllNoData,
				   double *geoAllNoData,double *terrainAllNoData,double *vegetationAllNoData,double *otherAllNoData,
				  double **climatesumAll, double **geologysumAll, double **terrainsumAll, double **vegetationsumAll, double **othersumAll,
				int *climatecounterAll, int *geologycounterAll, int *terraincounterAll, int *vegetationcounterAll, int *othercounterAll)//wtf
   {
	   AttributeRule* climateRules = new AttributeRule[ClimateLyrCnt];
		AttributeRule* geologyRules = new AttributeRule[GeogLyrCnt];
		AttributeRule* terrainRules = new AttributeRule[TerrainLyrCnt];
		AttributeRule* vegetationRules = new AttributeRule[VegeLyrCnt];
		AttributeRule* otherRules = new AttributeRule[OtherLyrCnt];
		int  numOfLyrs = ClimateLyrCnt + GeogLyrCnt + TerrainLyrCnt + VegeLyrCnt + OtherLyrCnt;
		int climateCnt =0;
		int geologyCnt = 0;
		int terrainCnt =0;
		int vegeCnt =0;
		int otherCnt =0;       
		for(int attriDataLyrIdx = 0; attriDataLyrIdx <numOfLyrs; attriDataLyrIdx++){

			string category = attriRules[attriDataLyrIdx].getCategory();

			if(category == "Climate"){
				 climateRules[climateCnt] = attriRules[attriDataLyrIdx];
				 climateCnt ++;
			}else if(category == "Geology"){

				  geologyRules[geologyCnt] = attriRules[attriDataLyrIdx];
				  geologyCnt ++;
			}else if(category == "Terrain"){
				  terrainRules[terrainCnt] = attriRules[attriDataLyrIdx];
				  terrainCnt ++;
			}else if(category == "Vegetation"){
				  vegetationRules[vegeCnt] = attriRules[attriDataLyrIdx];
				  vegeCnt ++;

			}else if(category == "Other"){
				   otherRules[otherCnt] = attriRules[attriDataLyrIdx];
				   otherCnt ++;
			}

		}

		double* CurrentCellAttributeSimilarities = new double[5];
		int* NumOfLayersInCategories = new int[5];

		double* ClimateSimilarities = new double[ClimateLyrCnt];
		double* GeologySimilarities = new double[GeogLyrCnt];
		double* TerrainSimilarities = new double[TerrainLyrCnt];
		double* VegetationSimilarities = new double[VegeLyrCnt];
		double* OthersSimilarities = new double[OtherLyrCnt];

		
//=================================================================================================================================
		for(int rowIdx = 0; rowIdx < Rows; rowIdx++){
			 for(int colIdx = 0; colIdx < Cols; colIdx++){
				double MaxSimilarity = 0;
				double SumSimilarity = 0;
				double SumSimilarityValue = 0;
				double MaxValue = -999999;
				for(int sampleIndex =0; sampleIndex< numOfSamples; sampleIndex ++){
					// climate
					for(int attriDataLyrIdx = 0; attriDataLyrIdx < ClimateLyrCnt; attriDataLyrIdx++){
					   double gridV = climate2Dvalues[attriDataLyrIdx][rowIdx * Cols + colIdx];
						if(abs(gridV - climateAllNoData[attriDataLyrIdx])< VERY_SMALL){
							ClimateSimilarities[attriDataLyrIdx] = -9999;
						}else{
								double sampleV = sampleClimateV[attriDataLyrIdx * numOfSamples + sampleIndex];
								double range = climateVRange[attriDataLyrIdx];
								ClimateSimilarities[attriDataLyrIdx] = climateRules[attriDataLyrIdx].getAttributeSimilarity(gridV,sampleV,range,climatesumAll[attriDataLyrIdx][sampleIndex],climatecounterAll[attriDataLyrIdx],climateStd[attriDataLyrIdx]);//,climatesum[attriDataLyrIdx][sampleIndex],climatecounter[attriDataLyrIdx],climateStd[attriDataLyrIdx]                                  
							}
					}
					double FinalClimateSimilarity = catIntegration.GetCategorySimilarity(ClimateSimilarities, ClimateLyrCnt);
					CurrentCellAttributeSimilarities[0] = FinalClimateSimilarity;
					NumOfLayersInCategories[0] = ClimateLyrCnt;
					// geology
					for(int attriDataLyrIdx = 0; attriDataLyrIdx < GeogLyrCnt; attriDataLyrIdx++){                      
					   double gridV = geo2Dvalues[attriDataLyrIdx][rowIdx * Cols + colIdx];
					   if(abs(gridV - geoAllNoData[attriDataLyrIdx])< VERY_SMALL){
							GeologySimilarities[attriDataLyrIdx] = -9999;
					   }else{

							double sampleV = sampleGeologyV[attriDataLyrIdx * numOfSamples + sampleIndex];
							double range = geologyVRange[attriDataLyrIdx];
							GeologySimilarities[attriDataLyrIdx] = geologyRules[attriDataLyrIdx].getAttributeSimilarity(gridV,sampleV,range,geologysumAll[attriDataLyrIdx][sampleIndex],geologycounterAll[attriDataLyrIdx],geoStd[attriDataLyrIdx]);//,//,geoStd[attriDataLyrIdx]
					   }

					}

					FinalClimateSimilarity = catIntegration.GetCategorySimilarity(GeologySimilarities, GeogLyrCnt);
					CurrentCellAttributeSimilarities[1] = FinalClimateSimilarity;
					NumOfLayersInCategories[1] = GeogLyrCnt;
					//terrain
					for(int attriDataLyrIdx = 0; attriDataLyrIdx < TerrainLyrCnt; attriDataLyrIdx++){
					   double gridV = terrain2Dvalues[attriDataLyrIdx][rowIdx * Cols + colIdx];
					   if(abs(gridV - terrainAllNoData[attriDataLyrIdx])< VERY_SMALL){
							TerrainSimilarities[attriDataLyrIdx] = -9999;							
					   }else{
							double sampleV = sampleTerrainV[attriDataLyrIdx * numOfSamples + sampleIndex];;
							double range = terrainVRange[attriDataLyrIdx];
							TerrainSimilarities[attriDataLyrIdx] = terrainRules[attriDataLyrIdx].getAttributeSimilarity(gridV,sampleV,range,terrainsumAll[attriDataLyrIdx][sampleIndex],terraincounterAll[attriDataLyrIdx],terrainStd[attriDataLyrIdx]);//,//,terrainStd[attriDataLyrIdx]
					   }


					}

					FinalClimateSimilarity = catIntegration.GetCategorySimilarity(TerrainSimilarities, TerrainLyrCnt);                    
					CurrentCellAttributeSimilarities[2] = FinalClimateSimilarity;
					NumOfLayersInCategories[2] = TerrainLyrCnt;
					//vegetation
					for(int attriDataLyrIdx = 0; attriDataLyrIdx < VegeLyrCnt; attriDataLyrIdx++){                 
					   double gridV = vege2Dvalues[attriDataLyrIdx][rowIdx * Cols + colIdx];
					   if(abs(gridV - vegetationAllNoData[attriDataLyrIdx])< VERY_SMALL){
							VegetationSimilarities[attriDataLyrIdx] = -9999;
					   }else{
							double sampleV = sampleVegetationV[attriDataLyrIdx * numOfSamples + sampleIndex];;
							double range = vegeVRange[attriDataLyrIdx];							
							VegetationSimilarities[attriDataLyrIdx] = vegetationRules[attriDataLyrIdx].getAttributeSimilarity(gridV,sampleV,range,vegetationsumAll[attriDataLyrIdx][sampleIndex],vegetationcounterAll[attriDataLyrIdx],vegetationStd[attriDataLyrIdx]);//,//,vegetationStd[attriDataLyrIdx]
					   }

					}

					FinalClimateSimilarity = catIntegration.GetCategorySimilarity(VegetationSimilarities, VegeLyrCnt);

					CurrentCellAttributeSimilarities[3] = FinalClimateSimilarity;
					NumOfLayersInCategories[3] = VegeLyrCnt;
					// others
					for(int attriDataLyrIdx = 0; attriDataLyrIdx < OtherLyrCnt; attriDataLyrIdx++){                        
					   double gridV = other2Dvalues[attriDataLyrIdx][rowIdx * Cols + colIdx];
					   if(abs(gridV - otherAllNoData[attriDataLyrIdx])< VERY_SMALL){
							OthersSimilarities[attriDataLyrIdx] = -9999;
					   }else{

							double sampleV = sampleOtherV[attriDataLyrIdx * numOfSamples + sampleIndex];;
							double range = otherVRange[attriDataLyrIdx];							
							OthersSimilarities[attriDataLyrIdx] = otherRules[attriDataLyrIdx].getAttributeSimilarity(gridV,sampleV,range,othersumAll[attriDataLyrIdx][sampleIndex],othercounterAll[attriDataLyrIdx],otherStd[attriDataLyrIdx]);//,/,otherStd[attriDataLyrIdx]
					   }

					}
					FinalClimateSimilarity = catIntegration.GetCategorySimilarity(OthersSimilarities, OtherLyrCnt);
					CurrentCellAttributeSimilarities[4] = FinalClimateSimilarity;
					NumOfLayersInCategories[4] = OtherLyrCnt;

					int count = 5;
					//单个样点与单个待推测点的相似度
					double SampleSimilarity = sampleIntegration.GetSampleSimilarity(CurrentCellAttributeSimilarities, count, NumOfLayersInCategories);                    			
					double PropertyValue = sampleValue[sampleIndex];
					if (SampleSimilarity >= MaxSimilarity)
					{
						MaxSimilarity = SampleSimilarity;
						MaxValue = PropertyValue;
					}
					SumSimilarity += SampleSimilarity;
					SumSimilarityValue += SampleSimilarity * PropertyValue;

				}// end of for(int sampleIndex =0; sampleIndex< numOfSamples; sampleIndex ++){

				double Uncertainty = 1 - MaxSimilarity;
				if ( Uncertainty < this->uncertaintyThreshold ){
					propertyVals[rowIdx][colIdx] = SumSimilarityValue/SumSimilarity;
					uncertaintyVals[rowIdx][colIdx] = Uncertainty;
					//if(propertyVals[rowIdx][colIdx]>40 || propertyVals[rowIdx][colIdx]<10){
					//   cout<<propertyVals[rowIdx][colIdx]<<endl;
					//}
				}else{

					propertyVals[rowIdx][colIdx] = noData;
					uncertaintyVals[rowIdx][colIdx] = noData;

				}		  

		//	climateSimipart[rowIdx][colIdx] = CurrentCellAttributeSimilarities[0];   //wtf
		//	terrainSimipart[rowIdx][colIdx] = CurrentCellAttributeSimilarities[2];   //wtf
		//	climate0Simipart[rowIdx][colIdx] = ClimateSimilarities[0];//wtf
			}//end of for
		}// end of for		
		delete[] climateRules;
		delete[] geologyRules;
		delete[] terrainRules;
		delete[] vegetationRules;
		delete[] otherRules;
		delete[] ClimateSimilarities;
		delete[] GeologySimilarities;
		delete[] TerrainSimilarities;
		delete[] VegetationSimilarities;
		delete[] OthersSimilarities;
		delete [] CurrentCellAttributeSimilarities;
		delete [] NumOfLayersInCategories;
   }

void ParallelSBMp::split(string s, string sep, vector<string> &flds){
	if (s.length() == 0)
		return;
	string fld;
	unsigned int i = 0, j = 0;
	do {
		j = s.find_first_of(sep,i);
		if (j > s.length())
			j = s.length();
		//fld = string(s,i,j-i);
		fld = s.substr(i,j-i);
		flds.push_back(fld);
		i = j + 1;
	}while(j < s.length());
}


