#ifndef PARALLELSBMP_H
#define PARALLELSBMP_H

#include <vector>
#include <mpi.h>
#include <omp.h>
#include "AscGrid.h"
#include "AttributeRule.h"
#include "CategoryIntegration.h"
#include "SampleIntegration.h"


#define VERY_SMALL 0.000001
using namespace std;

class ParallelSBMp
{
private:

	vector<AttributeRule> attriRules;
	CategoryIntegration catIntegration;
	SampleIntegration sampleIntegration;
	double uncertaintyThreshold;

	//
	double *sampleValue;

	double * sampleClimateV;
	double * sampleGeologyV;
	double * sampleTerrainV;
	double * sampleVegetationV;
	double * sampleOtherV;

	double * climateVRange;
	double * geologyVRange;
	double * terrainVRange;
	double * vegeVRange;
	double * otherVRange;

	double string_to_double( const std::string& s );

public:
	ParallelSBMp();
	ParallelSBMp(double * sampleValue);
	ParallelSBMp(double * sampleValue,double * sampleClimateV, double * sampleGeologyV, double * sampleTerrainV,
		double * sampleVegetationV, double * sampleOtherV,
		vector<string> &AttriRules,
		double * climateVRange,double * geologyVRange,double * terrainVRange,double * vegeVRange,double * otherVRange,
		string catIntegrationMethod, string sampleIntegrationMethod,
		string uncertaintyThreshold);

	void getPropertyMap(double *climateStd,double *geoStd,double *terrainStd,double *vegetationStd,double *otherStd,double ** climate2Dvalues, double ** geo2Dvalues,
		double ** terrain2Dvalues, double ** vege2Dvalues, double ** other2Dvalues,
		double ** propertyVals, double ** uncertaintyVals,
		int numOfSamples,
		int ClimateLyrCnt, int GeogLyrCnt, int TerrainLyrCnt, int VegeLyrCnt, int OtherLyrCnt,
		int Rows, int Cols,double noData, double *climateAllNoData,double *geoAllNoData,double *terrainAllNoData,
		double *vegetationAllNoData,double *otherAllNoData,
		double **climatesumAll, double **geologysumAll, double **terrainsumAll, double **vegetationsumAll, double **othersumAll,
		int *climatecounterAll, int *geologycounterAll, int *terraincounterAll, int *vegetationcounterAll, int *othercounterAll);//wtf
	void split(string s, string sep, vector<string> &flds);
};

#endif
