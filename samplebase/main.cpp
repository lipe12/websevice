#include <cstdlib>
#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <mpi.h>
//#include <omp.h>
#include "ParallelSBMp.h"
#include <vector>
#include "AscGrid.h"
#include "AttributeRule.h"
#include "CategoryIntegration.h"
#include "SampleIntegration.h"
#include <sstream>

#include <fstream>
#include <iomanip>
#define VERY_SMALL 0.000001
#define MEMORY_SIZE 5000

using namespace std;

double string_to_double( const std::string& s )
{
	std::istringstream i(s);
	double x;
	if (!(i >> x))
		return 0;
	return x;
}

static double string_to_double1( const std::string& s )
{
   std::istringstream i(s);
   double x;
   if (!(i >> x))
     return 0;
   return x;
 }

void parseStr(string str, char c, vector<string>& tokens){
	unsigned int posL = 0;
	unsigned int posR = 0;
	while(posR < str.length()-1){
		posR = str.find_first_of(c,posL);
		string sub = str.substr(posL,posR-posL);
		tokens.push_back(sub);
		posL = posR + 1;
	}
}
/*去除字符串中的换行符
*/
string& trimLineBreak(string &s)
{
	if (s.empty())
	{
		return s;
	}
	s.erase(0,s.find_first_not_of("\n"));
	s.erase(s.find_last_not_of("\n") + 1);
	return s;
}

// the sum value for gaussian
double **getSum(double ** category2Dvalues, int categoryLayerCount, double *categoryAllNoData, double * sampleCategoryV,
				int numOfSamples, int Rows, int Cols)//wtf
{
	double** categorySum;//the value of categorySum will be divided into four parts as for the parallel computing

	categorySum = new double*[categoryLayerCount];
	for(int attriDataLyrIdx = 0; attriDataLyrIdx < categoryLayerCount; attriDataLyrIdx++)
	  {
		  categorySum[attriDataLyrIdx]=new double[numOfSamples];
	  }
	double** categorySumAll;// 	it is the sum of categorySum
		categorySumAll = new double*[categoryLayerCount];
		for(int attriDataLyrIdx = 0; attriDataLyrIdx < categoryLayerCount; attriDataLyrIdx++)
		  {
			  categorySumAll[attriDataLyrIdx]=new double[numOfSamples];
		  }
	for(int attriDataLyrIdx = 0; attriDataLyrIdx < categoryLayerCount; attriDataLyrIdx++)
	{
	  for(int sampleIndex =0; sampleIndex< numOfSamples; sampleIndex ++)
	  {
		  categorySum[attriDataLyrIdx][sampleIndex]=0;
		  categorySumAll[attriDataLyrIdx][sampleIndex]=0;
		  double sampleV = sampleCategoryV[attriDataLyrIdx * numOfSamples + sampleIndex];
		   for(int rowIdx = 0; rowIdx < Rows; rowIdx++)
		   {
				for(int colIdx = 0;colIdx < Cols; colIdx++)
				{
					double gridV = category2Dvalues[attriDataLyrIdx][rowIdx * Cols + colIdx];
					if(abs(gridV - categoryAllNoData[attriDataLyrIdx])>= VERY_SMALL)
					{
						categorySum[attriDataLyrIdx][sampleIndex]+=pow(gridV-sampleV,2);
					}
				}
		  }
		  MPI_Allreduce(&categorySum[attriDataLyrIdx][sampleIndex],&categorySumAll[attriDataLyrIdx][sampleIndex],1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
		}
	}
	//MPI_Allreduce(categorySum,categorySumAll,categoryLayerCount*numOfSamples,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	return categorySumAll;
}
// the count value for gaussian
int *getCount(double ** category2Dvalues, int categoryLayerCount, double *categoryAllNoData, int Rows, int Cols)//wtf
{
	int *categoryCounter;
	categoryCounter=new int[categoryLayerCount];
	int *categoryCounterAll;
		categoryCounterAll=new int[categoryLayerCount];
	for(int attriDataLyrIdx = 0; attriDataLyrIdx < categoryLayerCount; attriDataLyrIdx++)
	{
	  categoryCounter[attriDataLyrIdx]=0;
	  categoryCounterAll[attriDataLyrIdx]=0;
	   for(int rowIdx = 0; rowIdx < Rows; rowIdx++)
	   {
			for(int colIdx = 0;colIdx < Cols; colIdx++)
			{
				double gridV = category2Dvalues[attriDataLyrIdx][rowIdx * Cols + colIdx];
				if(abs(gridV - categoryAllNoData[attriDataLyrIdx])>= VERY_SMALL)
				{
					categoryCounter[attriDataLyrIdx]+=1;
				}
			}
	  }
	   MPI_Allreduce(&categoryCounter[attriDataLyrIdx],&categoryCounterAll[attriDataLyrIdx],1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
	}
	//MPI_Allreduce(categoryCounter,categoryCounterAll,categoryLayerCount,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
	return categoryCounterAll;
}

int main(int argc, char *argv[]){

	GDALAllRegister();
	MPI::Init(argc, argv);

	double starttime, endtime, IOReadTime, IOWriteTime, readEnd, writBegin;
	starttime = MPI_Wtime();
	int rank, size;
	size = MPI::COMM_WORLD.Get_size();
	rank = MPI::COMM_WORLD.Get_rank();
    // prepare parameters begin
    char *  environLyrsPath =  argv[1];// enviroment layers
    char * samplePath = argv[2]; //样点文件
	char * attriRules = argv[3]; //Geology?Boolean#Climate?Gower#Terrain?Gower#Climate?Gower
    char * threshold = argv[4]; //阈值
    char * proLayrPath = argv[5]; //推测结果土壤属性图的存放路径
    char * unLayrPath = argv[6]; //推测结果土壤不确定性图的存放路径

    char * xName = argv[7]; //样点文件中坐标x列的名称
    char * yName = argv[8];
    char * propertyName = argv[9];//样点文件中要推测的属性列的名称
    char * categoryMethod = argv[10];//Limit或者Average
    char * sampleMethod = argv[11];

  //  char * cliLayrPath = argv[7];//wtf
  //  char * terLayrPath = argv[8];//wtf
  //  char * cli0LayrPath = argv[9];//wtf

	char* tablePath = samplePath;

	//const char* categoryIntegrationMethod = "Limit";
	//const char* sampleIntegrationMethod = "Limit";
	char* uncertaintyThreshold = threshold;
    char* propertyResultPath = proLayrPath;
	char* uncertaintyResultPath = unLayrPath;

	char* xCoor = xName;
	char* yCoor = yName;
	char* property = propertyName;
	char* categoryIntegrationMethod = categoryMethod;
	char* sampleIntegrationMethod = sampleMethod;

	// prepare parameter end
	int ClimateLyrCnt = 0;
	int GeogLyrCnt = 0;
	int TerrainLyrCnt = 0;
	int VegeLyrCnt =0;
	int OtherLyrCnt =0;
	int totalRows = 0;
	int totalCols = 0;
	int numOfSamples = 0;

	// analys the attributeRules
	vector<string> AttributeRules;
	//attriRules = argv[3];
	string attriRulesList = string(attriRules);
	parseStr(attriRulesList,'#',AttributeRules);   // Geology?Boolean#Climate?Gower#Terrain?Gower#Climate?Gower
	// count the numbers of enviromental types
	vector<string> type0fLayers;

	for(int i =0; i< AttributeRules.size(); i++){

	   vector<string> ruleParameterArray;
	   parseStr(AttributeRules[i],'?',ruleParameterArray);
	   string type = ruleParameterArray[0];
	   type0fLayers.push_back(type);

	   if(type == "Climate"){
	       ClimateLyrCnt ++;
	   }else if(type == "Geology"){
	       GeogLyrCnt ++;
	   }else if(type == "Terrain"){
	        TerrainLyrCnt ++;
	   }else if(type == "Vegetation"){
	         VegeLyrCnt ++;
	   }else if(type == "Other"){
	         OtherLyrCnt ++;
	   }

	}

	AscGrid * climateLyrs = NULL;
	AscGrid * geologyLyrs = NULL;
	AscGrid * terrainLyrs = NULL;
	AscGrid * vegetationLyrs = NULL;
	AscGrid * otherLyrs = NULL;

	double * climateVRange = NULL;
	double * geologyVRange = NULL;
	double * terrainVRange = NULL;
	double * vegeVRange = NULL;
	double * otherVRange = NULL;

	double * sampleValues = NULL;

	double * sampleClimateV = NULL;
	double * sampleGeologyV = NULL;
	double * sampleTerrainV = NULL;
	double * sampleVegetationV = NULL;
	double * sampleOtherV  = NULL;

	double noData = -9999;

	double lowerLeftY;
	double lowerLeftX;
	double cellSize;

	int climateCnt = 0;
	int geologyCnt = 0;
	int terrainCnt = 0;
	int vegeCnt = 0;
	int otherCnt = 0;

    double *climateAllNoData = NULL;
    double *geoAllNoData = NULL;
    double *terrainAllNoData = NULL;
    double *vegetationAllNoData = NULL;
    double *otherAllNoData = NULL;//wtf
	//=====================================================
	double *climateStd = NULL;//wtf
	double *geoStd = NULL;
	double *terrainStd = NULL;
	double *vegetationStd = NULL;
	double *otherStd = NULL;
	//=====================================================

	if (rank == 0){
		// load all enviroment layers
		vector<string> environLyrs;
		parseStr(string(environLyrsPath),'#',environLyrs);
		int numOflayers = environLyrs.size();
	    if(ClimateLyrCnt >0){
		   climateLyrs = new AscGrid[ClimateLyrCnt];
		   climateVRange = new double[ClimateLyrCnt];
		   climateAllNoData = new double[ClimateLyrCnt];//wtf
		   climateStd = new double[ClimateLyrCnt];//wtf
		}
		if(GeogLyrCnt > 0){
		   geologyLyrs = new AscGrid[GeogLyrCnt];
		   geologyVRange = new double[GeogLyrCnt];
		   geoAllNoData = new double[GeogLyrCnt];//wtf
		   geoStd = new double[GeogLyrCnt];//wtf
		}
		if(TerrainLyrCnt > 0){
		   terrainLyrs = new AscGrid[TerrainLyrCnt];
           terrainVRange = new double[TerrainLyrCnt];
           terrainAllNoData = new double[TerrainLyrCnt];//wtf
           terrainStd = new double[TerrainLyrCnt];//wtf
		}
		if(VegeLyrCnt > 0){
		   vegetationLyrs = new AscGrid[VegeLyrCnt];
		   vegeVRange = new double[VegeLyrCnt];
		   vegetationAllNoData = new double[VegeLyrCnt];//wtf
		   vegetationStd = new double[VegeLyrCnt];//wtf
		}if(OtherLyrCnt > 0){
		   otherLyrs = new AscGrid[OtherLyrCnt];
		   otherVRange = new double[OtherLyrCnt];
		   otherAllNoData = new double[OtherLyrCnt];//wtf
		   otherStd = new double[OtherLyrCnt];//wtf
		}

		// range = max - min, will be broadCast to other processes(jing cheng)

		for(int i = 0; i < numOflayers; i ++){
			AscGrid attributeLyr;
			attributeLyr.readAscGridGDAL(environLyrs[i]);
			string type = type0fLayers[i];

			double range = attributeLyr.getMax() - attributeLyr.getMin();

			//cout<<"stdvalue:"<<attributeLyr.getStdValue()<<endl;//wtf

		    if(type == "Climate"){
					climateLyrs[climateCnt] = attributeLyr;
					climateVRange[climateCnt] = range;
					climateStd[climateCnt]=attributeLyr.getStdValue();
					climateCnt ++;
		    }else if(type == "Geology"){
					geologyLyrs[geologyCnt] = attributeLyr;
					geologyVRange[geologyCnt] = range;
					geoStd[geologyCnt]=attributeLyr.getStdValue();
					geologyCnt ++;
		    }else if(type == "Terrain"){
				   terrainLyrs[terrainCnt] = attributeLyr;
					terrainVRange[terrainCnt] = range;
					terrainStd[terrainCnt]=attributeLyr.getStdValue();
					terrainCnt ++;
		    }else if(type == "Vegetation"){
					vegetationLyrs[vegeCnt] = attributeLyr;
					vegeVRange[vegeCnt] = range;
					vegetationStd[vegeCnt]=attributeLyr.getStdValue();
					vegeCnt ++;
		    }else if(type == "Other"){
					otherLyrs[otherCnt] = attributeLyr;
					otherVRange[otherCnt] = range;
					otherStd[otherCnt]=attributeLyr.getStdValue();
					otherCnt ++;
		    }

		}

        for(int i=0;i<ClimateLyrCnt;i++)
        {
            climateAllNoData[i] = climateLyrs[i].getNodaVal();
        }
        for(int i=0;i<GeogLyrCnt;i++)
        {
            geoAllNoData[i]= geologyLyrs[i].getNodaVal();
        }
        for(int i=0;i<TerrainLyrCnt;i++)
        {
            terrainAllNoData[i] = terrainLyrs[i].getNodaVal();
        }
        for(int i=0;i<VegeLyrCnt;i++)
        {
            vegetationAllNoData[i] = vegetationLyrs[i].getNodaVal();
        }
        for(int i=0;i<OtherLyrCnt;i++)
        {
            otherAllNoData[i] = otherLyrs[i].getNodaVal();
        }

	    // such as follow double * will be sended to other processes(jing cheng)


        // load sample file
        vector<vector<string> > sampleList;
		vector<string> fields;//存储的是已根据逗号分割好的每一行的数据

		FILE *pfin = NULL;
		pfin = fopen(tablePath,"r");
		if(pfin == NULL)
			cerr<<"fail to open sample file"<<endl;
		char row[500];
		while (fgets(row,500,pfin)!= NULL){
			string line = string(row);
			parseStr(line,',',fields);
			sampleList.push_back(fields);
			fields.clear();
		}
		fclose(pfin);
		// compute the sample point's value in  each raster
		if(ClimateLyrCnt>0){
		    totalRows = climateLyrs[0].getNumOfRows();
			totalCols = climateLyrs[0].getNumOfCols();
			lowerLeftY = climateLyrs[0].getYCor();
			lowerLeftX = climateLyrs[0].getXCor();
			cellSize = climateLyrs[0].getCellSize();
		}else if(GeogLyrCnt>0){
			totalRows = geologyLyrs[0].getNumOfRows();
			totalCols = geologyLyrs[0].getNumOfCols();
			lowerLeftY = geologyLyrs[0].getYCor();
			lowerLeftX = geologyLyrs[0].getXCor();
			cellSize = geologyLyrs[0].getCellSize();

		}else if(TerrainLyrCnt>0){
			totalRows = terrainLyrs[0].getNumOfRows();
			totalCols = terrainLyrs[0].getNumOfCols();
			lowerLeftY = terrainLyrs[0].getYCor();
			lowerLeftX = terrainLyrs[0].getXCor();
			cellSize = terrainLyrs[0].getCellSize();

		}else if(VegeLyrCnt>0){
			totalRows = vegetationLyrs[0].getNumOfRows();
			totalCols = vegetationLyrs[0].getNumOfCols();
			lowerLeftY = vegetationLyrs[0].getYCor();
			lowerLeftX = vegetationLyrs[0].getXCor();
			cellSize = vegetationLyrs[0].getCellSize();
		}else if(OtherLyrCnt>0){
		   totalRows = otherLyrs[0].getNumOfRows();
			totalCols = otherLyrs[0].getNumOfCols();
			lowerLeftY = otherLyrs[0].getYCor();
			lowerLeftX = otherLyrs[0].getXCor();
			cellSize = otherLyrs[0].getCellSize();
		}

		numOfSamples = sampleList.size() -1; // the firt row in sampefile is title
		int columnNums = sampleList[0].size();// 样点文件列数,sampleList[0]即样点文件中的第一行title
		int xIndex = 0, yIndex = 0, propertyIndex = 0;
	/*
	图像中以左上角为起点，起点向右为x轴，逐渐增大，起点向下为y轴，逐渐减小。本程序中（见AscGrid.cpp）xCor，yCor为左上角的坐标。
	其中yCor在AscGrid.cpp中
	readAscGridGDAL等相关函数中被转换为了左下角的y坐标（xCor = pTransform[0];   yCor = pTransform[3] - pTransform[1]*totalRows;// 由左上角变为左下角）
	这时候我们的起点就是笛卡尔坐标系中的原点了，向右为x轴，逐渐增大，向上为y轴，逐渐增大。又左上角和左下角的x坐标是一样的，
	因此，xCor,yCor就可以当做是左下角坐标了即lowerLeftX，lowerLeftY
	我们在将样点的地理坐标转为图像中的行列位置时就以左下角为起点参考了。
	*/

		int* sampleRows = new int[numOfSamples];
		int* sampleCols = new int[numOfSamples];
		for (int j=0; j<columnNums; j++)
        {
			//注意去除换行符
            if (strcmp(xCoor, trimLineBreak(sampleList[0][j]).c_str()) == 0)
            {
                xIndex = j;
            }

            if (strcmp(yCoor, trimLineBreak(sampleList[0][j]).c_str()) == 0)
            {
                yIndex = j;
            }

            if (strcmp(property, trimLineBreak(sampleList[0][j]).c_str()) == 0)
            {
                propertyIndex = j;
            }


        }
		for(int i = 1; i<= numOfSamples;i++){
			double curX = string_to_double(sampleList[i][xIndex]);//x属性列的值
			double curY = string_to_double(sampleList[i][yIndex]);
			sampleRows[i-1] = totalRows - (int)((curY - lowerLeftY)/cellSize)-1;
			sampleCols[i-1] = (int)((curX	- lowerLeftX)/cellSize);
		}
		//such  as follow double * will be broadCast to other processes(jing cheng)

		sampleValues = new double[numOfSamples];
		if(ClimateLyrCnt>0){
		   sampleClimateV = new double[numOfSamples * ClimateLyrCnt];
		}
		if(GeogLyrCnt>0){
		  sampleGeologyV = new double[numOfSamples * GeogLyrCnt];
		}
	    if(TerrainLyrCnt>0){
		  sampleTerrainV = new double[numOfSamples * TerrainLyrCnt];
		}
		if(VegeLyrCnt>0){
			sampleVegetationV = new double[numOfSamples * VegeLyrCnt];
		}
	    if(OtherLyrCnt>0){
			sampleOtherV = new double[numOfSamples * OtherLyrCnt];
		}

        for(int i =1; i<= numOfSamples; i++){
		    sampleValues[i -1] = string_to_double(sampleList[i][propertyIndex]);//推测属性列的值
		}

		for(int i =0; i< ClimateLyrCnt; i ++){
		   for(int j =0; j<numOfSamples; j++){
		       int rowIndex = sampleRows[j];
			   int colIndex = sampleCols[j];
		       sampleClimateV[i*numOfSamples + j] = climateLyrs[i].values[rowIndex][colIndex];
		   }
		}//end of for

		for(int i =0; i< GeogLyrCnt; i ++){
		   for(int j =0; j<numOfSamples; j++){
		       int rowIndex = sampleRows[j];
			   int colIndex = sampleCols[j];
		       sampleGeologyV[i*numOfSamples + j] = geologyLyrs[i].values[rowIndex][colIndex];
		   }
		}//end of for

		for(int i =0; i< TerrainLyrCnt; i ++){
		   for(int j =0; j<numOfSamples; j++){
		       int rowIndex = sampleRows[j];
			   int colIndex = sampleCols[j];
		       sampleTerrainV[i*numOfSamples + j] = terrainLyrs[i].values[rowIndex][colIndex];
		   }
		}//end of for

		for(int i =0; i< VegeLyrCnt; i ++){
		   for(int j =0; j<numOfSamples; j++){
		       int rowIndex = sampleRows[j];
			   int colIndex = sampleCols[j];
		       sampleVegetationV[i*numOfSamples + j] = vegetationLyrs[i].values[rowIndex][colIndex];
		   }
		}//end of for

		for(int i =0; i< OtherLyrCnt; i ++){
		   for(int j =0; j<numOfSamples; j++){
		       int rowIndex = sampleRows[j];
			   int colIndex = sampleCols[j];
		       sampleOtherV[i*numOfSamples + j] = otherLyrs[i].values[rowIndex][colIndex];
		   }
		}//end of for


	   delete[] sampleRows;
	   delete[] sampleCols;

	}

	//MPI 广播语句 将ClimateLyrCnt，GeogLyrCnt，TerrainLyrCnt，VegeLyrCnt，OtherLyrCnt，totalRows，totalCols广播到其他进程
	double computeTime1Begin = MPI_Wtime();
	MPI::COMM_WORLD.Bcast(&ClimateLyrCnt, 1, MPI_INT, 0);
	MPI::COMM_WORLD.Bcast(&GeogLyrCnt, 1, MPI_INT, 0);
	MPI::COMM_WORLD.Bcast(&TerrainLyrCnt, 1, MPI_INT, 0);
	MPI::COMM_WORLD.Bcast(&VegeLyrCnt, 1, MPI_INT, 0);
	MPI::COMM_WORLD.Bcast(&OtherLyrCnt, 1, MPI_INT, 0);
	MPI::COMM_WORLD.Bcast(&totalRows, 1, MPI_INT, 0);
	MPI::COMM_WORLD.Bcast(&totalCols, 1, MPI_INT, 0);
	MPI::COMM_WORLD.Bcast(&numOfSamples, 1, MPI_INT, 0);
//=================================================================================================================
	MPI::COMM_WORLD.Bcast(&noData, 1, MPI_DOUBLE, 0);//wtf
//================================================================================================================
	MPI::COMM_WORLD.Bcast(&lowerLeftY, 1, MPI_DOUBLE, 0);
	MPI::COMM_WORLD.Bcast(&lowerLeftX, 1, MPI_DOUBLE, 0);
	MPI::COMM_WORLD.Bcast(&cellSize, 1, MPI_DOUBLE, 0);
	double computeTime1End = MPI_Wtime();
	double computeTime1 = computeTime1End - computeTime1Begin;

	//其他进程重新分配内存

	int numOfsamClimate = numOfSamples * ClimateLyrCnt;
	int numOfsamGeology = numOfSamples * GeogLyrCnt;
	int numOfsamTerrain = numOfSamples * TerrainLyrCnt;
	int numOfsamVegetation = numOfSamples * VegeLyrCnt;
	int numOfsamOther = numOfSamples * OtherLyrCnt;
	if (rank != 0)
	{
		sampleValues = new double[numOfSamples];
		if(ClimateLyrCnt>0){
	      	sampleClimateV = new double[numOfsamClimate];
		}
		if(GeogLyrCnt>0){
		  sampleGeologyV = new double[numOfsamGeology];
		}
	    if(TerrainLyrCnt>0){
		  sampleTerrainV = new double[numOfsamTerrain];
		}
		if(VegeLyrCnt>0){
		  sampleVegetationV = new double[numOfsamVegetation];
		}
		if(OtherLyrCnt>0){
		  sampleOtherV = new double[numOfsamOther];
		}
	}

	// MPI 广播样点的值sampleValues, sampleClimateV, sampleGeologyV, sampleTerrainV, sampleVegetationV, sampleOtherV
	double computeTime2Begin = MPI_Wtime();
	MPI::COMM_WORLD.Bcast(sampleValues, numOfSamples, MPI_DOUBLE, 0);
	if(numOfsamClimate>0){
	    MPI::COMM_WORLD.Bcast(sampleClimateV, numOfsamClimate, MPI_DOUBLE, 0);
	}
    if(numOfsamGeology>0){
		MPI::COMM_WORLD.Bcast(sampleGeologyV, numOfsamGeology, MPI_DOUBLE, 0);
	}
    if(numOfsamTerrain>0){
	    MPI::COMM_WORLD.Bcast(sampleTerrainV, numOfsamTerrain, MPI_DOUBLE, 0);
	}
	if(numOfsamVegetation>0){
	   MPI::COMM_WORLD.Bcast(sampleVegetationV, numOfsamVegetation, MPI_DOUBLE, 0);
	}
	if(numOfsamOther>0){
	   MPI::COMM_WORLD.Bcast(sampleOtherV, numOfsamOther, MPI_DOUBLE, 0);
	}
	double computeTime2End = MPI_Wtime();
	double computeTime2 = computeTime2End - computeTime2Begin;
	//cout<<"Bcast sample value in raster"<<endl;

	// MPI 广播环境变量的范围 climateVRange, geologyVRange, terrainVRange, vegeVRange, otherVRange
	if(rank !=0){
	    if(ClimateLyrCnt>0){
		 climateVRange = new double[ClimateLyrCnt];
		 climateAllNoData = new double[ClimateLyrCnt];//wtf
		 climateStd = new double[ClimateLyrCnt];//wtf
		}
	    if(GeogLyrCnt>0){
		geologyVRange = new double[GeogLyrCnt];
		geoAllNoData = new double[GeogLyrCnt];//wtf
		geoStd = new double[GeogLyrCnt];//wtf
		}
		if(TerrainLyrCnt>0){
		  terrainVRange = new double[TerrainLyrCnt];
		  terrainAllNoData = new double[TerrainLyrCnt];//wtf
		  terrainStd = new double[TerrainLyrCnt];//wtf
		}
		if(VegeLyrCnt>0){
	      vegeVRange = new double[VegeLyrCnt];
	      vegetationAllNoData = new double[VegeLyrCnt];//wtf
	      vegetationStd = new double[VegeLyrCnt];//wtf
		}
		if(OtherLyrCnt>0){
		  otherVRange = new double[OtherLyrCnt];
		  otherAllNoData = new double[OtherLyrCnt];//wtf
		  otherStd = new double[OtherLyrCnt];//wtf
		}

	}

	double computeTime3Begin = MPI_Wtime();
	//=============================================================================
	if(ClimateLyrCnt>0){
	 MPI::COMM_WORLD.Bcast(climateVRange, ClimateLyrCnt, MPI_DOUBLE, 0);
	 MPI::COMM_WORLD.Bcast(climateAllNoData, ClimateLyrCnt, MPI_DOUBLE, 0);//wtf
	 MPI::COMM_WORLD.Bcast(climateStd, ClimateLyrCnt, MPI_DOUBLE, 0);//wtf
	}
	if(GeogLyrCnt>0){
	 MPI::COMM_WORLD.Bcast(geologyVRange, GeogLyrCnt, MPI_DOUBLE, 0);
	 MPI::COMM_WORLD.Bcast(geoAllNoData, GeogLyrCnt, MPI_DOUBLE, 0);//wtf
	 MPI::COMM_WORLD.Bcast(geoStd, GeogLyrCnt, MPI_DOUBLE, 0);//wtf
	}
	if(TerrainLyrCnt>0){
	  MPI::COMM_WORLD.Bcast(terrainVRange, TerrainLyrCnt, MPI_DOUBLE, 0);
	  MPI::COMM_WORLD.Bcast(terrainAllNoData, TerrainLyrCnt, MPI_DOUBLE, 0);//wtf
	  MPI::COMM_WORLD.Bcast(terrainStd, TerrainLyrCnt, MPI_DOUBLE, 0);//wtf
	}
	if(VegeLyrCnt){
	 MPI::COMM_WORLD.Bcast(vegeVRange, VegeLyrCnt, MPI_DOUBLE, 0);
	 MPI::COMM_WORLD.Bcast(vegetationAllNoData, VegeLyrCnt, MPI_DOUBLE, 0);//wtf
	 MPI::COMM_WORLD.Bcast(vegetationStd, VegeLyrCnt, MPI_DOUBLE, 0);//wtf
	}
	if(OtherLyrCnt>0){
	 MPI::COMM_WORLD.Bcast(otherVRange, OtherLyrCnt, MPI_DOUBLE, 0);
	 MPI::COMM_WORLD.Bcast(otherAllNoData, OtherLyrCnt, MPI_DOUBLE, 0);//wtf
	 MPI::COMM_WORLD.Bcast(otherStd, OtherLyrCnt, MPI_DOUBLE, 0);//wtf
	}
    //==============================================================================
    double computeTime3End = MPI_Wtime();
    double computeTime3 = computeTime3End - computeTime3Begin;

	int col_min = 0;
	int col_max = totalCols;
	int row_min = 0;
	int row_max = totalRows;
	//divide the whole extent
	int remainRows = (row_max - row_min) % size;
	int step = (row_max - row_min) / size;
	row_min = row_min + rank * step;
	row_max = row_min + step;
	if (rank == size - 1){
		row_max += remainRows;
	}

	int row_size = row_max - row_min;
	int col_size = col_max - col_min;

	double** climate2Dvalues = NULL;
	double** geo2Dvalues  = NULL;
	double** terrain2Dvalues = NULL;
    double** vege2Dvalues = NULL;
	double** other2Dvalues = NULL;

	if(ClimateLyrCnt>0)
		climate2Dvalues = new double * [ClimateLyrCnt];

	if(GeogLyrCnt>0)
		geo2Dvalues = new double * [GeogLyrCnt];

	if(TerrainLyrCnt>0)
		terrain2Dvalues = new double * [TerrainLyrCnt];

	if(VegeLyrCnt>0)
	   vege2Dvalues = new double * [VegeLyrCnt];

	if(OtherLyrCnt>0)
	   other2Dvalues = new double * [OtherLyrCnt];

	climateCnt = 0;
    geologyCnt = 0;
    terrainCnt = 0;
    vegeCnt = 0;
    otherCnt = 0;

	// load all enviroment layers
	//environLyrsPath = argv[1];
	vector<string> environLyrs;
	parseStr(string(environLyrsPath),'#',environLyrs);

	for(int i =0; i<type0fLayers.size(); i++){
	   string type = type0fLayers[i];
	   AscGrid attributeLyr;
	   double * rawData = attributeLyr.readAscGridGDAL_Block(environLyrs[i],0,row_min,row_size,col_size);
	   if(type == "Climate"){
	        climate2Dvalues[climateCnt] = rawData;
	        climateCnt ++;
	   }else if(type == "Geology"){
	        geo2Dvalues[geologyCnt] = rawData;

			geologyCnt ++;
	   }else if(type == "Terrain"){
	        terrain2Dvalues[terrainCnt] = rawData;

			terrainCnt ++;
	   }else if(type == "Vegetation"){
	        vege2Dvalues[vegeCnt] = rawData;

			vegeCnt ++;
	   }else if(type == "Other"){
	        other2Dvalues[otherCnt] = rawData;

			otherCnt ++;
	   }

	}
	// calculate the sum and count which is essential for Gaussian

	double** climatesumAll;
	climatesumAll=new double*[ClimateLyrCnt];
	for(int i = 0; i < ClimateLyrCnt; i++)
	  {
		  climatesumAll[i]=new double[numOfSamples];
	  }
	int *climatecounterAll;
	climatecounterAll=new int[ClimateLyrCnt];

	double** geologysumAll ;//单个图层单个样点有一个sum值
	geologysumAll=new double*[GeogLyrCnt];
	for(int i = 0; i < GeogLyrCnt; i++)
	  {
		  geologysumAll[i]=new double[numOfSamples];
	  }
	int *geologycounterAll ;//单个图层有一个counter值
	geologycounterAll=new int[GeogLyrCnt];

	double** terrainsumAll ;//单个图层单个样点有一个sum值
	terrainsumAll=new double*[TerrainLyrCnt];
	for(int i = 0; i < TerrainLyrCnt; i++)
	  {
		  terrainsumAll[i]=new double[numOfSamples];
	  }
	int *terraincounterAll;
	terraincounterAll = new int[TerrainLyrCnt];

	double** vegetationsumAll ;//单个图层单个样点有一个sum值
	vegetationsumAll=new double*[VegeLyrCnt];
	for(int i = 0; i < VegeLyrCnt; i++)
	  {
		  vegetationsumAll[i]=new double[numOfSamples];
	  }
	int *vegetationcounterAll ;//单个图层有一个counter值
	vegetationcounterAll=new int[VegeLyrCnt];

	double** othersumAll ;//单个图层单个样点有一个sum值
	othersumAll=new double*[OtherLyrCnt];
	for(int i = 0; i < OtherLyrCnt; i++)
	  {
		  othersumAll[i]=new double[numOfSamples];
	  }
	int *othercounterAll ;//单个图层有一个counter值
	othercounterAll=new int[OtherLyrCnt];

	climatesumAll = getSum(climate2Dvalues, ClimateLyrCnt, climateAllNoData, sampleClimateV ,numOfSamples, row_size, col_size);
	geologysumAll = getSum(geo2Dvalues, GeogLyrCnt, geoAllNoData, sampleGeologyV,numOfSamples, row_size, col_size);
	terrainsumAll = getSum(terrain2Dvalues, TerrainLyrCnt, terrainAllNoData, sampleTerrainV,numOfSamples, row_size, col_size);

	vegetationsumAll = getSum(vege2Dvalues, VegeLyrCnt, vegetationAllNoData, sampleVegetationV,numOfSamples, row_size, col_size);
	othersumAll = getSum(other2Dvalues, OtherLyrCnt, otherAllNoData, sampleOtherV,numOfSamples, row_size, col_size);

	climatecounterAll = getCount(climate2Dvalues, ClimateLyrCnt, climateAllNoData, row_size, col_size);
	geologycounterAll = getCount(geo2Dvalues, GeogLyrCnt, geoAllNoData, row_size, col_size);
	terraincounterAll = getCount(terrain2Dvalues, TerrainLyrCnt, terrainAllNoData, row_size, col_size);
	vegetationcounterAll = getCount(vege2Dvalues, VegeLyrCnt, vegetationAllNoData, row_size, col_size);
	othercounterAll = getCount(other2Dvalues, OtherLyrCnt, otherAllNoData, row_size, col_size);

	// calculate the property and uncertainty
	double **propertyValues;
	double **uncertaintyValues;
	//double **climateSimipart;//wtf
	//double **terrainSimipart;
	//double **climate0Simipart;

	double* propertyValuesStorage = new double[row_size*col_size];
	double* uncertaintyValuesStorage = new double[row_size*col_size];
	//double* climateSimipartStorage = new double[row_size*col_size];//wtf
	//double* terrainSimipartStorage = new double[row_size*col_size];
	//double* climate0SimipartStorage = new double[row_size*col_size];
	propertyValues = new double*[row_size];
	//climateSimipart = new double*[row_size];
	//climate0Simipart = new double*[row_size];
	//terrainSimipart = new double*[row_size];//wtf
	uncertaintyValues = new double*[row_size];
	for(int i = 0; i < row_size; i ++){     //wtf
		propertyValues[i] = &propertyValuesStorage[i*col_size];
		uncertaintyValues[i] = &uncertaintyValuesStorage[i*col_size];
		//climateSimipart[i] = &climateSimipartStorage[i*col_size];
		//climate0Simipart[i] = &climate0SimipartStorage[i*col_size];
		//terrainSimipart[i] = &terrainSimipartStorage[i*col_size];
	}
	readEnd = MPI_Wtime();
	IOReadTime = readEnd - starttime;

    double computeTime4Begin =  MPI_Wtime();
	ParallelSBMp inf(sampleValues,sampleClimateV, sampleGeologyV, sampleTerrainV,
		sampleVegetationV, sampleOtherV,
		AttributeRules,
		climateVRange,geologyVRange,terrainVRange,vegeVRange,otherVRange,
		categoryIntegrationMethod, sampleIntegrationMethod,
		uncertaintyThreshold
		);

	inf.getPropertyMap(climateStd,geoStd,terrainStd,vegetationStd,otherStd,climate2Dvalues, geo2Dvalues,
		terrain2Dvalues, vege2Dvalues, other2Dvalues,propertyValues, uncertaintyValues,	numOfSamples,
		ClimateLyrCnt, GeogLyrCnt, TerrainLyrCnt, VegeLyrCnt, OtherLyrCnt,	row_size, col_size,noData,
		climateAllNoData,geoAllNoData,terrainAllNoData,vegetationAllNoData,otherAllNoData,
		climatesumAll, geologysumAll, terrainsumAll, vegetationsumAll, othersumAll,
		climatecounterAll, geologycounterAll, terraincounterAll, vegetationcounterAll, othercounterAll);//wtf

     double computeTime4End =  MPI_Wtime();
     double computeTime4 = 	 computeTime4End - computeTime4Begin;
	 //writBegin = MPI_Wtime();
	//integration
	double ** propertyAll;
	double* propertyAllStorage;

	//double **climateAll;//wtf
	//double* climateAllStorage;//wtf

	//double **climate0All;//wtf
	//double* climate0AllStorage;//wtf

	//double **terrainAll;//wtf
	//double* terrainAllStorage;//wtf

	double ** uncertaintyAll;
	double* uncertaintyAllStorage;

	if (rank == size - 1){
		propertyAllStorage = new double[totalRows*col_size];
		//climateAllStorage = new double[totalRows*col_size];//wtf
		//climate0AllStorage = new double[totalRows*col_size];//wtf
		//terrainAllStorage = new double[totalRows*col_size];//wtf
		uncertaintyAllStorage = new double[totalRows*col_size];
		propertyAll = new double*[totalRows];
		//climateAll = new double*[totalRows];//wtf
		//climate0All = new double*[totalRows];//wtf
		//terrainAll = new double*[totalRows];//wtf
		uncertaintyAll = new double*[totalRows];
		for (int i = 0; i <totalRows; i ++){
			propertyAll[i] = &propertyAllStorage[i*col_size];
			//climateAll[i] = &climateAllStorage[i*col_size];//wtf
			//climate0All[i] = &climate0AllStorage[i*col_size];//wtf
			//terrainAll[i] = &terrainAllStorage[i*col_size];//wtf
			uncertaintyAll[i] = &uncertaintyAllStorage[i*col_size];
		}
	}

	int displ=0;
	int recvCnt = row_size * col_size;
	int *recvCntRoot = new int[size];
	int *displRoot = new int[size];

	double computeTime5Begin = MPI_Wtime();

	MPI::COMM_WORLD.Gather(&recvCnt,1,MPI::INT,recvCntRoot,1,MPI::INT,size - 1);
	MPI::COMM_WORLD.Bcast(recvCntRoot,size,MPI::INT,size - 1);

	for (int k = 0; k < rank; k ++){
		displ += recvCntRoot[k];
	}
	MPI::COMM_WORLD.Gather(&displ,1,MPI::INT,displRoot,1,MPI::INT,size - 1);
	MPI::COMM_WORLD.Bcast(displRoot,size,MPI::INT,size - 1);
	MPI::COMM_WORLD.Gatherv(propertyValuesStorage,row_size*col_size,
		MPI_DOUBLE,propertyAllStorage,recvCntRoot,displRoot,MPI_DOUBLE,size-1);

	MPI::COMM_WORLD.Gatherv(uncertaintyValuesStorage,row_size*col_size,
		MPI_DOUBLE,uncertaintyAllStorage,recvCntRoot,displRoot,MPI_DOUBLE,size-1);

	//MPI::COMM_WORLD.Gatherv(climateSimipartStorage,row_size*col_size,
	//	MPI_DOUBLE,climateAllStorage,recvCntRoot,displRoot,MPI_DOUBLE,size-1);//wtf

//	MPI::COMM_WORLD.Gatherv(climate0SimipartStorage,row_size*col_size,
	//	MPI_DOUBLE,climate0AllStorage,recvCntRoot,displRoot,MPI_DOUBLE,size-1);//wtf

	//MPI::COMM_WORLD.Gatherv(terrainSimipartStorage,row_size*col_size,
		//MPI_DOUBLE,terrainAllStorage,recvCntRoot,displRoot,MPI_DOUBLE,size-1);//wtf

	double computeTime5End = MPI_Wtime();
	double computeTime5 = computeTime5End - computeTime5Begin;
	writBegin = MPI_Wtime();
	if (rank == size - 1){

		AscGrid propertyMap(totalCols,totalRows,lowerLeftX,lowerLeftY,cellSize,noData,propertyAll);//wtf no_data
		//AscGrid climateMap(totalCols,totalRows,lowerLeftX,lowerLeftY,cellSize,noData,climateAll);//wtf
		//AscGrid climate0Map(totalCols,totalRows,lowerLeftX,lowerLeftY,cellSize,noData,climate0All);//wtf
		//AscGrid terrainMap(totalCols,totalRows,lowerLeftX,lowerLeftY,cellSize,noData,terrainAll);//wtf
		AscGrid uncertaintyMap(totalCols,totalRows,lowerLeftX,lowerLeftY,cellSize,noData,uncertaintyAll);//wtf
        //	string climateFile = climateResultPath;//wtf
	//	string climate0File = climate0ResultPath;//wtf
		//string terrainFile = terrainResultPath;//wtf
		string propertyFile = propertyResultPath;
		string uncertaintyFile = uncertaintyResultPath;

		//climateMap.createAscGridGADL(environLyrs[0],climateFile);//wtf
		//climateMap.writeAscGridGDAL(climateFile); //wtf

		//climate0Map.createAscGridGADL(environLyrs[0],climate0File);//wtf
		//climate0Map.writeAscGridGDAL(climate0File); //wtf


		//terrainMap.createAscGridGADL(environLyrs[0],terrainFile);//wtf
		//terrainMap.writeAscGridGDAL(terrainFile); //wtf

		propertyMap.createAscGridGADL(environLyrs[0],propertyFile);
		propertyMap.writeAscGridGDAL(propertyFile);
		uncertaintyMap.createAscGridGADL(environLyrs[0],uncertaintyFile);
		uncertaintyMap.writeAscGridGDAL(uncertaintyFile);

	}
	delete[] propertyValues;
	delete[] uncertaintyValues;
	delete[] propertyValuesStorage;
	delete[] uncertaintyValuesStorage;
	//delete[] climateSimipart;//wtf
	//delete[] climateSimipartStorage;//wtf
	//delete[] climate0Simipart;//wtf
//	delete[] climate0SimipartStorage;//wtf
	//delete[] terrainSimipart;//wtf
	//delete[] terrainSimipartStorage;//wtf
	delete[] climateAllNoData;
	delete[] geoAllNoData;
	delete[] terrainAllNoData;
	delete[] vegetationAllNoData;
	delete[] otherAllNoData;
	delete[] climateStd;
	delete[] geoStd;
	delete[] terrainStd;
	delete[] vegetationStd;
	delete[] otherStd;

	if(rank == size -1){
		delete[] propertyAll;
		delete[] uncertaintyAll;
		delete[] propertyAllStorage;
		delete[] uncertaintyAllStorage;
	//	delete[] climateAll;//wtf
	//	delete[] climateAllStorage;//wtf
	//	delete[] climate0All;//wtf
	//	delete[] climate0AllStorage;//wtf
	//	delete[] terrainAll;//wtf
	//	delete[] terrainAllStorage;//wtf

        // delete []climateAllNoData;//wtf
        // delete []geoAllNoData;
        // delete []terrainAllNoData;
        // delete []vegetationAllNoData;
        // delete []otherAllNoData;

        // delete []climateStd;
        // delete []geoStd;
        // delete []terrainStd;
        // delete []vegetationStd;
        // delete []otherStd;
	}
	delete[] recvCntRoot;
	delete[] displRoot;

	// release memory

	if(climateLyrs != NULL)
		delete[] climateLyrs;
	if(geologyLyrs != NULL)
		delete[] geologyLyrs;
	if(terrainLyrs != NULL)
		delete[] terrainLyrs;
	if(vegetationLyrs != NULL)
		delete[] vegetationLyrs;
	if(otherLyrs != NULL)
		delete[] otherLyrs;

	if(climateVRange != NULL)
		delete[] climateVRange;
	if(geologyVRange != NULL)
		delete[] geologyVRange;
	if(terrainVRange != NULL)
		delete[] terrainVRange;
	if(vegeVRange != NULL)
		delete[] vegeVRange;
	if(otherVRange != NULL)
		delete[] otherVRange;

	if(sampleValues != NULL)
		delete[] sampleValues;

	if(sampleClimateV != NULL)
		delete[] sampleClimateV;
	if(sampleGeologyV != NULL)
		delete[] sampleGeologyV;
	if(sampleTerrainV != NULL)
		delete[] sampleTerrainV;
	if(sampleVegetationV != NULL)
		delete[] sampleVegetationV;
	if(sampleOtherV != NULL)
		delete[] sampleOtherV;

	endtime = MPI_Wtime();
	cout<<"hello world"<<endl;
	IOWriteTime = endtime - writBegin;
	if(rank == size - 1)
	{
		//cout<<"[DEBUG][TIMESPAN][IO]"<< IOReadTime+IOWriteTime  << endl;
		//cout<<"[DEBUG][TIMESPAN][COMPUTING]"<< writBegin - readEnd << endl;
		//cout<<"[DEBUG][TIMESPAN][Comunicate1]"<< computeTime1  << endl;
		//cout<<"[DEBUG][TIMESPAN][Comunicate2]"<< computeTime2 << endl;
		//cout<<"[DEBUG][TIMESPAN][Comunicate3]"<< computeTime3<< endl;
		//cout<<"[DEBUG][TIMESPAN][Comunicate4]"<< computeTime5 << endl;
		cout<<"[DEBUG][TIMESPAN][COMPUTING]"<< computeTime1 + computeTime2 + computeTime3 +computeTime4 +computeTime5 << endl;
		cout<<"[DEBUG][TIMESPAN][IO]"<< endtime-starttime - computeTime1 - computeTime2 - computeTime3 -computeTime4 - computeTime5<< endl;
	}

	MPI::Finalize();
	return 0;
}

