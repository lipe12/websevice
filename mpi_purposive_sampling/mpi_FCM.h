
#include <iostream>
#include <stdio.h>
#include <string>
#include <sstream>
#include <math.h>
#include <vector>
#include <mpi.h>
#include<time.h>
#include<gdal_priv.h>
using namespace std;
#define Eps 0.00001
class mpi_FCM
{
public:
	mpi_FCM(int number,int maxIteration,float error);
	~mpi_FCM();

public:
	int FCM_Processing(char* inputfile, char* outputfile, float error);
	void strParsing(char* inputfile);   //字符串解析
	void readInfo(const char* filename);  //解析信息
	void readDataFile(const char* filename);
	void writeDataFile(const char* filename,int tag);
	void createFile(const char *destfilename,int tag);
	void CreateRandom(int k,int n,int *centerIndex); //产生可以随输出控制的 k与n （可舍弃） 
	void InitDegree();//初始化隶属度
	void fnClusterMembershipValues();   //使每一个样本与各个聚类中心的隶属度之和为1
	double fnDistance(int i, int j); //计算距离
	void fnDegreeMembership(); //计算隶属度
	void fnClusterCentroid(); //计算聚类中心
	double fnObjectiveFunction(); //计算目标函数
public:
	int row; //行
	int column; //列
	int imageNum; //影像数目
	int block; //分块数据量
	int number; //分类数目
	int rank; //进程号
	int numProcs; //进程数目
	int interval;//°ŽÐÐ·ÖœâÊýŸÝ,ÐÐÊý
	int subRow; //分块行数
	int *subData;
	int maxIteration; //最大迭代次数
	double Fuzzyness;
	double tolerance;
	vector<char *>allPathfiles; //数据
	vector<char *>outPathfiles;
	string format;//ÎÄŒþžñÊœ
	string projection;
	double noData;
	double* pMyarray; //ž÷žöœø³ÌµÃµœµÄÊý×é
	vector <double *> vecAllData; //ÊýŸÝ
	double*center; //中心
	double *sumNumerator; //分子求和
	double *sumDenominator; //分母求和
	double *toalNumerator; //分子归并
	double *toalDenominator; //分母归并
	double* pTransform;
	int* pClass;
	double objectValue;
	double objectValueNew;
	double** degree;//隶属度
	double partitionCoef;
	double entropy;
	int samplingPixelNum;
};
