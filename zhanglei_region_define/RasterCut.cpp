/*
* 功能： 将多个栅格数据进行叠加求得公共部分，输出公共部分，并输出公共部分的外边界矢量SHP文件
* 说明：
* 输入参数：4个
*	第1个：输入的栅格数据文件名，用'#'分割；
*	第2个：输出的栅格数据文件名，用'#'分割；
*	第3个：输出的临时栅格数据文件名，该临时文件最终会自动被删除；
*	第4个：输出的外边界SHP数据文件名。
*
* 函数返回值：
*	返回 0：  表示程序正确运行
*	返回 444：表示文件读取错误（输入或输出的文件名有误等情况）
*	返回 999：表示输入的参数格式不正确
*	其他：    表示存在图层之间没有交集的情况，返回值为与之前图层没有交集的图层序号（例如：第三个输入图层与前两个图层无交集，则输出第三个图层的序号:2）
*/

#include <stdio.h>
#include <iostream>
#include <vector>
#include <string>
#include "gdal_priv.h"
#include "ogrsf_frmts.h" //for ogr
#include "gdal_alg.h"	 //for GDALPolygonize
#include "cpl_conv.h"	 //for CPLMalloc()

using namespace std;

// 解析输入参数
void parseStr(string str, char c, vector<string>& tokens) {
	unsigned int posL = 0;
	unsigned int posR = 0;
	while(posR < str.length()-1) {
		posR = str.find_first_of(c, posL);
		string sub = str.substr(posL, posR-posL);
		tokens.push_back(sub);
		posL = posR + 1;
	}
}

// 地理坐标转行列号
bool Projection2ImageRowCol(double *adfGeoTransform, double dProjX, double dProjY, int &iCol, int &iRow)  
{  
	try  
	{  
		double dTemp = adfGeoTransform[1]*adfGeoTransform[5] - adfGeoTransform[2]*adfGeoTransform[4];  
		double dCol = 0.0, dRow = 0.0;  
		dCol = (adfGeoTransform[5]*(dProjX - adfGeoTransform[0]) -   
			adfGeoTransform[2]*(dProjY - adfGeoTransform[3])) / dTemp + 0.5;  
		dRow = (adfGeoTransform[1]*(dProjY - adfGeoTransform[3]) -   
			adfGeoTransform[4]*(dProjX - adfGeoTransform[0])) / dTemp + 0.5;  

		iCol = static_cast<int>(dCol);  
		iRow = static_cast<int>(dRow);  
		return true;  
	}  
	catch(...)  
	{  
		return false;  
	}  
}

// DEM数据公共部分裁剪
int DEMCut(vector<string> inputRasterFileNames,vector<string> outputRasterFileNames)
{
	GDALAllRegister();

	// 读入数据
	vector<GDALDataset *> datasets;
	for(int i = 0; i < inputRasterFileNames.size(); i++)
	{
		GDALDataset *ds;
		ds = (GDALDataset*)GDALOpen(inputRasterFileNames[i].c_str(), GA_ReadOnly);
		if(ds != NULL)
		{
			datasets.push_back(ds);
		}
		else
		{
			return 444;		// 文件读取有误，返回444
		}
	}

	// 获取公共部分的范围 left top right bottom
	// adfGeoTransform[6]  数组adfGeoTransform保存的是仿射变换中的一些参数，分别含义见下  
	// adfGeoTransform[0]  左上角x坐标   
	// adfGeoTransform[1]  东西方向分辨率  
	// adfGeoTransform[2]  旋转角度, 0表示图像 "北方朝上"  
	// adfGeoTransform[3]  左上角y坐标   
	// adfGeoTransform[4]  旋转角度, 0表示图像 "北方朝上"  
	// adfGeoTransform[5]  南北方向分辨率  
	double left = 0;
	double right = 0;
	double top = 0;
	double bottom = 0;
	for(int i = 0; i < datasets.size(); i++)
	{
		GDALDataset *ds;
		ds = datasets[i];
		double geoTransform[6];
		ds->GetGeoTransform(geoTransform);
		if(i == 0)
		{
			left = geoTransform[0];
			top = geoTransform[3];
			right = geoTransform[0] + geoTransform[1] * ds->GetRasterXSize() + geoTransform[2] * ds->GetRasterYSize();
			bottom = geoTransform[3] + geoTransform[4] * ds->GetRasterXSize() + geoTransform[5] * ds->GetRasterYSize();
		}
		else
		{
			double left_temp = geoTransform[0];
			double top_temp = geoTransform[3];
			double right_temp = geoTransform[0] + geoTransform[1] * ds->GetRasterXSize() + geoTransform[2] * ds->GetRasterYSize();
			double bottom_temp = geoTransform[3] + geoTransform[4] * ds->GetRasterXSize() + geoTransform[5] * ds->GetRasterYSize();

			// 判断是否有交集
			if(left_temp >= right || top_temp <= bottom || right_temp <= left || bottom_temp >= top)
			{
				return i;	// 结束该函数，第i个图层与之前的图层没有交集
			}

			// 有交集，继续执行
			if(left_temp > left)	// 取left最大值
			{
				left = left_temp;
			}
			if(right_temp < right)	// 取right最小值
			{
				right = right_temp;
			}
			if(top_temp < top)		// 取top最小值
			{
				top = top_temp;
			}
			if(bottom_temp > bottom)// 取bottom最大值
			{
				bottom = bottom_temp;
			}
		}
	}

	// 根据边界获取每个栅格文件需要裁剪出的行列范围，并输出裁剪后的栅格数据
	vector<float*> pDataList;
	int pixelCount = 0;
	float noDataValue = datasets[0]->GetRasterBand(1)->GetNoDataValue();
	for(int i = 0; i < datasets.size(); i++)
	{
		GDALDataset *ds;
		ds = datasets[i];
		double geoTransform[6];
		ds->GetGeoTransform(geoTransform);
		int row_start, col_start, row_end, col_end, rowCount, colCount;
		Projection2ImageRowCol(geoTransform, left, top, col_start, row_start);	// 根据地理坐标转换得到起点左上角的行列号
		Projection2ImageRowCol(geoTransform, right, bottom, col_end, row_end);	// 根据地理坐标转换得到终点右下角的行列号
		rowCount = row_end - row_start;// + 1;
		colCount = col_end - col_start;// + 1;
		pixelCount = colCount*rowCount;

		// 读取数据
		float *pData = (float*)CPLMalloc(sizeof(float)*colCount*rowCount);
		GDALDataType dataType = ds->GetRasterBand(1)->GetRasterDataType();
		ds->RasterIO(GF_Read, col_start, row_start, colCount, rowCount, pData, colCount, rowCount, GDT_Float32, 1, 0, 0, 0, 0);

		pDataList.push_back(pData);
	}
	float *pBoundaryData = (float*)CPLMalloc(sizeof(float)*pixelCount);
	for(int i = 0; i < pixelCount; i++)
	{
		bool flag = true;
		for(int j = 0; j < pDataList.size(); j++)
		{
			float value = pDataList[j][i];
			if(value == noDataValue)
			{
				flag = false;
				break;
			}
		}
		if(flag == true)	// 有数据
		{
			pBoundaryData[i] = 1;
		}
		else	// 没有数据,此栅格值为noDataValue
		{
			pBoundaryData[i] = 0;
		}
	}
	
	for(int i = 0; i < datasets.size(); i++)
	{
		GDALDataset *ds;
		ds = datasets[i];
		double geoTransform[6];
		ds->GetGeoTransform(geoTransform);
		int row_start, col_start, row_end, col_end, rowCount, colCount;
		Projection2ImageRowCol(geoTransform, left, top, col_start, row_start);	// 根据地理坐标转换得到起点左上角的行列号
		Projection2ImageRowCol(geoTransform, right, bottom, col_end, row_end);	// 根据地理坐标转换得到终点右下角的行列号
		rowCount = row_end - row_start;// + 1;
		colCount = col_end - col_start;// + 1;

		// 创建输出数据
		GDALDataset *output;
		GDALDriver  *pDriver;
		pDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
		if(pDriver == NULL)
		{
			cout<<"GDAL DriverManager Error!\n";
			return 444;
		}
		output = pDriver->Create(outputRasterFileNames[i].c_str(), colCount, rowCount, 1, GDT_Float32, NULL);
		if(output == NULL)
		{
			cout<<"GDAL Create Error!\n";
			return 444;
		}

		for(int j = 0; j < pixelCount; j++)
		{
			if(pBoundaryData[j] == 0)
			{
				pDataList.at(i)[j] = noDataValue;
			}
		}

		double newGeotransform[6];
		newGeotransform[0] = left;
		newGeotransform[1] = geoTransform[1];
		newGeotransform[2] = geoTransform[2];
		newGeotransform[3] = top;
		newGeotransform[4] = geoTransform[4];
		newGeotransform[5] = geoTransform[5];
		output->SetGeoTransform(newGeotransform);
		output->GetRasterBand(1)->RasterIO(GF_Write, 0, 0, colCount, rowCount, pDataList.at(i), colCount, rowCount, GDT_Float32, 0, 0);
		output->SetProjection(ds->GetProjectionRef());
		output->GetRasterBand(1)->SetNoDataValue(ds->GetRasterBand(1)->GetNoDataValue());

		GDALClose((GDALDatasetH)ds);		// 关闭数据集
		GDALClose((GDALDatasetH)output);	// 关闭数据集
	}
	return 0;
}

// 构造外轮廓的栅格文件
int CreateBoundaryRaster(vector<string> inputRasterFileNames, char *outputFileName)
{
	// 读入数据
	vector<GDALDataset *> datasets;
	for(int i = 0; i < inputRasterFileNames.size(); i++)
	{
		GDALDataset *ds;
		ds = (GDALDataset*)GDALOpen(inputRasterFileNames[i].c_str(), GA_ReadOnly);
		if(ds != NULL)
		{
			datasets.push_back(ds);
		}
		else
		{
			cout<<"文件读取错误！\n";
			return 444;
		}
	}
	int colCount = datasets[0]->GetRasterXSize();
	int rowCount = datasets[0]->GetRasterYSize();
	GDALDataset *output;
	GDALDriver  *pDriver;
	pDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
	if(pDriver == NULL)
	{
		cout<<"GDAL DriverManager Error!\n";
		return 444;
	}
	output = pDriver->Create(outputFileName, colCount, rowCount, 1, GDT_Float32, NULL);
	if(output == NULL)
	{
		cout<<"GDAL Create Error!\n";
		return 444;
	}
	
	float *pFinalData = (float*)CPLMalloc(sizeof(float)*colCount*rowCount);	// 存放最终输出栅格的数据
	vector<float *> pDataList;	// 存放每个栅格图层的数据
	for(int i = 0; i < datasets.size(); i++)	// 存储到pDataList中
	{
		GDALDataset *ds;
		ds = datasets[i];
		float *pData = (float*)CPLMalloc(sizeof(float)*colCount*rowCount);
		ds->RasterIO(GF_Read, 0, 0, colCount, rowCount, pData, colCount, rowCount, GDT_Float32, 1, 0, 0, 0, 0);
		pDataList.push_back(pData);
	}

	float noDataValue = datasets[0]->GetRasterBand(1)->GetNoDataValue();
	for(int i = 0; i < colCount*rowCount; i++)
	{
		bool flag = true;	// 判断是否有数据
		for(int j = 0; j < pDataList.size(); j++)
		{
			float value = pDataList[j][i];
			if(value == noDataValue)
			{
				flag = false;
				break;
			}
		}
		if(flag == true)	// 有数据
		{
			pFinalData[i] = 1;
		}
		else	// 没有数据,此栅格值为noDataValue
		{
			pFinalData[i] = 0;
		}
	}
	
	GDALDataset *ds = datasets[0];
	output->GetRasterBand(1)->RasterIO(GF_Write, 0, 0, colCount, rowCount, pFinalData, colCount, rowCount, GDT_Float32, 0, 0);
	double geoTransform[6];
	ds->GetGeoTransform(geoTransform);
	output->SetGeoTransform(geoTransform);
	output->SetProjection(ds->GetProjectionRef());
	output->GetRasterBand(1)->SetNoDataValue(0);

	for(int i = 0; i < datasets.size(); i++)	// 存储到pDataList中
	{
		GDALDataset *ds = datasets[i];
		GDALClose((GDALDatasetH)ds);		// 关闭数据集
	}
	
	GDALClose((GDALDatasetH)output);	// 关闭数据集
	return 0;
}

// 将外轮廓栅格文件矢量化
int ImagePolygonize(char *inputFileName, char* outputFileName, const char* pszFormat)
{
	GDALAllRegister();
	OGRRegisterAll();	// 添加驱动注册
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");

	GDALDataset* poSrcDS=(GDALDataset*)GDALOpen(inputFileName, GA_ReadOnly);
	if(poSrcDS==NULL)
	{
		return 444;
	}
	// 创建输出矢量文件
	OGRSFDriver *poDriver;
	poDriver = (OGRSFDriver*)OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName( pszFormat );
	if (poDriver == NULL)
	{  
		GDALClose((GDALDatasetH)poSrcDS); 
		return 444;
	}
	poDriver->DeleteDataSource(outputFileName);	// 若该文件存在，则先删除已有的shp文件
	//根据文件名创建输出矢量文件
	OGRDataSource* poDstDS = poDriver->CreateDataSource(outputFileName, NULL);
	if (poDstDS==NULL)
	{
		GDALClose((GDALDatasetH)poSrcDS);
		return 444;
	}
	// 定义空间参考，与输入图像相同;
	OGRSpatialReference *poSpatialRef = new OGRSpatialReference(poSrcDS->GetProjectionRef());
	OGRLayer* poLayer = poDstDS->CreateLayer("boundary", poSpatialRef, wkbPolygon, NULL);
	if (poDstDS == NULL)
	{
		GDALClose((GDALDatasetH)poSrcDS); 
		OGRDataSource::DestroyDataSource(poDstDS); 
		delete poSpatialRef; 
		poSpatialRef = NULL; 
		return 444;
	}
	OGRFieldDefn ofieldDef("Segment", OFTInteger);	//创建属性表，只有一个字段即“Segment”，里面保存对应的栅格的像元值
	poLayer->CreateField(&ofieldDef);
	GDALRasterBandH hSrcBand = (GDALRasterBandH) poSrcDS->GetRasterBand(1);		//获取图像的第一个波段
	GDALPolygonize(hSrcBand, NULL, (OGRLayerH)poLayer, 0, NULL, NULL, NULL);	//调用栅格矢量化

	// 删除noDataValue部分的要素
	OGRFeature *poFeature;
	poLayer->ResetReading();
	while( (poFeature = poLayer->GetNextFeature()) != NULL )	// 获得要素
	{
		int value = poFeature->GetFieldAsInteger("Segment");
		if(value == 0)	// 需要删除的noData要素
		{
			long fid = poFeature->GetFID();
			poLayer->DeleteFeature(fid);
		}
	}
	OGRFeatureDefn *pFeatureDefn = poLayer->GetLayerDefn();
	std::string strLayerName = pFeatureDefn->GetName();		// 读取该图层的名称
	std::string strSQL = "REPACK " + strLayerName;
	poDstDS->ExecuteSQL(strSQL.c_str(), NULL, "");
	
	GDALClose(poSrcDS); // 关闭文件
	OGRDataSource::DestroyDataSource(poDstDS);
	return 0;
}

// 删除临时文件
void DeleteTempFile(char *filename)
{
	remove(filename);
}

// 主要的处理函数
int DEMProcessing(vector<string> inputRasterFileNames, vector<string> outputRasterFileNames,
				   char *tempRasterFileName, char *shpFileName)
{
	int result = DEMCut(inputRasterFileNames, outputRasterFileNames);		// DEM数据公共部分裁剪
	if(result != 0)
	{
		return result;		// 结束操作，存在图层之间没有交集的情况，返回图层序号
	}
	CreateBoundaryRaster(outputRasterFileNames, tempRasterFileName);		// 构造外轮廓的栅格文件
	ImagePolygonize(tempRasterFileName, shpFileName, "ESRI Shapefile");		// 将外轮廓栅格文件矢量化
	DeleteTempFile(tempRasterFileName);										// 删除临时文件
	return 0;
}


int main(int argc, char *argv[])
{
	// ******************* 数据准备--载入测试数据 ******************* //

	// 参数示例
	//argc = 5;
	//argv[1] = "D:/data/demcut/dem1.tif#D:/data/demcut/dem2.tif";
	//argv[2] = "D:/RasterCut/output/dem1.tif#D:/RasterCut/output/dem2.tif";
	//argv[3] = "D:/RasterCut/output/tempRaster.tif";
	//argv[4] = "D:/RasterCut/output/boundary.shp";
	//char *argument = "D:/data/demcut/twi.tif#D:/data/demcut/slope.tif#D:/data/demcut/plan.tif#D:/data/demcut/Elevation.asc D:/RasterCut/output/twi.tif#D:/RasterCut/output/slope.tif#D:/RasterCut/output/plan.tif#D:/RasterCut/output/Elevation.tif D:/RasterCut/output/tempRaster.tif D:/RasterCut/output/boundary.shp";
	//char *argument = "D:/data/demcut/dem1.tif#D:/data/demcut/dem2.tif D:/RasterCut/output/dem1.tif#D:/RasterCut/output/dem2.tif D:/RasterCut/output/tempRaster.tif D:/RasterCut/output/boundary.shp";

	if(argc-1 != 4) // 判断输入参数个数是否为4
	{
		cout<<"您输入的参数有误！\n";
		return 999;		// 输入参数有误，返回999
	}

	char *str_inputRasterFileNames = argv[1];
	char *str_outputRasterFileNames = argv[2];
	char *str_tempRasterFileName = argv[3];
	char *str_shpFileName = argv[4];

	vector<string> inputRasterFileNames;
	parseStr(str_inputRasterFileNames, '#', inputRasterFileNames);
	vector<string> outputRasterFileNames;
	parseStr(str_outputRasterFileNames, '#', outputRasterFileNames);
	char *tempRasterFileName = str_tempRasterFileName;
	char *shpFileName = str_shpFileName;

	// ******************* 数据准备--完毕 ******************* //

	// 主要的处理函数
	int result = DEMProcessing(inputRasterFileNames, outputRasterFileNames, tempRasterFileName, shpFileName);
	if(result != 0)
	{
		if(result == 444)
		{
			cout<<"文件读取错误！\n";
		}
		else
		{
			cout<<"您输入的第"<<result+1<<"个图层与之前的图层无重叠部分，请输入正确的图层数据。\n";	// 存在图层之间没有交集的情况，返回与之前图层没有交集的图层序号
		}
		return result;
	}

	cout<<"\n--------DONE!--------\n";
	//system("pause");
	return 0;				// 运行正常没有错误，返回值为0
}