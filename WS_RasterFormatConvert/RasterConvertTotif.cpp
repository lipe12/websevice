#include <iostream>
#include <strstream>
#include <string>
#include <vector>
#include "gdal.h"
#include "gdal_priv.h"

using namespace std;

bool RasterConvert(char* inputFilename, char* outputFilename, char* outputType)
{
	GDALAllRegister();
	GDALDataset *pDatasetRead;
	GDALDataset *pDatasetSave;
	GDALDriver  *pDriver;
	pDriver = GetGDALDriverManager()->GetDriverByName(outputType);
	if(pDriver == NULL)
	{
		return false;
	}
	pDatasetRead = (GDALDataset*)GDALOpen(inputFilename, GA_ReadOnly);
	if(pDatasetRead == NULL)
	{
		return false;
	}
	pDatasetSave = pDriver->CreateCopy(outputFilename, pDatasetRead, FALSE, NULL, NULL, NULL);
	if(pDatasetSave == NULL)
	{
		return false;
	}
	GDALClose((GDALDatasetH)pDatasetRead);	// 关闭数据集
	GDALClose((GDALDatasetH)pDatasetSave);	// 关闭数据集
	return true;
}

int main(int argc, char *argv[])
{
	char *rasterFile = argv[1];
	char* outputFilename = argv[2];
	char* outputType = "GTiff";
	RasterConvert(rasterFile, outputFilename, outputType);
	return 0;
}